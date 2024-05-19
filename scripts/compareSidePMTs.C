// THIS SCRIPT'S IS BASED IN BECCA'S compareConfigurations.C
// and it's repurposed to analize the PMT on the Sides Geometry.

#include "simEvent.h"
#include <vector>
#include <map>
#include "TFile.h"
#include "TChain.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TGraph.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TLine.h"
#include "THStack.h"
#include "TGraphErrors.h"

// Vectors here are fancy "arrays"
std::vector<TChain*> vfin;
std::vector<TH2*> vhist; //an array of histograms
std::vector<TH2*> vhistD; //an array of histograms
std::vector<TH2*> vhistCUT; //an array of histograms
std::vector<TH2*> vhistO; //an array of histograms
std::vector<TH2*> vhistPos; //*New: array of histograms
std::vector<TH1*> vhistScalars; //an array of histograms
std::vector<TH1*> vhistDProj; //an array of projections of histograms
std::vector<TH1*> vhistOProj; //an array of projections of histograms
std::vector<TH1*> vhistPrecision;


//std::vector<std::string> fnames = {"data/test2/Sim_D2ODetector*.root"};
std::vector<std::string> fnames = {"data/test7_no_events_in_PMTs/Sim_D2ODetector*.root"};


//std::vector<std::string> fnames = {"data/Tank_WaterTC_config/ROOT_Plots/10K_Events_30_MeV_New_Geometry_12_PMTs/Sim_D2ODetector*.root","data/Tank_PMTs_Sides_config/30_MeV_electrons/Sim_D2ODetector*.root"}; //For 30 MeV events only
//std::vector<std::string> fnames = {"data/Tank_WaterTC_config/Sim_D2ODetector*.root","data/Tank_PMTs_Sides_config/0-55_MeV_electrons/Sim_D2ODetector*.root"};

//std::vector<const char*> config = {"D2O_PMT_Top_30_MeV"};
std::vector<const char*> config = {"PMTs_on_Sides"};


//std::vector<const char*> config = {"D2O_PMT_Top_30MeV","D2O_PMT_Sides_30MeV"};
//std::vector<const char*> config = {"D2O_PMT_Top_","D2O_PMT_Sides_"};

std::vector<double> d2oScale = {0.687};
// std::vector<double> d2oScale = {0.574,0.574};
//std::vector<double> d2oScale = {0.687,0.687};

double xsntotal_dcc = 5.5E-41;
double xsntotal_Occ = 8.6696988e-42;


double cutEnergy = 30;
//double cutEnergy = 50;
// double allowance = 0.01; //fraction of cut Energy allowed

//New from Becca's PlotFiducial():
int firstCutValue = 1;
int maxCutValue = 55;
int cutSteps = 1;
int cutEnergyF = firstCutValue; //Modified: cutEnergy->cutEnergyF for PlotFiducial
double allowance = 0.01; //fraction of cut Energy allowed
//

int color = 1;
   std::map<double, double> e0Sigma; //map for Standard Deviation
   std::map<double, double> hitsPerMeV; //map for Hits per MeV
//   std::vector<std::map<double, double>> e0Sigma; //New from Becca's PlotFiducial()
//   std::vector<std::map<double, double>> hitsPerMeV; //New from Becca's PlotFiducial()

double SNSHoursPerYear = 5000.0;
double chargePerPulse = 21.59E-6;
double POTperhour = 1.4/1.3*chargePerPulse/1.602E-19*60*3600;
//double POTperday = 7E20;
double POTperPulse = 1.4/1.3*chargePerPulse/1.602E-19;
double POTperYear = POTperhour*SNSHoursPerYear;
double pionsperPOT = 0.09; // @1010MeV/p
double enuPerYear = POTperYear*pionsperPOT;
double areaAt20m = 4*TMath::Pi()*TMath::Power(20.0*100.0,2.0);
double enuFluxat20m = enuPerYear/areaAt20m*2.0/*SNS-years*/;
// double enuFluxat20m = enuPerYear/areaAt20m*4.0/*SNS-years*/;
double deuteronsIn1ton = 1E6/(4.0+16.0)*TMath::Na()*2.0;
double oxygenIn1ton = 1E6/(4.0+16.0)*TMath::Na()*1.0;
double totalDCCPerTon = deuteronsIn1ton*xsntotal_dcc*enuFluxat20m;
double totalOCCPerTon = oxygenIn1ton*xsntotal_Occ*enuFluxat20m;
// double totalOCCPerTon = deuteronsIn1ton*xsntotal_Occ*enuFluxat20m;

bool logBinning = false;

//generating bins;
double* genLogBins(int nbins, double axismin, double axismax)
{
  double* dbins = new double[nbins+1];
  for(int ibin = 0; ibin <= nbins; ibin++)
  {
    dbins[ibin] = axismin*pow(axismax/axismin, ((double)ibin)/nbins);
  }
  return dbins;
}


void CleanUp(){
  for (int i = 0; i < vfin.size(); i++){
    delete vfin[i];
  }
  vfin.clear();
  
  for (int i = 0; i < vhist.size(); i++){
    delete vhist[i];
  }
  vhist.clear();

  for (int i = 0; i < vhistD.size(); i++){
    delete vhistD[i];
  }
  vhistD.clear();

  for (int i = 0; i < vhistCUT.size(); i++){
    delete vhistCUT[i];
  }
  vhistCUT.clear();

  for (int i = 0; i < vhistO.size(); i++){
    delete vhistO[i];
  }
  vhistO.clear();

  for (int i = 0; i < vhistScalars.size(); i++){
    delete vhistScalars[i];
  }
  vhistScalars.clear();

  for (int i = 0; i < vhistDProj.size(); i++){
    delete vhistDProj[i];
  }
  vhistDProj.clear();

  for (int i = 0; i < vhistOProj.size(); i++){
    delete vhistOProj[i];
  }
  vhistOProj.clear();

  for (int i = 0; i < vhistPrecision.size(); i++){
    delete vhistPrecision[i];
  }
  vhistPrecision.clear();

  return;
}

void compareConfig(){ //Reading in the data...
  int ival = 0;
  
  for(auto fname : fnames){ //For all listed configurations do all of the following:
    auto Sim_Tree = new TChain("Sim_Tree","Sim_Tree"); //("Default name of TTree","Title")
    std::cout << fname << ":\t Adding a " << config[ival] << " config\n"; //Prints the name of the configuration we are working with
    Sim_Tree->Add(fname.c_str()); // fname.c_str() returns pointer to fname: this is, it returns the name of the configuration. "Add" adds a branch to the tree (??)
    // Sim_Tree->Print();

    simEvent *ev = new simEvent(); //A new object (ev) from class simEvent. "ev" represents an event
    TVector3 pos = ev -> position0;//*New
    Sim_Tree -> SetBranchAddress("eventData", &ev); //Adding "eventData" branch to Sim_Tree (????). &ev is the address of the object "ev"
    
    //Defining the HIstograms. NOT filling them yet:
    
    TString hNameLY(Form("hLYPE_%sConfig",config[ival]));
    TH2* hLY = 0; //Defining a new pointer
    //hLY = new TH2D(hNameLY.Data(),hNameLY.Data(), 50, 0, 50, 200, 0, 2000); //A 2D histogram object named hLY 
    // Name="hNameLY.Data()", Title="hNameLY.Data()", X dims: 50 bins, from 0 to 50. Y dims:200 bins from 0 to 2000.
    hLY = new TH2D(hNameLY.Data(),hNameLY.Data(), 55, 0, 55, 200, 0, 2000);
    vhist.push_back(hLY); // This inserts the object hLY into the vector vhist on its last position and the size increases by 1.
    
    TString hNameD(Form("hLYPED_%sConfig",config[ival]));
    TH2* hLYD = 0;
    //hLYD = new TH2D(hNameD.Data(),hNameD.Data(),50,0,50,200,0,2000);
    hLYD = new TH2D(hNameD.Data(),hNameD.Data(),55,0,55,200,0,2000);
    vhistD.push_back(hLYD);

    TString hNameCUT(Form("hLYPECUT_%sConfig",config[ival]));
    TH2* hLYCUT = 0;
    //hLYCUT = new TH2D(hNameCUT.Data(),hNameCUT.Data(), 50,0,50,100,0,2000);
    hLYCUT = new TH2D(hNameCUT.Data(),hNameCUT.Data(), 55,0,55,100,0,2000);
    vhistCUT.push_back(hLYCUT);

    TString hNameO(Form("hLYPEO_%sConfig",config[ival]));
    TH2* hLYO = 0;
    //hLYO = new TH2D(hNameO.Data(),hNameO.Data(),50,0,50,200,0,2000);
    hLYO = new TH2D(hNameO.Data(),hNameO.Data(),55,0,55,200,0,2000);
    vhistO.push_back(hLYO);
    
     // New*
    TString hNamePos(Form("hLYPos_%sConfig",config[ival]));
     TH2* hLYPos = 0;
     hLYPos = new TH2D(hNamePos.Data(),hNamePos.Data(),100,-1000,1000,200,0,2000);
     vhistPos.push_back(hLYPos);

    //
     
    // New*: Plotting initial positions
    TString hNameInitPos(Form("hLYInitPos_%sConfig",config[ival]));
     TH2* hLYInitPos = 0;
     hLYInitPos = new TH2D(hNameInitPos.Data(),hNameInitPos.Data(),2000,-1000,1000,1000,-500,500);
     vhistPos.push_back(hLYInitPos);

    //
     

    TString hNameScalar(Form("hScalars_%sConfig",config[ival]));
    TH1* hScalar = new TH1D(hNameScalar.Data(),hNameScalar.Data(),10,0,10.0); //10 bins of size 1
    vhistScalars.push_back(hScalar);
    
    hLYD->SetDirectory(gROOT); // These objects are owned by ROOT (not sure what this means ???)
    hLYO->SetDirectory(gROOT);
    hLYPos->SetDirectory(gROOT);
    hLYInitPos->SetDirectory(gROOT);
    hScalar->SetDirectory(gROOT);
    
    int ievent = 0; //Filling out these 5 histograms!
    while(ievent < Sim_Tree -> GetEntries()){ //Sim_Tree -> GetEntries() is the number of generated events 
      Sim_Tree -> GetEntry(ievent); //Read entry (event) "ievent" from Sim_Tree
      int ndetpe = 0;
      hScalar->Fill(ev->vol0); //put (accumulate) ev->vol0 in hScalar histogram (1 dimensional).
      //vol0=1 if the event was in d2o volume, vol0=2 if it was in h2o volume.
      
     // Filling out the CUT histogram:
      if (ev -> sourceParticleEnergy > (cutEnergy - allowance*cutEnergy) && ev -> sourceParticleEnergy < (cutEnergy + allowance*cutEnergy)){ // If the energy of the event is between 30 and the allowance, then put this event into this 2D histogram.
	hLYCUT -> Fill(ev -> sourceParticleEnergy, ev -> numHits); //put [ev -> sourceParticleEnergy] in X, [ev -> numHits] in Y
      }
//*****************************CHECK THISSSSSSS!!! : 
    //
      if(ev->vol0==1){//if event happened in d2o volume...
        hLY->Fill(ev->sourceParticleEnergy,ev->numHits); //... put it in this histogram
        hLYD->Fill(ev->sourceParticleEnergy,ev->numHits);// and in this one too.
        // hLYD->Fill(ev->sourceParticleEnergy,ev->numHits,ev->WeightD()); //same but WEIGHTED by the nu-D cross section
      }
      
      if(ev->vol0>0){// because oxygen is in BOTH d2o and h2o. Recall that vol0=2 if event happened in h2o volume.
        hLYO->Fill(ev->sourceParticleEnergy,ev->numHits);
        // hLYO->Fill(ev->sourceParticleEnergy,ev->numHits,ev->WeightO()); //same but WEIGHTED by the nu-O crss section
        hLYPos->Fill(ev->position0.Z(),ev->numHits); //*New
        hLYInitPos->Fill(ev->position0.X(),ev->position0.Y()); //*New
      }
      
      ievent++; //number of events counter
    }
    
    ival++;  //number of configurations counter  
  } //end for
  
  //Putting all of this into a File: (Can we change the name of this file???? What's the RECREATE for?)
  
  TFile* fOut = TFile::Open("TCStudyHalfV2.root","RECREATE");
  for(auto h : vhist){
    h->SetDirectory(fOut);
  }
  for(auto h : vhistD){
    h->SetDirectory(fOut);
  }
  for(auto h : vhistCUT){
    h->SetDirectory(fOut);
  }
  for(auto h : vhistO){
    h->SetDirectory(fOut);
  }
  for(auto h : vhistPos){ //*New
    h->SetDirectory(fOut);
  }
  for(auto h : vhistScalars){
    h->SetDirectory(fOut);
  }

  fOut->Write();
  CleanUp();
  fOut->Close();
}// end of compareConfig()


/*
void ObtainResolution(){ //Energy Resolution stuff...
  for (int ival = 0; ival < config.size(); ival++){ //config.size() returns the number of configs being analized here!
    //Determine how to convert PE to MeV
    // cLightEReco -> cd();
    vhistD[ival]->FitSlicesY(); //Projects slices of vhistD along Y (all bins in x) AND FITS THEM WITH A GAUSSIAN (by default). See https://root.cern.ch/doc/master/classTH2.html#ace3d272713ed83badf620b2942e2e896 for more. 
    TH1* dMean = (TH1*)gDirectory->Get(Form("%s_1",vhistD[ival]->GetName())); //
    dMean -> SetXTitle("Source Energy (MeV)");
    dMean -> SetYTitle("Number of Hits in Event");
    // dMean -> SetLineColor(ival + 1);
    // dMean -> Draw("same");
    TFitResultPtr peToMeV = dMean -> Fit("pol1", "S0"); //RECONSTRUCTED ENERGY????

    vhistCUT[ival] -> ProjectionY(); //Takes Y projection of CUT histogram
    TH1* hRes = (TH1*)gDirectory->Get(Form("%s_py", vhistCUT[ival]->GetName()));
    TFitResultPtr bestFit = hRes -> Fit("gaus", "S0", "", 0.85*cutEnergy*peToMeV -> Parameter(1), 1.15*cutEnergy*peToMeV -> Parameter(1)); //FITTING WITH A GAUSSIAN?????

    hitsPerMeV[cutEnergy] = bestFit -> Parameter(1); //Adds this cut Energy to the map
    e0Sigma[cutEnergy] = bestFit -> Parameter(2); //Adds this cut Energy to the map
    
  }
} //end of ObtainResolution
*/



void LoadHists(){ //Loading histograms produced by compareConfig():
  TFile* fIn = 0;
  fIn = TFile::Open("TCStudyHalfV2.root"); //Preparing the file to get stuff from it
  
  int ival = 0;
  for(auto tc : config){ //do this for each configuration
    TString hNamePE(Form("hLYPE_%sConfig",config[ival]));
    auto hlype = (TH2*)fIn->Get(hNamePE.Data());
    //hlype->Scale(1.0,"width"); //multiply histogram by 1 in width??
    vhist.push_back(hlype);
    
    TString hNameD(Form("hLYPED_%sConfig",config[ival]));
    auto hlyped = (TH2*)fIn->Get(hNameD.Data());
    //hlyped->Scale(1.0,"width");
    vhistD.push_back(hlyped);

    TString hNameCUT(Form("hLYPECUT_%sConfig",config[ival]));
    auto hlypeCUT = (TH2*)fIn->Get(hNameCUT.Data());
    //hlypeCUT->Scale(1.0,"width");
    std::cout << hlypeCUT -> GetEntries()<< "\t" << hlypeCUT -> Integral() << "\n\n"; //for instance: "152(entries of the histogram hlypeCUT)	7.6 (its integral)"
    vhistCUT.push_back(hlypeCUT);
    
    TString hNameO(Form("hLYPEO_%sConfig",config[ival]));
    auto hlypeo = (TH2*)fIn->Get(hNameO.Data());
    //hlypeo->Scale(1.0,"width");
    vhistO.push_back(hlypeo);
    
    //*New:
    TString hNamePos(Form("hLYPos_%sConfig",config[ival]));
    auto hlypos = (TH2*)fIn->Get(hNamePos.Data());
    //hlypos->Scale(1.0,"width");
    vhistPos.push_back(hlypos);
    //*
    
    //*New:
    TString hNameInitPos(Form("hLYPos_%sConfig",config[ival]));
    auto hlyinitpos = (TH2*)fIn->Get(hNameInitPos.Data());
    //hlyinitpos->Scale(1.0,"width");
    vhistPos.push_back(hlyinitpos);
    //*
    
    TString hNameScalar(Form("hScalars_%sConfig",config[ival]));
    auto hscalars = (TH1*)fIn->Get(hNameScalar.Data());
    vhistScalars.push_back(hscalars);
    
    ival++;
  }//end for


} // end of LoadHists()

void PlotResolution(){
    
    double lowfit;
    double upfit;
    
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //Creating new canvases (Name,Title):
  TCanvas* cLightYield = new TCanvas("cLightYield","cLightYield");
  TCanvas* cSingle;
  TCanvas* cEnergyResolution = new TCanvas("cERes","cERes");

  TF1 *gaussian; //New 1D function named "gaussian"
  TLegend *tl = new TLegend(0.7, 0.6, 0.9, 0.9); //Creating legend box (x1,y1,x2,y2) item named "tl"
  TLegend *tl2 = new TLegend(0.7, 0.6, 0.9, 0.9);//Creating legend box item named "tl2"
  
//----- Producing the data for Energy Resolution Study:------//
  
  for (int ival = 0; ival < config.size(); ival++){ //for each configuration
    //-----Determine how to convert PE to MeV:-----//
    // cLightEReco -> cd();
    vhistD[ival]->FitSlicesY(); //Projects slices of vhistD along Y (all bins in x) AND FITS THEM WITH A GAUSSIAN (by default).
    TH1* dMean = (TH1*)gDirectory->Get(Form("%s_1",vhistD[ival]->GetName())); //A histogram named "dMean" to be filled with values of a TF1 function; in this case the projections of vhistD obtained by FitSlicesY(), which are denoted by vhistD_1.
    // Form(): Formats a string using a printf style format descriptor
    dMean -> SetXTitle("Source Energy (MeV)");
    dMean -> SetYTitle("Number of Hits in Event");
    // dMean -> SetLineColor(ival + 1);
    // dMean -> Draw("same");
    TFitResultPtr peToMeV = dMean -> Fit("pol1", "S0"); //Fits dMean to a Polynomial of grade 1 (pol1), then Saves (S) the result of the fit of dMean in an object of the type TFitResultPtr (and name peToMeV) without plotting it (0).
    
    // pol1 = new TF1(Form("fit_%s", vhistD[ival] -> GetName()), "pol1", 0, 60);
    // pol1 -> SetParameters(peToMeV -> Parameter(0), peToMeV -> Parameter(1));
    // pol1 -> SetLineColor(ival + 1);
    // pol1 -> Draw("same");

    cSingle = new TCanvas(Form("cLightYield_%s", config[ival]), Form("cLightYield_%s", config[ival])); //A new canvas (name, title) for Lighth Yield. 
    vhistCUT[ival] -> ProjectionY(); //Projects vhistCUT into Y axis to get light yield.
    TH1* hRes = (TH1*)gDirectory->Get(Form("%s_py", vhistCUT[ival]->GetName())); //A Histogram named "hRes" to store the projection of vhistCUT (vhistCUT_py)
    
    // Normalize by number of events (since just comparing widths)
    // hRes -> Scale(1./hRes -> GetMaximum());
    // hRes -> Scale(1./hRes -> Integral());
    // hRes -> Scale(1./hRes -> GetEntries());
    
    hRes -> SetLineColor(kRed); //Light Yield will be RED
    hRes -> SetLineWidth(2);
    // hRes -> GetYaxis() -> SetRangeUser(-0.001, 1.1);
    
    //Now let's fit our Light Yield to a nice Gaussian:
    //Fits will look different for both Water and Acrylic.

    //First, the fit ranges:
   /* if(ival == 0){//Acrylic
        
        double lowfit = 0.85;
        
    }
    
    if(ival == 1){//Water: we need a dual gaussian fit (*pending)
        double lowfit = 0.85;
    }*/

     TFitResultPtr bestFit = hRes -> Fit("gaus", "S0", "", 0.85*cutEnergy*peToMeV -> Parameter(1), 1.15*cutEnergy*peToMeV -> Parameter(1)); //Performing the gaussian fit: the result will be stored in "bestFit", and the fit will be applied from 0.85*cutEnergy*peToMeV, to 1.15*cutEnergy*peToMeV.
     
     gaussian = new TF1(Form("fit_%s", hRes -> GetName()), "gaus", 0.65*cutEnergy*peToMeV -> Parameter(1), 2000); //creating a new Gaussian Function called "gaussian". (Name, Title, xmin, xmax). Gaussian=(constant term, mean, sigma)
     
    /*TFitResultPtr bestFit = hRes -> Fit("gaus", "S0", "", lowfit*cutEnergy*peToMeV -> Parameter(1), 1.15*cutEnergy*peToMeV -> Parameter(1));
     gaussian = new TF1(Form("fit_%s", hRes -> GetName()), "gaus", lowfit*cutEnergy*peToMeV -> Parameter(1), 2000); */
     
     //Setting the parameters of our function gaussian:
     gaussian -> SetParameters(bestFit -> Parameter(0), bestFit -> Parameter(1), bestFit -> Parameter(2)); 
     gaussian -> SetLineWidth(2); 
     gaussian -> SetLineColor(kBlack);
 
     //What the hell is Parameter(0)???? 
     hitsPerMeV[cutEnergy] = bestFit -> Parameter(1); //Adds this cut Energy to the map
     e0Sigma[cutEnergy] = bestFit -> Parameter(2); //Adds this cut Energy to the map

    cSingle -> cd(); //Let's go to cSingle..
    hRes -> GetXaxis() -> SetRangeUser(0, 1200); //Our x-range (number of photoelectrons)
    hRes -> SetXTitle("Num. of Photoelectrons");
    hRes -> SetYTitle("Num. of Events");
    hRes -> Draw(); //Draw it in our canvas!
    gaussian -> Draw("same"); //Draw "gaussian" in the same plot as Light Yield.

  //   TF1* fGausEmg2 = new TF1("fGausEmg2","[3]/([1]*sqrt(2*pi))*exp(-0.5*((x-[0])/[1])**2)+[4]*[3]*0.5*[2]*TMath::Exp(0.5*[2]*(-2.0*[0]+[2]*[1]**2.0+2.0*x))*TMath::Erfc((-[0]+[2]*[1]**2+x)/(2.0**0.5*[1]))",0.7,1.2);
  // fGausEmg2->SetParNames("Mean","Sigma","Lambda","Amplitude","PileUpFraction");

    cLightYield -> cd();//Change current directory to cLightYield canvas
    TH1D *hScaled = (TH1D*)(hRes -> Clone()); //hScaled is a 1D histogram and we will clone hRes into hScaled
    
    //Setting hScaled details:
    hScaled -> SetLineColor(kWhite);
    hScaled -> SetXTitle("Reconstructed Energy (MeV)");
    hScaled -> GetXaxis() -> SetRangeUser(20, 40);
    hScaled -> GetYaxis() -> SetRangeUser(-0.001, 1.01);
    hScaled -> Draw("same");
    
     TF1 *gaus2 = new TF1(Form("fit_%s", hRes -> GetName()), "gaus", 0, 60); //A new 1D function called "gaus2": (name, title, xmin, xmax)
     gaus2 -> SetParameters(1., 30., bestFit -> Parameter(2)/peToMeV -> Parameter(1)); //Setting parameters of function gaus2: Parameter(2) corresponds to the Standard Deviation of our gaussian, peToMeV is the Mean of the distribution.
     
     if (color == 3 || color == 5 || color == 7){
       color++;
     } //Different colors for different configurations??
     
     gaus2 -> SetLineColor(color);
     gaus2 -> Draw("SAME A L");
 
     vhist[ival]->FitSlicesY(); //Now let's project slices of vhist along Y (all bins in x) AND FIT THEM WITH A GAUSSIAN (by default).
     TH1* dM = (TH1*)gDirectory->Get(Form("%s_1",vhist[ival]->GetName())); // A histogram named dM to store the values of mean???
     TH1* dSigma = (TH1*)gDirectory->Get(Form("%s_2",vhist[ival]->GetName()));// A histogram named dSigma to store the values of standard deviation???
     
     TH1* dSigma_frac = (TH1*)dSigma->Clone(Form("%s_Frac",dSigma->GetName())); //Get the contents of dSigma and put it into dSigma_frac
     dSigma_frac->Divide(dM); //Divide Standard Deviation (dSigma) by Mean 
     
     cEnergyResolution->cd(); //Go to Energy Resolution direction/canvas!
     //Giving a nice format to 
     dSigma_frac->SetLineColor(color);
     dSigma_frac->SetLineWidth(2);
     color++;
     dSigma_frac->Draw(ival==0?"":"SAME");
     dSigma_frac->SetXTitle("Energy (MeV)");
 
     tl -> AddEntry(gaus2, Form("%s Config", config[ival])); //Adding an ntry to TLegend tl. (Object to be added, Label associated to this object)   
     tl2 -> AddEntry(dSigma_frac, Form("%s Config", config[ival]));    
  } //end ival

  // Putting a title on our Legend Boxes:
  cLightYield -> cd();
  tl -> SetHeader(Form("Events with E_{sim} = %3.1f #pm %2.1f %%", cutEnergy, allowance*100),"C");
  tl -> Draw();

  cEnergyResolution -> cd();
  tl2 -> SetHeader("Fractional Energy Resolution","C");
  tl2 -> Draw();
  color = 1;

}

//-----------------------------------ADDED FROM 02-20 Branch--------------------------------------------------------------//

void PlotIntegrated(){//Produces a lot of plots...
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0); 
  std::vector<double> pePerMeV; //A vector for number of photoelectrons per MeV
  pePerMeV.resize(config.size()); //Resizes the pePerMeV vector to the size of config vector

  TF1* dPEdE = new TF1("dPEdE","pol1(0)"); //A function named dPdE. This function is pol1(0)
  //Defining our new Canvases:
  TCanvas* cEnergyResolution = new TCanvas("cEnergyResolution","Energy Resolution");
  TCanvas* cRecoEnergy = new TCanvas("cRecoEnergy","Reconstructed Energy");
  TCanvas* cRecoEnergyIntegral = new TCanvas("cRecoEnergyIntegral","Integral Above Threshold Reconstructed Energy");
  TCanvas* cSignalToNoise = new TCanvas("cSignalToNoise","Signal to noise");
  TCanvas* cFavorites = new TCanvas("cFavorites", "cFavorites"); 
  cFavorites -> Divide(2, 2); //more than one type of plot on the canvas
  TCanvas* cOxygenAsSignal = new TCanvas("cOxygenAsSignal","Oxygen As Signal");
  cOxygenAsSignal -> Divide(2, 2); //more than one type of plot on the canvas
  
  for(int ival = 0; ival < config.size(); ival++){ //for all configurations
      
//***************** Getting Energy Resolution:  *****************//
    vhistD[ival]->FitSlicesY(); //(Projects slices along Y) Produces "ival" TH1D histograms
    TH1* dMean = (TH1*)gDirectory->Get(Form("%s_1",vhistD[ival]->GetName())); //Accessing the generated Histograms
    TH1* dSigma = (TH1*)gDirectory->Get(Form("%s_2",vhistD[ival]->GetName()));
    dPEdE->SetRange(dMean->GetXaxis()->GetXmin(),dMean->GetXaxis()->GetXmax()); //Setting range of dPEdE
    auto fitR = dMean->Fit(dPEdE,"S0");
    pePerMeV[ival]=fitR->Parameter(1);
    // std::cout<< ival <<" Config("<<config[ival]<<" "<<pePerMeV[ival]<<std::endl;

    TH1* dSigma_frac = (TH1*)dSigma->Clone(Form("%s_Frac",dSigma->GetName()));
    dSigma_frac->Divide(dMean);
    cEnergyResolution->cd();
    dSigma_frac->SetLineColor(kBlack+ival);
    dSigma_frac->Draw(ival==0?"":"SAME");
//*************************************************************************************//
    TH1* projD = vhistD[ival]->ProjectionY(); //projD is the projection in Y of the histogram vhistD for the configuration ival
    double *xBinsEPrime = 0;
    //binning:
    if(logBinning){//if logBinning=true...
      xBinsEPrime = genLogBins(projD->GetNbinsX(),
                               projD->GetXaxis()->GetXmin()/pePerMeV[ival],
                               projD->GetXaxis()->GetXmax()/pePerMeV[ival]);
      
      projD->SetBins(projD->GetNbinsX(),
                     xBinsEPrime); //... it uses the function genLogBins to generate the binning.
    } else{//if logBinning=false...
      projD->SetBins(projD->GetNbinsX(),
                     projD->GetXaxis()->GetXmin()/pePerMeV[ival],
                     projD->GetXaxis()->GetXmax()/pePerMeV[ival]); //it DOESN'T use such function to generate the binning.
    }
    
    //Giving some format to projD histogram:
    
    projD->SetLineColor(kBlack+ival);
    projD->SetLineStyle(1);
    projD->SetLineWidth(2);
    projD->GetXaxis()->SetTitle(Form("Reconstructed Energy [MeV] for %s Config",config[ival]));
    projD->SetTitle(Form("CC #nu-d (Config:%s)", config[ival]));
    projD->GetXaxis()->SetTitleOffset(1.0);
    vhistDProj.push_back(projD); //add projD to this "array of histograms"
    cRecoEnergy->cd(); //pointer to "this" directory???
    projD->Draw(ival==0?"L NOHIST":"L NOHIST SAME");
    
    double totalD2OEntries = vhistScalars[ival]->GetBinContent(1+1); //getting the number of events in d2o from histogram vhistScalars
    projD->Scale(1.0/projD->Integral("width")); //scaling: number of entries in each bin is divided by width of the bin???
    projD->Scale(totalDCCPerTon*d2oScale[ival]); //also, we normalize by multiplying by totalDCCPerTon and mass of the volume (d2oScale) 
    // projD->Scale(totalDCCPerTon*0.67);
    std::cout<<"Integral for "<<projD->GetName()<<" "<<projD->Integral("width")<<std::endl; //prints the integral of projD histogram (for instance "Integral for hLYPED_box_dodConfig_py 2416.9")

    //Getting counts above threshold (????):
    auto projD_Integral = (TH1*)projD->Clone(Form("%s_integral",projD->GetName()));
    projD_Integral->Reset();
    for(int ibin = projD_Integral->GetNbinsX(); ibin>=1; ibin--){
      projD_Integral->SetBinContent(ibin,
                                    projD->Integral(ibin,projD->GetNbinsX(),"width"));
    }
    cRecoEnergyIntegral->cd();
    projD_Integral->GetYaxis()->SetTitle("Counts above threshold");
    projD_Integral->Draw(ival==0?"C NOHIST":"C NOHIST SAME");
    
    // Now we mess with the Oxygen info:
    
    TH1* projO = vhistO[ival]->ProjectionY(); //projection in Y os the histogram vhistO = projO
    //binning:
    if(logBinning){
      projO->SetBins(projO->GetNbinsX(),
                     xBinsEPrime);
    }else{
      projO->SetBins(projO->GetNbinsX(),
                     projO->GetXaxis()->GetXmin()/pePerMeV[ival],
                     projO->GetXaxis()->GetXmax()/pePerMeV[ival]);
    }
    
    //Giving some format to projO histogram:
    
    projO->GetXaxis()->SetTitle(Form("Reconstructed Energy [MeV] for %s Config",config[ival]));
    projO->GetXaxis()->SetTitleOffset(1.0);
    projO->SetTitle(Form("CC #nu-O (Config: %s)",config[ival]));
    
    double totalH2OEntries = vhistScalars[ival]->GetBinContent(2+1)+totalD2OEntries; //getting the number of events in h2o from histogram vhistScalars (WHY IS IT ADDING D2O ENTRIES???????)
    std::cout<<"Total D2O/H2O Entries "
             <<totalD2OEntries<<"/"<<totalH2OEntries
             <<"="<<double(totalD2OEntries)/totalH2OEntries<<std::endl; //for instance, for a 12,500 run: Total D2O/H2O Entries 6003/12500 = 0.48024 
    projO->Scale(1.0/projO->Integral("width"));//scaling: number of entries in each bin is divided by width of the bin???
    projO->Scale(totalOCCPerTon*d2oScale[ival]*totalH2OEntries/totalD2OEntries);//also, we normalize by multiplying by totalDCCPerTon and mass of the volume (d2oScale) 
    // projO->Scale(totalOCCPerTon*0.67*totalH2OEntries/totalD2OEntries);
    std::cout<<"Integral for "<<projO->GetName()<<" "<<projO->Integral("width")<<std::endl;// For instance "Integral for hLYPEO_box_dodConfig_py 396.653"

    vhistOProj.push_back(projO); // add projO to vhistOProj histogram array
    projO->SetLineColor(kBlack+ival); //giving nice format to the histogram
    projO->SetLineStyle(2);
    projO->SetLineWidth(2);
    cRecoEnergy->cd();
    projO->Draw("L NOHIST SAME");

    //Looks like getting counts above threshold for projD but now for Oxygen(????):
    auto projO_Integral = (TH1*)projO->Clone(Form("%s_integral",projO->GetName()));
    projO_Integral->Reset();
    for(int ibin = projO_Integral->GetNbinsX(); ibin>=1; ibin--){
      projO_Integral->SetBinContent(ibin,
                                    projO->Integral(ibin,projO->GetNbinsX(),"width"));
    }
    cRecoEnergyIntegral->cd();
    projO_Integral->Draw("C NOHIST SAME");

    //Signal to Noise stuff:
    cSignalToNoise->cd();
    //Defining objects and putting labels:
    auto projD_SNR = (TH1*)projD->Clone(Form("%s_SNR",projD->GetName()));
    projD_SNR->Reset();
    projD_SNR->GetYaxis()->SetTitle("#nu_{e}-d Precision");
    projD_SNR->SetTitle(Form("Config: %s #sigma/S=#sqrt{S+B}/S",config[ival]));
    auto projD_SNR2 = (TH1*)projD->Clone(Form("%s_SNR2",projD->GetName()));
    projD_SNR2->Reset();
    projD_SNR2->SetLineStyle(2);
    projD_SNR2->GetYaxis()->SetTitle("#nu_{e}-d Precision");
    projD_SNR2->SetTitle(Form("Config: %s #sigma/S=#sqrt{S+2B}/S",config[ival]));
    //Filling out these plots:
    for(int ibin=1; ibin<=projD_SNR->GetNbinsX(); ibin++){
      if((projD_Integral->GetBinContent(ibin)+projO_Integral->GetBinContent(ibin))<=0) continue;
      projD_SNR->SetBinContent(ibin, 1.0/(projD_Integral->GetBinContent(ibin)/sqrt(projD_Integral->GetBinContent(ibin)+projO_Integral->GetBinContent(ibin))));
      projD_SNR2->SetBinContent(ibin, 1.0/(projD_Integral->GetBinContent(ibin)/sqrt(projD_Integral->GetBinContent(ibin)+2*projO_Integral->GetBinContent(ibin))));
    }
    //Giving format to the graph:
    projD_SNR->GetYaxis()->SetTitleOffset(1.5);
    projD_SNR->GetXaxis()->SetRangeUser(0.1,45);
    projD_SNR->Draw(ival==0?"C NOHIST":"C NOHIST SAME");
    projD_SNR2->GetXaxis()->SetRangeUser(0.1,45);
    projD_SNR2->Draw("C NOHIST SAME");
    vhistPrecision.push_back(projD_SNR2);


    //This chunk is meant to show the favorite, but I'm showing all of them for the moment.
    //Change ival -> ifav for showing a single plot later.
    //SO THIS PLOT IS SHOWING COUNTS FOR EACH INTERACTION (deuterium,Oxygen) :
    cFavorites -> cd(1 + ival);
    TH1* hD = vhistDProj[ival]->DrawCopy(""); //using projD histogram
    TH1* hO = vhistOProj[ival]->DrawCopy(" SAME"); //using projO histogram
    TH1* hPrec = vhistPrecision[ival]->DrawCopy("C NOHIST SAME"); // using projD_SNR2
    
    //Setting Errors:
    for(int ibin = 1; ibin <= hD->GetNbinsX(); ibin++){
      hD->SetBinError(ibin,sqrt(hD->GetBinContent(ibin)*hD->GetXaxis()->GetBinWidth(ibin)));
      hO->SetBinError(ibin,sqrt(hO->GetBinContent(ibin)*hD->GetXaxis()->GetBinWidth(ibin)));
    }
    
    //Giving some format to the labels:
    hD->GetXaxis()->SetLabelFont(62);
    hD->GetXaxis()->SetTitleFont(62);
    hD->GetYaxis()->SetLabelFont(62);
    hD->GetYaxis()->SetTitleFont(62);
    hD->GetXaxis()->SetLabelSize(0.04);
    hD->GetXaxis()->SetTitleSize(0.04);
    hD->GetYaxis()->SetLabelSize(0.04);
    hD->GetYaxis()->SetTitleSize(0.04);
    
    //Giving some format to the plots:
    hD->GetYaxis()->SetTitle("Counts/MeV 4 SNS#bulletYears");
    hD->GetXaxis()->SetRangeUser(1.0,60.0);
    hD->GetYaxis()->SetRangeUser(0.0,120.0);
    hD->SetLineColor(kBlue);
    hO->SetLineColor(kRed);
    hD->SetTitle("CC #nu_{e}-d");
    hO->SetTitle("CC #nu_{e}-O");
    hO->SetLineStyle(1);
    
    hPrec->SetLineColor(kBlack);
    hPrec->SetLineStyle(1);
    hPrec->SetTitle("Precision (Above Threshold)");
    hPrec->GetXaxis()->SetRangeUser(1.0,60.0);
    hPrec->SetName("PrecCopy");
    double aScaling = 1000.0;
    hPrec->Scale(aScaling); //IDK why the scaling is 1,000 here...
    gPad -> BuildLegend(0.7,0.7,0.9,0.9);
    cFavorites->Modified();
    cFavorites->Update();
    
    //***************** I DON'T KNOW WHAT THIS CHUNK OF CODE IS FOR *****************//
    double PrecisionMin = 1.0;
    double PrecisionE = 0.0;
    for(int ibin = hPrec->GetXaxis()->FindBin(1.0); ibin<hPrec->GetXaxis()->FindBin(50.0); ibin++){
      if(hPrec->GetBinContent(ibin)/aScaling<PrecisionMin){
	PrecisionMin=hPrec->GetBinContent(ibin)/aScaling;
	PrecisionE= hPrec->GetXaxis()->GetBinCenter(ibin);
      }
    }
    std::cout<<"Precision Minimum:"<<PrecisionMin<<" at Energy Cut "<<PrecisionE<<std::endl; //For instance: Precision Minimum:0.0237453 at Energy Cut 0.836214
    TLine* tl = new TLine(gPad->GetUxmin(),PrecisionMin*aScaling,gPad->GetUxmax(),PrecisionMin*aScaling);
    // TLine* tl = new TLine(cFavorites->GetUxmin(),PrecisionMin*aScaling,cFavorites->GetUxmax(),PrecisionMin*aScaling);
    tl->SetLineWidth(2);
    tl->SetLineColor(kBlack);
    tl->SetLineStyle(2);
    tl->Draw();
    
    //TF1 *f2=new TF1("f2","x",0.0,cFavorites->GetUymax()/aScaling);
    TGaxis* gAPrecision = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0.0,gPad->GetUymax()/aScaling,510,"+L");
    // TGaxis* gAPrecision = new TGaxis(cFavorites->GetUxmax(),cFavorites->GetUymin(),cFavorites->GetUxmax(),cFavorites->GetUymax(),0.0,cFavorites->GetUymax()/aScaling,510,"+L");
    gAPrecision->SetLabelOffset(0.01);
    gAPrecision->SetTitle("Precision");
    gAPrecision->SetTitleOffset(1.2);
    gAPrecision->CenterTitle();
    //TGaxis* gAPrecision = new TGaxis(90,0,90,45,"f2",510,"-");
    gAPrecision->Draw();
    //*************************************************************************************//
    
    auto hOxygenAsSignal = (TH1*)hO->Clone("hOxygenAsSignal");
    for(int xbin = 1; xbin<=hOxygenAsSignal->GetNbinsX(); xbin++){
      hOxygenAsSignal->SetBinError(xbin,sqrt(hOxygenAsSignal->GetBinContent(xbin)+2*hD->GetBinContent(xbin))); //Setting bin errors
      
    }
    
    cOxygenAsSignal -> cd(1 + ival);
    //Setting labels
    hOxygenAsSignal->GetXaxis()->SetTitle(Form("Reconstructed Energy [MeV] for %s Config",config[ival]));
    hOxygenAsSignal->GetYaxis()->SetTitle("Counts/MeV 4 SNS#bulletYears");
    //hOxygenAsSignal->Draw();
    //Stacking Plots:
    THStack* hCCOAsSignal = new THStack("hCCOAsSignal","hCCOAsSignal");
    auto hDAsCCOBackground = (TH1*)hD->Clone("hDAsCCOBackground");
    auto hOAsCCOSignal = (TH1*)hO->Clone("hOAsCCOSignal");
    //Giving it a nice format: NEED TO KNOW WHAT EACH OF THESE GUYS IS!
    hDAsCCOBackground->SetFillColor(kGreen-6);
    hDAsCCOBackground->SetLineColor(kGreen-6);
    hDAsCCOBackground->GetXaxis()->CenterTitle();
    hDAsCCOBackground->GetYaxis()->CenterTitle();
    hOAsCCOSignal->SetFillColor(kYellow-6);
    hOAsCCOSignal->SetLineColor(kYellow-6);
    hCCOAsSignal->Add(hDAsCCOBackground);
    hCCOAsSignal->Add(hOAsCCOSignal);
    hCCOAsSignal->Draw("HIST");
    hCCOAsSignal->GetXaxis()->SetTitle(Form("Reconstructed Energy [MeV] for %s Config",config[ival]));
    cOxygenAsSignal -> Modified();

    //Setting errors for all these:
    TGraphErrors* grCCOAsSignal = new TGraphErrors(hOAsCCOSignal->GetNbinsX());
    for(int ibin = 0; ibin<grCCOAsSignal->GetN();ibin++){
        grCCOAsSignal->SetPoint(ibin,hDAsCCOBackground->GetXaxis()->GetBinCenter(ibin+1),
                                hDAsCCOBackground->GetBinContent(ibin+1)
                                + hOAsCCOSignal->GetBinContent(ibin+1)
                                );
        grCCOAsSignal->SetPointError(ibin,0.01,
                                     sqrt(2.0*hDAsCCOBackground->GetBinContent(ibin+1)
                                          + hOAsCCOSignal->GetBinContent(ibin+1)));
        grCCOAsSignal->SetMarkerStyle(20);
        grCCOAsSignal->Draw("P");
    }    

  } //end for ival
  
  //Setting Legends:
  cRecoEnergy->BuildLegend(0.6,0.6,0.85,0.85)->SetLineWidth(0);
  cRecoEnergyIntegral->BuildLegend(0.6,0.6,0.85,0.85)->SetLineWidth(0);
  cSignalToNoise->BuildLegend(0.2,0.6,0.45,0.85)->SetLineWidth(0);
  
}

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


//  void PlotFiducial(){//Position studies
//    gStyle->SetOptStat(0);
//    gStyle->SetOptTitle(0);
//    
//    TCanvas* cFid;
//    TCanvas* cIndividual; //New from Becca's PlotFiducial()
//  
//    // for (int ival = 2; ival < 3; ival++){
//    for (int ival = 0; ival < config.size(); ival++){
//      auto Sim_Tree = new TChain("Sim_Tree","Sim_Tree");
//      std::cout << fnames[ival] << ":\t Adding a " << config[ival] << " config\n";
//      Sim_Tree->Add(fnames[ival].c_str());
//      
//      simEvent *ev = new simEvent();
//      Sim_Tree -> SetBranchAddress("eventData", &ev);
//  
//      cFid = new TCanvas(Form("cFiducial_%s_Config", config[ival]), Form("cFiducial_%s_Config", config[ival]), 1500, 900); 
//      cFid -> Divide(2);
//      
//      cIndividual = new TCanvas(Form("cIndividual_%s_Config", config[ival]), Form("cIndividual_%s_Config", config[ival]), 1500, 900);//New from Becca's PlotFiducial()
//      cIndividual -> Divide(2, 2);//New from Becca's PlotFiducial() 
//     
//      std::vector<double> xIn;
//      std::vector<double> yIn;
//      std::vector<double> zIn;
//      
//      std::vector<double> xOut;
//      std::vector<double> yOut;
//      std::vector<double> zOut;
//      
//      // for (int iEvent = 0; iEvent < 500; iEvent++){ //Loop over events in config
//      for (int iEvent = 0; iEvent < Sim_Tree -> GetEntries(); iEvent++){ //Loop over events in config
//          
//  //*******************************ENERGY RESOLUTION + POSITION STUFF*******************************        
//        Sim_Tree -> GetEntry(iEvent);
//        cutEnergyF = firstCutValue; //New from Becca's PlotFiducial()
//        
//     while (cutEnergyF < maxCutValue){ //New from Becca's PlotFiducial()
//         
//        if (ev -> sourceParticleEnergy > (cutEnergyF - allowance*cutEnergyF) && ev -> sourceParticleEnergy < (cutEnergyF + allowance*cutEnergyF)){
//            
//          double eReconstructed = cutEnergyF*(ev -> numHits)/hitsPerMeV[cutEnergyF];
// //       double eReconstructed = cutEnergyF*(ev -> numHits)/hitsPerMeV[ival][cutEnergyF]; //New from Becca's PlotFiducial()
//          TVector3 pos = ev -> position0; // starting position of the event
//  	
//               if (cutEnergyF*(1 - 2*e0Sigma[cutEnergyF]/hitsPerMeV[cutEnergyF]) < eReconstructed && eReconstructed < cutEnergyF*(1 + 2*e0Sigma[cutEnergyF]/hitsPerMeV[cutEnergyF])){
// //            if (cutEnergyF*(1 - e0Sigma[cutEnergyF]/hitsPerMeV[cutEnergyF]) < eReconstructed && eReconstructed < cutEnergyF*(1 + e0Sigma[cutEnergyF]/hitsPerMeV[cutEnergyF])){ //New from Becca's PlotFiducial()
//                 xIn.push_back(pos.X());
//                 yIn.push_back(pos.Y());
//                 zIn.push_back(pos.Z());
//             } else{
//             xOut.push_back(pos.X());
//             yOut.push_back(pos.Y());
//             zOut.push_back(pos.Z());
//                 }//endelse
//  	      
//        }//endif
//        
//     cutEnergyF += cutSteps; //New from Becca's PlotFiducial()
//     
//         }//endwhile //New from Becca's PlotFiducial()
//     
//       }//END for(iEvent) 
//  
//      
//      //New from Becca's PlotFiducial():
//     TH1D *hInX = new TH1D(Form("hInX_%s_Config", config[ival]), Form("hInX_%s_Config", config[ival]), 100, -1000, 1000);
//     hInX -> SetLineColor(kRed);
//     
//      TH1D *hOutX = new TH1D(Form("hOutX_%s_Config", config[ival]), Form("hOutX_%s_Config", config[ival]), 100, -1000, 1000);
//      hOutX -> SetXTitle("X Position (mm)");
//      hOutX -> SetYTitle("Num. of Events");
// 
//     TH1D *hInY = new TH1D(Form("hInY_%s_Config", config[ival]), Form("hInY_%s_Config", config[ival]), 100, -1000, 1000);
//     hInY -> SetLineColor(kRed);
//     
//      TH1D *hOutY = new TH1D(Form("hOutY_%s_Config", config[ival]), Form("hOutY_%s_Config", config[ival]), 100, -1000, 1000);
//      hOutY -> SetXTitle("Y Position (mm)");
//      hOutX -> SetYTitle("Num. of Events");
// 
// 
//     TH1D *hInZ = new TH1D(Form("hInZ_%s_Config", config[ival]), Form("hInZ_%s_Config", config[ival]), 200, -2000, 2000);
//     hInZ -> SetLineColor(kRed);
//     
//      TH1D *hOutZ = new TH1D(Form("hOutZ_%s_Config", config[ival]), Form("hOutZ_%s_Config", config[ival]), 200, -2000, 2000);
//      hOutZ -> SetXTitle("Z Position (mm)");
//      hOutX -> SetYTitle("Num. of Events");
//     
//  
//      int numOut = xOut.size();
//      double Xout[numOut], Yout[numOut], Zout[numOut];
//      for (int i = 0; i < numOut; i++){
//       hOutX -> Fill(xOut[i]);//New from Becca's PlotFiducial()
//       hOutY -> Fill(yOut[i]);//New from Becca's PlotFiducial()
//       hOutZ -> Fill(zOut[i]);//New from Becca's PlotFiducial()
//        Xout[i] = xOut[i];
//        Yout[i] = yOut[i];
//        Zout[i] = zOut[i];
//       }
//      
//      int numIn = xIn.size();
//      double Xin[numIn], Yin[numIn], Zin[numIn];
//      for (int i = 0; i < numIn; i++){
//       hInX -> Fill(xIn[i]);//New from Becca's PlotFiducial()
//       hInY -> Fill(yIn[i]);//New from Becca's PlotFiducial()
//       hInZ -> Fill(zIn[i]);//New from Becca's PlotFiducial()
//        Xin[i] = xIn[i];
//        Yin[i] = yIn[i];
//        Zin[i] = zIn[i];
//      }
// 
//      cFid -> cd(1);
//       TGraph *xyOut = new TGraph(numOut, Xout, Yout);
//       xyOut -> SetMarkerColor(kBlue);
//       xyOut -> SetName("xyOut");//New from Becca's PlotFiducial()
//       xyOut -> SetTitle("XY; X pos (mm); Y pos (mm)");//New from Becca's PlotFiducial()
//       xyOut -> Draw("A P *");
//      TGraph *xyIn = new TGraph(numIn, Xin, Yin);
//      xyIn -> SetMarkerColor(kRed);
//      xyIn -> SetName("xyIn");//New from Becca's PlotFiducial()
//      xyIn -> SetTitle("XY; X pos (mm); Y pos (mm)");//New from Becca's PlotFiducial()
//      xyIn -> Draw("P * same");
//      xyIn -> Draw("A P *");
//      
//  
//      cFid -> cd(2);
//       TGraph *xzOut = new TGraph(numOut, Xout, Zout);
//       xzOut -> SetName("xzOut");//New from Becca's PlotFiducial()
//       xzOut -> SetTitle("XZ; X pos (mm); Z pos (mm)");//New from Becca's PlotFiducial()
//       xzOut -> SetMarkerColor(kBlue);
//       xzOut -> Draw("A P *");
//      TGraph *xzIn = new TGraph(numIn, Xin, Zin);
//      xzIn -> SetName("xzIn"); //New from Becca's PlotFiducial()
//      xzIn -> SetTitle("XZ; X pos (mm); Z pos (mm)"); //New from Becca's PlotFiducial()
//      xzIn -> SetMarkerColor(kRed);
//      xzIn -> Draw("P * same");
//      xzIn -> Draw("A P *");
//      
//      //New from Becca's PlotFiducial()
//      cIndividual -> cd(1);
//      hOutX -> Draw();
//      hInX -> Draw("same");
//      hInX -> Draw();
// 
//      cIndividual -> cd(2);
//      hOutY -> Draw();
//      hInY -> Draw("same");
//      hInY -> Draw();
// 
//      cIndividual -> cd(3);
//      hOutZ -> Draw();
//      hInZ -> Draw("same");
//      hInZ -> Draw();
//     
// 
//    }//END for(ival)
//  }//END PlotFiducial() 

 void PlotFiducial(){//Position studies (OLD, WORKING VERSION)
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   
   TCanvas* cFid;
 
   // for (int ival = 2; ival < 3; ival++){
   for (int ival = 0; ival < config.size(); ival++){
     auto Sim_Tree = new TChain("Sim_Tree","Sim_Tree");
     std::cout << fnames[ival] << ":\t Adding a " << config[ival] << " config\n";
     Sim_Tree->Add(fnames[ival].c_str());
     
     simEvent *ev = new simEvent();
     Sim_Tree -> SetBranchAddress("eventData", &ev);
 
     cFid = new TCanvas(Form("cFiducial_%s_Config", config[ival]), Form("cFiducial_%s_Config", config[ival]));
     cFid -> Divide(2);
     
     std::vector<double> xIn;
     std::vector<double> yIn;
     std::vector<double> zIn;
     
     std::vector<double> xOut;
     std::vector<double> yOut;
     std::vector<double> zOut;
     
     // for (int iEvent = 0; iEvent < 500; iEvent++){ //Loop over events in config
     for (int iEvent = 0; iEvent < Sim_Tree -> GetEntries(); iEvent++){ //Loop over events in config
       Sim_Tree -> GetEntry(iEvent);
       if (ev -> sourceParticleEnergy > (cutEnergy - allowance*cutEnergy) && ev -> sourceParticleEnergy < (cutEnergy + allowance*cutEnergy)){
 	double eReconstructed = cutEnergy*(ev -> numHits)/hitsPerMeV[cutEnergy];
 	TVector3 pos = ev -> position0; // starting position of the event
 	if (cutEnergy*(1 - 2*e0Sigma[cutEnergy]/hitsPerMeV[cutEnergy]) < eReconstructed && eReconstructed < cutEnergy*(1 + 2*e0Sigma[cutEnergy]/hitsPerMeV[cutEnergy])){
 	  xIn.push_back(pos.X());
 	  yIn.push_back(pos.Y());
 	  zIn.push_back(pos.Z());
 	} else{
 	  xOut.push_back(pos.X());
 	  yOut.push_back(pos.Y());
 	  zOut.push_back(pos.Z());
 	}
       }
     }
 
     int numOut = xOut.size();
     double Xout[numOut], Yout[numOut], Zout[numOut];
     for (int i = 0; i < numOut; i++){
       Xout[i] = xOut[i];
       Yout[i] = yOut[i];
       Zout[i] = zOut[i];
     }
     int numIn = xIn.size();
     double Xin[numIn], Yin[numIn], Zin[numIn];
     for (int i = 0; i < numIn; i++){
       Xin[i] = xIn[i];
       Yin[i] = yIn[i];
       Zin[i] = zIn[i];
     }
 
     cFid -> cd(1);
     TGraph *xyOut = new TGraph(numOut, Xout, Yout);
     xyOut -> SetMarkerColor(kBlue);
     xyOut -> Draw("A P *");
     TGraph *xyIn = new TGraph(numIn, Xin, Yin);
     xyIn -> SetMarkerColor(kRed);
     xyIn -> Draw("P * same");
 
     cFid -> cd(2);
     TGraph *xzOut = new TGraph(numOut, Xout, Zout);
     xzOut -> SetMarkerColor(kBlue);
     xzOut -> Draw("A P *");
     TGraph *xzIn = new TGraph(numIn, Xin, Zin);
     xzIn -> SetMarkerColor(kRed);
     xzIn -> Draw("P * same");
 
   }
 }
