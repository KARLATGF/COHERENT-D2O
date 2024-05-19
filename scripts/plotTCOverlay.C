#include "simEvent.h"
#include <vector>
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

//gSystem->AddIncludePath(" -I/Users/frn/Documents/Projects/COHERENT/d2o/sim/G4d2o/simEvent  ")
std::vector<TChain*> vfin;
std::vector<TH2*> vhist;
std::vector<TH2*> vhistD;
std::vector<TH2*> vhistO;
std::vector<TH1*> vhistScalars;
std::vector<TH1*> vhistDProj;
std::vector<TH1*> vhistOProj;
std::vector<TH1*> vhistPrecision;

std::vector<std::string> fnames = {"data/prodTC5cm/Sim_D2ODetector00*.root",
  "data/prodTC10cm/Sim_D2ODetector00*.root",
  "data/prodTC15cm/Sim_D2ODetector00*.root",
  "data/prodTC20cm/Sim_D2ODetector00*.root"};

std::vector<int> tcval = {5, 10, 15, 20};

double xsntotal_dcc = 5.5E-41;
double xsntotal_Occ = 8.6696988e-42;

double SNSHoursPerYear = 5000.0;
double chargePerPulse = 21.59E-6;
double POTperhour = 1.4/1.3*chargePerPulse/1.602E-19*60*3600;
//double POTperday = 7E20;
double POTperPulse = 1.4/1.3*chargePerPulse/1.602E-19;
double POTperYear = POTperhour*SNSHoursPerYear;
double pionsperPOT = 0.09; // @1010MeV/p
double enuPerYear = POTperYear*pionsperPOT;
double areaAt20m = 4*TMath::Pi()*TMath::Power(20.0*100.0,2.0);
double enuFluxat20m = enuPerYear/areaAt20m*2.0;
double deuteronsIn1ton = 1E6/(4.0+16.0)*TMath::Na()*2.0;
double oxygenIn1ton = 1E6/(4.0+16.0)*TMath::Na()*1.0;
double totalDCCPerTon = deuteronsIn1ton*xsntotal_dcc*enuFluxat20m;
double totalOCCPerTon = deuteronsIn1ton*xsntotal_Occ*enuFluxat20m;

bool logBinning = false;

double* genLogBins(int nbins, double axismin, double axismax)
{
  double* dbins = new double[nbins+1];
  for(int ibin = 0; ibin <= nbins; ibin++)
  {
    dbins[ibin] = axismin*pow(axismax/axismin, ((double)ibin)/nbins);
  }
  return dbins;
}

void plotTCOverlay(){
  int ival = 0;
  int nBinsPE = 50;
  int nBinsE = 50;

  for(auto fname : fnames){
    auto Sim_Tree = new TChain("Sim_Tree","Sim_Tree");
    Sim_Tree->Add(fname.c_str());
    Sim_Tree->Print();
    TTreeReader t(Sim_Tree);
    TTreeReaderValue<simEvent> ev(t,"eventData");

    double *xBinsPE = genLogBins(nBinsPE,0.1,2000);
    double *xBinsE = genLogBins(nBinsE,0.1,55);

    TString hNameLY(Form("hLYPE_TC%dcm",tcval[ival]));
    TH2* hLY = 0;
    if(logBinning){
      hLY = new TH2D(hNameLY.Data(),hNameLY.Data(),nBinsE,xBinsE,nBinsPE,xBinsPE);
    }else{
      hLY = new TH2D(hNameLY.Data(),hNameLY.Data(),100,0,50,100,0,2000);
    }
    vhist.push_back(hLY);
    
    TString hNameD(Form("hLYPED_TC%dcm",tcval[ival]));
    TH2* hLYD = 0;
    if(logBinning){
      hLYD = new TH2D(hNameD.Data(),hNameD.Data(),nBinsE,xBinsE,nBinsPE,xBinsPE);
    }else{
      hLYD = new TH2D(hNameD.Data(),hNameD.Data(),100,0,50,100,0,2000);
    }
    vhistD.push_back(hLYD);
    
    TString hNameO(Form("hLYPEO_TC%dcm",tcval[ival]));
    TH2* hLYO = 0;
    if(logBinning){
      hLYO = new TH2D(hNameO.Data(),hNameO.Data(),nBinsE,xBinsE,nBinsPE,xBinsPE);
    }else{
      hLYO = new TH2D(hNameO.Data(),hNameO.Data(),100,0,50,100,0,2000);
    }
    vhistO.push_back(hLYO);
    TString hNameScalar(Form("hScalars_TC%dcm",tcval[ival]));
    TH1* hScalar = new TH1D(hNameScalar.Data(),hNameScalar.Data(),10,0,10.0);
    vhistScalars.push_back(hScalar);
    
    hLYD->SetDirectory(gROOT);
    hLYO->SetDirectory(gROOT);
    hScalar->SetDirectory(gROOT);
    
//    Sim_Tree->Draw(Form("numHits>>%s",hNameD.Data()),"WeightD()*(vol0==1)","",200000);
//    Sim_Tree->Draw(Form("numHits>>%s",hNameO.Data()),"WeightO()*(vol0>0)","",200000);
    
    int maxevents = 200000;
    int ievent = 0;
    
    while(t.Next()&&ievent<maxevents){
      ievent++;
      //if( !(40.0<sourceParticleEnergy&&sourceParticleEnergy<50) ) continue;
      int ndetpe = 0;
      hScalar->Fill(ev->vol0);
      if(ev->vol0==1){
        hLY->Fill(ev->sourceParticleEnergy,ev->numHits);
        hLYD->Fill(ev->sourceParticleEnergy,ev->numHits,ev->WeightD());
      }
      if(ev->vol0>0){
        hLYO->Fill(ev->sourceParticleEnergy,ev->numHits,ev->WeightO());
      }
    }

//    vhist.push_back((TH1*)vfin.back()->Get(Form("hPEperEvent_%dcmLW",tcval[ival])));
//    vhist.back()->SetLineColor(kBlack+ival);
//    vhist.back()->SetLineWidth(2);
//    vhist.back()->Rebin(5);
//    vhist.back()->SetTitle(Form("LW %d cm",tcval[ival]));
//    vhist.back()->Draw(ival==0?"":"same");
    ival++;
    
  }
  
  TFile* fOut = TFile::Open("TCStudyV2.root","RECREATE");
  for(auto h : vhist){
    h->SetDirectory(fOut);
  }
  for(auto h : vhistD){
    h->SetDirectory(fOut);
  }
  for(auto h : vhistO){
    h->SetDirectory(fOut);
  }
  for(auto h : vhistScalars){
    h->SetDirectory(fOut);
  }

  fOut->Write();
  fOut->Close();
//  int bin3540 = vhistD[0]->GetXaxis()->FindBin(37.0);
//  for(int i = 0; i<vhistD.size(); i++){
//    TH1* hy = vhistD[i]->ProjectionY(Form("%s_py%d",vhistD[i]->GetName(),bin3540),bin3540,bin3540);
//    hy->SetLineColor(kBlack+i);
//    hy->Draw(i==0?"":"same");
//  }
  
  
}

void LoadHists(){
  TFile* fIn = 0;
  if(logBinning){
    fIn = TFile::Open("TCStudyV1Log.root");
  }else{
    fIn = TFile::Open("TCStudyV1.root");
  }
  
  int ival = 0;
  for(auto tc : tcval){
    TString hNamePE(Form("hLYPE_TC%dcm",tcval[ival]));
    auto hlype = (TH2*)fIn->Get(hNamePE.Data());
    hlype->Scale(1.0,"width");
    vhist.push_back(hlype);
    
    TString hNameD(Form("hLYPED_TC%dcm",tcval[ival]));
    auto hlyped = (TH2*)fIn->Get(hNameD.Data());
    hlyped->Scale(1.0,"width");
    vhistD.push_back(hlyped);
    
    TString hNameO(Form("hLYPEO_TC%dcm",tcval[ival]));
    auto hlypeo = (TH2*)fIn->Get(hNameO.Data());
    hlypeo->Scale(1.0,"width");
    vhistO.push_back(hlypeo);
    
    TString hNameScalar(Form("hScalars_TC%dcm",tcval[ival]));
    auto hscalars = (TH1*)fIn->Get(hNameScalar.Data());
    vhistScalars.push_back(hscalars);
    
    ival++;
  }
}

double getMaxValueInRange(TF1* tg, double xmin, double xmax){
    double mval=tg->Eval(xmin);
    double mx=xmin;
    double ix = xmin;
    for(ix = xmin; ix<=xmax;ix+=(xmax-xmin)/1000.0){
        if(tg->Eval(ix)>mval){
            mval = tg->Eval(ix);
            mx =ix;
        }
    }
    return mx;
}

void PlotResolution(){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TCanvas* cLightLightYield = new TCanvas("cLightLightYield","cLightLightYield");
  TCanvas* cLightEReco = new TCanvas("cLightEReco","cLightEReco");

  TF1* fGausEmg2 = new TF1("fGausEmg2","[3]/([1]*sqrt(2*pi))*exp(-0.5*((x-[0])/[1])**2)+[4]*[3]*0.5*[2]*TMath::Exp(0.5*[2]*(-2.0*[0]+[2]*[1]**2.0+2.0*x))*TMath::Erfc((-[0]+[2]*[1]**2+x)/(2.0**0.5*[1]))",0.7,1.2);
  fGausEmg2->SetParNames("Mean","Sigma","Lambda","Amplitude","PileUpFraction");
  double fgeMean = 1000.0;
  fGausEmg2->SetParameters(fgeMean,0.1*fgeMean,0.5/fgeMean,0.001,0.1);
  
  TF1* fGausEmg3 = new TF1("fGausEmg2","[3]/([1]*sqrt(2*pi))*exp(-0.5*((x-[0])/[1])**2)+[4]*[3]*0.5*[2]*TMath::Exp(0.5*[2]*(-2.0*[0]+[2]*[1]**2.0+2.0*x))*TMath::Erfc((-[0]+[2]*[1]**2+x)/(2.0**0.5*[1]))",0.7,1.2);
  fGausEmg3->SetParNames("Mean","Sigma","Lambda","Amplitude","PileUpFraction");

  for(int ival = 0; ival < tcval.size(); ival++){
    int pbin = vhist[ival]->GetXaxis()->FindBin(45.0);

    TH1* hLYProj = vhist[ival]->ProjectionY(Form("%s_py%d",vhist[ival]->GetName(),pbin),pbin,pbin);
    hLYProj->SetLineColor(kBlack+ival);
    hLYProj->GetXaxis()->SetTitle("Light Yield [photoelectrons]");
    hLYProj->GetXaxis()->SetTitleOffset(1.0);
    hLYProj->GetXaxis()->CenterTitle(1.0);
    hLYProj->SetTitle(Form("TC %d cm",tcval[ival]));
    hLYProj->Sumw2();
    hLYProj->Rebin(2);
    hLYProj->Scale(1.0/hLYProj->Integral("width"));
    cLightLightYield->cd();
    hLYProj->Draw(ival==0?"":"SAME");

    TH1* hLYProj_cal = (TH1*)hLYProj->Clone(Form("%s_cal",hLYProj->GetName()));

    fGausEmg2->SetLineColor(hLYProj->GetLineColor());
    hLYProj->Fit(fGausEmg2);
    //double fitmean = fGausEmg2->GetParameter(0);
    double fitmean = getMaxValueInRange(fGausEmg2,600.0,1400.0);
      std::cout<<"DistMean: " << fGausEmg2->GetParameter(0) << " Peak Position "<<fitmean<<std::endl;
    double calPEperkeV = fitmean/vhist[ival]->GetXaxis()->GetBinCenter(pbin);
    hLYProj_cal->SetBins(hLYProj->GetNbinsX(),
                         hLYProj->GetXaxis()->GetXmin()/calPEperkeV,
                         hLYProj->GetXaxis()->GetXmax()/calPEperkeV);
    cLightEReco->cd();
    hLYProj_cal->Scale(1.0/hLYProj_cal->Integral("width"));
    hLYProj_cal->GetXaxis()->SetTitle("Reconstructed Energy [MeV]");
    //hLYProj_cal->Draw(ival==0?"L NOHIST":"L NOHIST SAME");
    hLYProj_cal->Draw(ival==0?"L NOHIST":"L NOHIST SAME");

    double fgeMean3 = vhist[ival]->GetXaxis()->GetBinCenter(pbin);
    fGausEmg3->SetParameters(fgeMean3,0.1*fgeMean3,0.5/fgeMean3,0.5,0.1);
    fGausEmg3->SetLineColor(hLYProj_cal->GetLineColor());
    hLYProj_cal->Fit(fGausEmg3);
    
  }
  
  TLegend* tl = gPad->BuildLegend();
  
  tl->SetLineWidth(0);
}

void PlotIntegrated(){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  std::vector<double> pePerMeV;
  pePerMeV.resize(tcval.size());

  TF1* dPEdE = new TF1("dPEdE","pol1(0)");

  TCanvas* cRecoEnergy = new TCanvas("cRecoEnergy","Reconstructed Energy");
  TCanvas* cRecoEnergyIntegral = new TCanvas("cRecoEnergyIntegral","Integral Above Threshold Reconstructed Energy");
  TCanvas* cEnergyResolution = new TCanvas("cEnergyResolution","Energy Resolution");

  TCanvas* cSignalToNoise = new TCanvas("cSignalToNoise","Signal to noise");

  for(int ival = 0; ival < tcval.size(); ival++){

    vhistD[ival]->FitSlicesY();
    TH1* dMean = (TH1*)gDirectory->Get(Form("%s_1",vhistD[ival]->GetName()));
    TH1* dSigma = (TH1*)gDirectory->Get(Form("%s_2",vhistD[ival]->GetName()));
    dPEdE->SetRange(dMean->GetXaxis()->GetXmin(),dMean->GetXaxis()->GetXmax());
    auto fitR = dMean->Fit(dPEdE,"S0");
    pePerMeV[ival]=fitR->Parameter(1);
    std::cout<<ival<<" TC("<<tcval[ival]<<" "<<pePerMeV[ival]<<std::endl;
    
    TH1* dSigma_frac = (TH1*)dSigma->Clone(Form("%s_Frac",dSigma->GetName()));
    dSigma_frac->Divide(dMean);
    cEnergyResolution->cd();
    dSigma_frac->SetLineColor(kBlack+ival);
    //dSigma_frac->Draw(ival==0?"C NOHIST":"C NOHIST SAME");
    dSigma_frac->Draw(ival==0?"":"SAME");
    TH1* projD = vhistD[ival]->ProjectionY();
    double *xBinsEPrime = 0;
    if(logBinning){
      xBinsEPrime = genLogBins(projD->GetNbinsX(),
                               projD->GetXaxis()->GetXmin()/pePerMeV[ival],
                               projD->GetXaxis()->GetXmax()/pePerMeV[ival]);

      projD->SetBins(projD->GetNbinsX(),
                     xBinsEPrime);
    }else{
      projD->SetBins(projD->GetNbinsX(),
                     projD->GetXaxis()->GetXmin()/pePerMeV[ival],
                     projD->GetXaxis()->GetXmax()/pePerMeV[ival]);
    }
    projD->SetLineColor(kBlack+ival);
    projD->SetLineStyle(1);
    projD->SetLineWidth(2);
    projD->GetXaxis()->SetTitle("Reconstructed Energy [MeV]");
    projD->SetTitle(Form("CC #nu-d (TC:%d cm)",tcval[ival]));
    projD->GetXaxis()->SetTitleOffset(1.0);
    vhistDProj.push_back(projD);
    cRecoEnergy->cd();
    projD->Draw(ival==0?"L NOHIST":"L NOHIST SAME");
    
    double totalD2OEntries = vhistScalars[ival]->GetBinContent(1+1);
    projD->Scale(1.0/projD->Integral("width"));
    projD->Scale(totalDCCPerTon*1.3);
    std::cout<<"Integral for "<<projD->GetName()<<" "<<projD->Integral("width")<<std::endl;

    auto projD_Integral = (TH1*)projD->Clone(Form("%s_integral",projD->GetName()));
    projD_Integral->Reset();
    for(int ibin = projD_Integral->GetNbinsX(); ibin>=1; ibin--){
      projD_Integral->SetBinContent(ibin,
                                    projD->Integral(ibin,projD->GetNbinsX(),"width"));
      //std::cout<<ibin<<" "<<projD->Integral(ibin,projD->GetNbinsX(),"width")<<std::endl;
    }
    cRecoEnergyIntegral->cd();
    projD_Integral->GetYaxis()->SetTitle("Counts above threshold");
    projD_Integral->Draw(ival==0?"C NOHIST":"C NOHIST SAME");
    
    
    TH1* projO = vhistO[ival]->ProjectionY();
    if(logBinning){
      projO->SetBins(projO->GetNbinsX(),
                     xBinsEPrime);
    }else{
      projO->SetBins(projO->GetNbinsX(),
                     projO->GetXaxis()->GetXmin()/pePerMeV[ival],
                     projO->GetXaxis()->GetXmax()/pePerMeV[ival]);
    }
    projO->GetXaxis()->SetTitle("Reconstructed Energy [MeV]");
    projO->GetXaxis()->SetTitleOffset(1.0);
    projO->SetTitle(Form("CC #nu-O (TC:%d cm)",tcval[ival]));
    double totalH2OEntries = vhistScalars[ival]->GetBinContent(2+1)+totalD2OEntries;
    std::cout<<"Total D2O/H2O Entries "
             <<totalD2OEntries<<"/"<<totalH2OEntries
             <<"="<<double(totalD2OEntries)/totalH2OEntries<<std::endl;
    projO->Scale(1.0/projO->Integral("width"));
    projO->Scale(totalOCCPerTon*1.3*totalH2OEntries/totalD2OEntries);
    std::cout<<"Integral for "<<projO->GetName()<<" "<<projO->Integral("width")<<std::endl;

    vhistOProj.push_back(projO);
    projO->SetLineColor(kBlack+ival);
    projO->SetLineStyle(2);
    projO->SetLineWidth(2);
    cRecoEnergy->cd();
    projO->Draw("L NOHIST SAME");

    auto projO_Integral = (TH1*)projO->Clone(Form("%s_integral",projO->GetName()));
    projO_Integral->Reset();
    for(int ibin = projO_Integral->GetNbinsX(); ibin>=1; ibin--){
      projO_Integral->SetBinContent(ibin,
                                    projO->Integral(ibin,projO->GetNbinsX(),"width"));
    }
    cRecoEnergyIntegral->cd();
    projO_Integral->Draw("C NOHIST SAME");

    cSignalToNoise->cd();

    auto projD_SNR = (TH1*)projD->Clone(Form("%s_SNR",projD->GetName()));
    projD_SNR->Reset();
    projD_SNR->GetYaxis()->SetTitle("#nu_{e}-d Precision");
    projD_SNR->SetTitle(Form("TC:%d cm #sigma/S=#sqrt{S+B}/S",tcval[ival]));
    auto projD_SNR2 = (TH1*)projD->Clone(Form("%s_SNR2",projD->GetName()));
    projD_SNR2->Reset();
    projD_SNR2->SetLineStyle(2);
    projD_SNR2->GetYaxis()->SetTitle("#nu_{e}-d Precision");
    projD_SNR2->SetTitle(Form("TC:%d cm #sigma/S=#sqrt{S+2B}/S",tcval[ival]));
    for(int ibin=1; ibin<=projD_SNR->GetNbinsX(); ibin++){
      if((projD_Integral->GetBinContent(ibin)+projO_Integral->GetBinContent(ibin))<=0) continue;
//      projD_SNR->SetBinContent(ibin,projD_Integral->GetBinContent(ibin)/sqrt(projD_Integral->GetBinContent(ibin)+projO_Integral->GetBinContent(ibin)));
//      projD_SNR2->SetBinContent(ibin,projD_Integral->GetBinContent(ibin)/sqrt(projD_Integral->GetBinContent(ibin)+2*projO_Integral->GetBinContent(ibin)));
      projD_SNR->SetBinContent(ibin,1.0/(projD_Integral->GetBinContent(ibin)/sqrt(projD_Integral->GetBinContent(ibin)+projO_Integral->GetBinContent(ibin))));
      projD_SNR2->SetBinContent(ibin,1.0/(projD_Integral->GetBinContent(ibin)/sqrt(projD_Integral->GetBinContent(ibin)+2*projO_Integral->GetBinContent(ibin))));
    }
    projD_SNR->GetYaxis()->SetTitleOffset(1.5);
    projD_SNR->GetXaxis()->SetRangeUser(0.1,45);
    projD_SNR->Draw(ival==0?"C NOHIST":"C NOHIST SAME");
    projD_SNR2->GetXaxis()->SetRangeUser(0.1,45);
    projD_SNR2->Draw("C NOHIST SAME");
    vhistPrecision.push_back(projD_SNR2);

  }
  cRecoEnergy->BuildLegend(0.6,0.6,0.85,0.85)->SetLineWidth(0);
  cRecoEnergyIntegral->BuildLegend(0.6,0.6,0.85,0.85)->SetLineWidth(0);
  cSignalToNoise->BuildLegend(0.2,0.6,0.45,0.85)->SetLineWidth(0);
  
  auto cFavorites = new TCanvas("cFavorites","cFavorites");
  int ifav = 1; // TC = 10cm
//  TH1* hD = vhistDProj[ifav]->DrawCopy("C NOHIST");
//  TH1* hO = vhistOProj[ifav]->DrawCopy("C NOHIST SAME");
  TH1* hD = vhistDProj[ifav]->DrawCopy("");
  TH1* hO = vhistOProj[ifav]->DrawCopy(" SAME");
  TH1* hPrec = vhistPrecision[ifav]->DrawCopy("C NOHIST SAME");
  
  for(int ibin = 1; ibin <= hD->GetNbinsX(); ibin++){
    hD->SetBinError(ibin,sqrt(hD->GetBinContent(ibin)*hD->GetXaxis()->GetBinWidth(ibin)));
    hO->SetBinError(ibin,sqrt(hO->GetBinContent(ibin)*hD->GetXaxis()->GetBinWidth(ibin)));
  }
  hD->GetXaxis()->SetLabelFont(62);
  hD->GetXaxis()->SetTitleFont(62);
  hD->GetYaxis()->SetLabelFont(62);
  hD->GetYaxis()->SetTitleFont(62);
  hD->GetXaxis()->SetLabelSize(0.04);
  hD->GetXaxis()->SetTitleSize(0.04);
  hD->GetYaxis()->SetLabelSize(0.04);
  hD->GetYaxis()->SetTitleSize(0.04);
  
  hD->GetYaxis()->SetTitle("Counts/MeV (1.3t D_{2}O) 2 SNS#bulletYears");
  hD->GetXaxis()->SetRangeUser(1.0,60.0);
  hD->GetYaxis()->SetRangeUser(0.0,120.0);
  hD->SetLineColor(kBlue);
  hO->SetLineColor(kRed);
  hD->SetTitle("CC #nu_{e}-d");
  hO->SetTitle("CC #nu_{e}-O");
  hO->SetLineStyle(1);
  hPrec->SetLineColor(kBlack);
  hPrec->SetLineStyle(1);
  hPrec->SetTitle("Precision (Above Threshold)");
  hPrec->GetXaxis()->SetRangeUser(1.0,60.0);
  hPrec->SetName("PrecCopy");
  double aScaling = 1000.0;
  hPrec->Scale(aScaling);
  cFavorites->Modified();
  cFavorites->Update();
  TLine* tl = new TLine(cFavorites->GetUxmin(),hPrec->GetMinimum(),cFavorites->GetUxmax(),hPrec->GetMinimum());
  tl->SetLineWidth(2);
  tl->SetLineColor(kBlack);
  tl->SetLineStyle(2);
  tl->Draw();
  
  //TF1 *f2=new TF1("f2","x",0.0,cFavorites->GetUymax()/aScaling);
  TGaxis* gAPrecision = new TGaxis(cFavorites->GetUxmax(),cFavorites->GetUymin(),cFavorites->GetUxmax(),cFavorites->GetUymax(),0.0,cFavorites->GetUymax()/aScaling,510,"+L");
  gAPrecision->SetLabelOffset(0.01);
  gAPrecision->SetTitle("Precision");
  gAPrecision->SetTitleOffset(1.2);
  gAPrecision->CenterTitle();
  //TGaxis* gAPrecision = new TGaxis(90,0,90,45,"f2",510,"-");
  gAPrecision->Draw();
    
  auto hOxygenAsSignal = (TH1*)hO->Clone("hOxygenAsSignal");
    for(int xbin = 1; xbin<=hOxygenAsSignal->GetNbinsX(); xbin++){
        hOxygenAsSignal->SetBinError(xbin,sqrt(hOxygenAsSignal->GetBinContent(xbin)+2*hD->GetBinContent(xbin)));
        
    }
    TCanvas* cOxygenAsSignal = new TCanvas("cOxygenAsSignal","Oxygen As Signal");
    hOxygenAsSignal->GetXaxis()->SetTitle("Reconstructed Energy [MeV]");
    hOxygenAsSignal->GetYaxis()->SetTitle("Counts/MeV (1.3t D_{2}O) 2 SNS#bulletYears");
    //hOxygenAsSignal->Draw();
    THStack* hCCOAsSignal = new THStack("hCCOAsSignal","hCCOAsSignal");
    auto hDAsCCOBackground = (TH1*)hD->Clone("hDAsCCOBackground");
    auto hOAsCCOSignal = (TH1*)hO->Clone("hOAsCCOSignal");
    hDAsCCOBackground->SetFillColor(kGreen-6);
    hDAsCCOBackground->SetLineColor(kGreen-6);
    hDAsCCOBackground->GetXaxis()->CenterTitle();
    hDAsCCOBackground->GetYaxis()->CenterTitle();
    hOAsCCOSignal->SetFillColor(kYellow-6);
    hOAsCCOSignal->SetLineColor(kYellow-6);
    hCCOAsSignal->Add(hDAsCCOBackground);
    hCCOAsSignal->Add(hOAsCCOSignal);
    hCCOAsSignal->Draw("HIST");
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
    
}
