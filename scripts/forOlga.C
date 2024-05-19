// THIS SCRIPT'S IS BASED IN BECCA'S compareConfigurations.C
// and it's repurposed to analize PMT on Top Geometry.

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


#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <iostream>
#include <fstream>
#include <array>
#include <TPad.h>
#include <TSystem.h>
#include <math.h>

#include <TAttLine.h>

#include "TPaveStats.h"
#include "TPaveLabel.h"


// Vectors here are fancy "arrays"
std::vector<TChain*> vfin;
std::vector<TH2*> vhist; //an array of histograms
std::vector<TH2*> vhistD; 
std::vector<TH2*> vhistO; 
std::vector<TH1*> vhistScalars; 


// //------Analyzing one single configuration:
std::vector<std::string> fnames = {"data/for_Olga/Sim_D2ODetector*.root"};

// //------Analyzing one single configuration:
 std::vector<const char*> config = {"PMTs_on_Top"};

// //------Analyzing one single configuration:
std::vector<double> d2oScale = {0.592};

// Create ROOT file to store histograms:
std::unique_ptr<TFile> MCanalysisFile( TFile::Open("light_yield_deuterium_top.root", "RECREATE") );


//------------------------------------------------------------------


void generatePlots(){ //Reading in the data...
    
    int ival = 0;
    
  for(auto fname : fnames){ //For all listed configurations do all of the following:
      
      
    auto Sim_Tree = new TChain("Sim_Tree","Sim_Tree"); //("Default name of TTree","Title")
    
    std::cout << fname << ":\t Adding a " << config[ival] << " config\n"; //Prints the name of the configuration we are working with
    Sim_Tree->Add(fname.c_str()); // fname.c_str() returns pointer to fname: this is, it returns the name of the configuration. "Add" adds a branch to the tree (??)
    // Sim_Tree->Print();

    simEvent *ev = new simEvent(); //A new object (ev) from class simEvent. "ev" represents an event
    TVector3 pos = ev -> position0;// INITIAL POSITION
    
    Sim_Tree -> SetBranchAddress("eventData", &ev); //Adding "eventData" branch to Sim_Tree (????). &ev is the address of the object "ev"

    
    //Defining the Histograms. NOT filling them yet:
    
    TString hNameD(Form("hLYPED_%s_Config",config[ival]));
    TH2* hLYD = 0;
    hLYD = new TH2D(hNameD.Data(),hNameD.Data(),55,0,55,200,0,2000);
    vhistD.push_back(hLYD);

    TString hNameO(Form("hLYPEO_%s_Config",config[ival]));
    TH2* hLYO = 0;
    hLYO = new TH2D(hNameO.Data(),hNameO.Data(),55,0,55,200,0,2000);
    vhistO.push_back(hLYO);
    
    TString hNameScalar(Form("hScalars_%s_Config",config[ival]));
    TH1* hScalar = new TH1D(hNameScalar.Data(),hNameScalar.Data(),10,0,10.0); //10 bins of size 1
    vhistScalars.push_back(hScalar);

    hLYD->SetDirectory(gROOT); // These objects are owned by ROOT (not sure what this means ???)
    hLYO->SetDirectory(gROOT);
    hScalar->SetDirectory(gROOT);
    
    
    int ievent = 0; //Filling out these 5 histograms!
    
    while(ievent < Sim_Tree -> GetEntries()){ //Sim_Tree -> GetEntries() is the number of generated events 
        
        Sim_Tree -> GetEntry(ievent); //Read entry (event) "ievent" from Sim_Tree
        int ndetpe = 0;
        hScalar->Fill(ev->vol0); //for initial volume
        
        //vol0=1 if the event was in d2o volume, vol0=2 if it was in h2o volume.
      
 
        //if event happened in d2o volume...
        if(ev->vol0==1) hLYD->Fill(ev->sourceParticleEnergy,ev->numHits);// ...fill histogram
        
        
      
        if(ev->vol0>=0) hLYO->Fill(ev->sourceParticleEnergy,ev->numHits);
        // because oxygen is in BOTH d2o and h2o. Recall that vol0=2 if event happened in h2o volume.
      
            
        ievent++; //number of events counter
        
    } // END OF While (ievent < #entries)
    
   
    MCanalysisFile->WriteObject(hLYD,"hLYD");
    ival++;  //number of configurations counter  
    
 } //END For (fnames = # of configurations)
 
 MCanalysisFile->Close();


  
}// END of compareConfig()




   








 
   



