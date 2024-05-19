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

//gSystem->AddIncludePath(" -I/Users/frn/Documents/Projects/COHERENT/d2o/sim/G4d2o/simEvent  ")

const double hc_evnm = 1.23984193 *1e3;

std::vector<double> qe_bialkali_y = {0.000,0.000,0.000,0.000,
  0.001,0.003,0.008,0.017,
  0.033,0.052,0.076,0.106,
  0.139,0.164,0.193,0.204,
  0.228,0.240,0.254,0.254,
  0.240,0.213,0.164,0.032,
  0.003,0.0};

std::vector<double> qe_bialkali_xgev = {1.63e-9,1.68e-9,1.72e-9,1.77e-9,
  1.82e-9,1.88e-9,1.94e-9,2.00e-9,
  2.07e-9,2.14e-9,2.21e-9,2.30e-9,
  2.38e-9,2.48e-9,2.58e-9,2.70e-9,
  2.82e-9,2.95e-9,3.10e-9,3.26e-9,
  3.44e-9,3.65e-9,3.88e-9,4.13e-9,
  4.43e-9,4.7e-9};

TGraph* gBialkali = 0;
TGraph* gBialkali_ev = 0;

double qe_nm(double wavelength_nm){
  if(!gBialkali) return 0.0;
  if(wavelength_nm>gBialkali->GetX()[0]
     || wavelength_nm<gBialkali->GetX()[gBialkali->GetN()-1])
    return 0.0;
  return gBialkali->Eval(wavelength_nm,0,"");
}

double qe_ev(double energy_ev){
  if(!gBialkali_ev) return 0.0;
  if(energy_ev<gBialkali_ev->GetX()[0]
     || energy_ev>gBialkali_ev->GetX()[gBialkali_ev->GetN()-1])
    return 0.0;
  return gBialkali_ev->Eval(energy_ev,0,"");
}



void plotLightYieldChain(){
  std::vector<double> qe_bialkali_xnm;
  for(auto xgev : qe_bialkali_xgev) qe_bialkali_xnm.push_back(hc_evnm/(xgev*1e9));
  
  gBialkali = new TGraph(qe_bialkali_xnm.size());
  gBialkali_ev = new TGraph(qe_bialkali_xgev.size());
  
  for(int ipoint = 0; ipoint<qe_bialkali_xnm.size(); ipoint++ ){
    gBialkali->SetPoint(ipoint, qe_bialkali_xnm[ipoint],qe_bialkali_y[ipoint]);
    gBialkali_ev->SetPoint(ipoint,qe_bialkali_xgev[ipoint]*1e9,qe_bialkali_y[ipoint]);
  }

  std::vector<int> rlist = {24};
  auto tch = new TChain("Sim_Tree","Sim_Tree");
  for(auto rn : rlist)
    tch->Add(Form("data/test/Sim_D2ODetector%03d.root",rn));
  TTreeReader t(tch);
  TTreeReaderValue<simEvent> ev(t,"eventData");
  
  TH3* hPEperEvent = new TH3D("hPEperEvent","hPEperEvent",200,-1.0,1.0,12,0,60,1000,0,2000);

  TProfile2D* hPEperEventProf = new TProfile2D("hPEperEventProf","hPEperEventProf",200,-1.0,1.0,12,0,60,"S");

  TProfile3D* hPEperEventPos = new TProfile3D("hPEperEventPos","hPEperEventPos",140,-700.0,700.0,60,-300,300,140,-700.0,700.0,"S");

  
  while(t.Next()){
    int ndetpe = 0;
    for(int ihit = 0; ihit < ev->numHits; ihit++){
      ndetpe++;
    }
    for(int ihit = 0; ihit < ev->numHitsArea; ihit++){
      ndetpe++;
    }
    hPEperEvent->Fill(ev->direction0.Y(),ev->sourceParticleEnergy,ndetpe);
    hPEperEventProf->Fill(ev->direction0.Y(),ev->sourceParticleEnergy,ndetpe);
    hPEperEventPos->Fill(ev->position0.X(),ev->position0.Y(),ev->position0.Z(),
                         ndetpe);
    //std::cout<<"TotalPEs "<<ev->numHits << " Total Detected PEs: "<<ndetpe<<std::endl;
    
    //if (hPEperEvent->GetEntries()>1000) break;
  }
  
}
