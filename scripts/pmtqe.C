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


void pmtqe(const char* fname=""){
  std::vector<double> qe_bialkali_xnm;
  for(auto xgev : qe_bialkali_xgev) qe_bialkali_xnm.push_back(hc_evnm/(xgev*1e9));

  gBialkali = new TGraph(qe_bialkali_xnm.size());
  gBialkali_ev = new TGraph(qe_bialkali_xgev.size());

  for(int ipoint = 0; ipoint<qe_bialkali_xnm.size(); ipoint++ ){
    gBialkali->SetPoint(ipoint, qe_bialkali_xnm[ipoint],qe_bialkali_y[ipoint]);
    gBialkali_ev->SetPoint(ipoint,qe_bialkali_xgev[ipoint]*1e9,qe_bialkali_y[ipoint]);
  }

  new TCanvas();
  gBialkali->Draw("ACP");
  new TCanvas();
  gBialkali_ev->Draw("ACP");
  if(TString(fname)=="") return;
  
  TFile* fin = TFile::Open(fname);
  TTreeReader t("Sim_Tree", fin);
  TTreeReaderValue<simEvent> ev(t,"eventData");
  TH1* hPEperEvent = new TH1D("hPEperEvent","hPEperEvent",1000,0,2000);
  TH1* hPEperEventEdge = new TH1D("hPEperEventEdge","hPEperEventEdge",1000,0,2000);
  while(t.Next()){
    int ndetpe = 0;
    for(int ihit = 0; ihit < ev->numHits; ihit++){
      if(gRandom->Uniform()<qe_ev(ev->GetHit(ihit)->photonEnergy*1e6))
      {
        ndetpe++;
      }
    }
    for(int ihit = 0; ihit < ev->numHitsArea; ihit++){
      if(gRandom->Uniform()<qe_ev(ev->GetAHit(ihit)->photonEnergy*1e6))
      {
        ndetpe++;
      }
    }
    hPEperEvent->Fill(ndetpe);
    std::cout<<"TotalPEs "<<ev->numHits << " Total Detected PEs: "<<ndetpe<<std::endl;

  }
  new TCanvas();
  hPEperEvent->Draw();
}

