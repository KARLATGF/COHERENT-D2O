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


void plotCosmic()
{
    
    double pePerMeV = 20.2; //for 10cm catcher
    double simDuration = 32.0;
    double SNSHoursPerYear = 5000.0;
    
    
    double totalLiveTimeSNSYear = SNSHoursPerYear * 60 * 3600 * 1e-5;
    
    
    std::string fname = "data/CosmicMicronPhC/Sim_D2ODetector001.root";
    auto Sim_Tree = new TChain("Sim_Tree","Sim_Tree");
    Sim_Tree->Add(fname.c_str());
    Sim_Tree->Print();
    TTreeReader t(Sim_Tree);
    TTreeReaderValue<simEvent> ev(t,"eventData");
    
    TString hNameLY(Form("hLY_Cosmic"));
    TH1* hLY = new TH1D(hNameLY.Data(),hNameLY.Data(),100,0,100);

    int maxevents = 1000000;
    int ievent = 0;
    
    while(t.Next()&&ievent<maxevents){
        ievent++;
        //if( !(40.0<sourceParticleEnergy&&sourceParticleEnergy<50) ) continue;
        int ndetpe = 0;
        hLY->Fill(ev->numHits/pePerMeV);
    }


    
    hLY->Scale(totalLiveTimeSNSYear/simDuration);
    hLY->Scale(1.0,"width");
    double maxbinx = 59.9;
    int maxbin = hLY->GetXaxis()->FindBin(maxbinx);
    auto projLY_Integral = (TH1*)hLY->Clone(Form("%s_integral",hLY->GetName()));
    projLY_Integral->Reset();
    for(int ibin = maxbinx; ibin>=1; ibin--){
        projLY_Integral->SetBinContent(ibin,
                                      hLY->Integral(ibin,maxbinx,"width"));
        //std::cout<<ibin<<" "<<projD->Integral(ibin,projD->GetNbinsX(),"width")<<std::endl;
    }
    hLY->Draw();
    projLY_Integral->Draw("same");
}
