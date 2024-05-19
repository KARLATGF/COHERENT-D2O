
Double_t gaus1d(Double_t *x, Double_t *par){
 
    return par[0]*TMath::Gaus(x[0],par[1],par[2]);
    
}

void fitCylindricalLightCollection(){
 
//    TString theSingleSidePMTFile("/Volumes/LaCie/coherent/simulatedData/D2O/cylindrical/singleSidePMTs/Sim_D2ODetector001.root");
//    TString theDoubleSidePMTFile("/Volumes/LaCie/coherent/simulatedData/D2O/cylindrical/doubleSidePMTs/Sim_D2ODetector001.root");
//    TString theSingleSidePMTFile("/Volumes/LaCie/coherent/simulatedData/D2O/cylindrical/singleSidePMTs/Sim_D2ODetector100.root");
//    TString theDoubleSidePMTFile("/Volumes/LaCie/coherent/simulatedData/D2O/cylindrical/doubleSidePMTs/Sim_D2ODetector100.root");
    TString theSingleSidePMTFile("/Volumes/LaCie/coherent/simulatedData/D2O/cylindrical/singleSidePMTs/Sim_D2ODetector200.root");
    TString theDoubleSidePMTFile("/Volumes/LaCie/coherent/simulatedData/D2O/cylindrical/doubleSidePMTs/Sim_D2ODetector200.root");

    TFile *fSingleSide = TFile::Open(theSingleSidePMTFile);
    TFile *fDoubleSide = TFile::Open(theDoubleSidePMTFile);

    TH2 *hSingleSide = (TH2*)fSingleSide->Get("histos/photonsPerEvent");
    TH2 *hDoubleSide = (TH2*)fDoubleSide->Get("histos/photonsPerEvent");

    hSingleSide->SetLineColor(kBlack);
    hDoubleSide->SetLineColor(kRed);
    
    hSingleSide->SetAxisRange(0,2000,"x");
    
    gStyle->SetOptStat(0);
    TCanvas *c1 = (TCanvas*)PlotTools::canvas("c1",1200,1200);
    c1->Divide(1,2);
    
    c1->cd(1);
    hSingleSide->DrawCopy();
    hDoubleSide->DrawCopy("same");

    c1->cd(2);
    
    Double_t maxSingleX = hSingleSide->GetBinCenter(hSingleSide->GetMaximumBin());
    TF1 *fGausSingle = new TF1("fGausSingle",gaus1d,maxSingleX-100, maxSingleX+300,3);
    fGausSingle->SetLineColor(kBlack);
    fGausSingle->SetNpx(1000);
    fGausSingle->SetParameters(hSingleSide->GetMaximum(),maxSingleX,50);
    hSingleSide->Fit("fGausSingle","r");

    Double_t maxDoubleX = hDoubleSide->GetBinCenter(hDoubleSide->GetMaximumBin());
    TF1 *fGausDouble = new TF1("fGausDouble",gaus1d,maxDoubleX-100, maxDoubleX+300,3);
    fGausDouble->SetLineColor(kRed);
    fGausDouble->SetNpx(1000);
    fGausDouble->SetParameters(hDoubleSide->GetMaximum(),maxDoubleX,50);
    hDoubleSide->Fit("fGausDouble","r");

    hSingleSide->DrawCopy();
    hDoubleSide->DrawCopy("same");

}
