void plotPMMA(){
    
    Double_t eV=1.0,m=100.0,mm=0.1;
    
    Double_t hc_evnm = 1.23984193 *1e3;
    
    Double_t photonEnergy[] =
            {2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
            2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
            2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
            2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
            2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
            2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
            2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
            3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
            3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
            3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};
    
    
    Int_t nentries = sizeof(photonEnergy)/sizeof(Double_t);
    Double_t* photonWavelength = new Double_t[nentries];
    for(int i = 0; i<nentries; i++){
        photonWavelength[i] = hc_evnm/photonEnergy[i];
    }
    
    Double_t absWLSfiber[] =
       {5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
            5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,
            5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,5.40*m,1.10*m,
            1.10*m,1.10*m,1.10*m,1.10*m,1.10*m,1.10*m, 1.*mm, 1.*mm, 1.*mm, 1.*mm,
             1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm, 1.*mm};
    std::cout<<nentries<<std::endl;
    TGraph* gr = new TGraph(nentries,photonWavelength,absWLSfiber);
    gr->SetMarkerStyle(21);
    gr->Draw("AL");
}
