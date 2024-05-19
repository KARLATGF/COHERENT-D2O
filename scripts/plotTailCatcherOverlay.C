void plotTailCatcherOverlay(){
  std::vector<std::string> fnames = {"hPEperEvent_5cmLW.root",
  "hPEperEvent_10cmLW.root",
  "hPEperEvent_15cmLW.root",
    "hPEperEvent_20cmLW.root"};
  
  std::vector<int> tcval = {5, 10, 15, 20};
  
  std::vector<TFile*> vfin;
  std::vector<TH1*> vhist;
  int ival = 0;
  for(auto fname : fnames){
    vfin.push_back(TFile::Open(fname.c_str()));
    vfin.back()->Print();
    vhist.push_back((TH1*)vfin.back()->Get(Form("hPEperEvent_%dcmLW",tcval[ival])));
    vhist.back()->SetLineColor(kBlack+ival);
    vhist.back()->SetLineWidth(2);
    vhist.back()->Rebin(5);
    vhist.back()->SetTitle(Form("LW %d cm",tcval[ival]));
    vhist.back()->Draw(ival==0?"":"same");
    ival++;
    
  }
  
}
