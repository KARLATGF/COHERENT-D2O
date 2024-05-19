void plotPMTPositions(const char* fname)
{
    TFile* fin = TFile::Open(fname);
    auto fd = (TDirectoryFile*)fin->Get("pmtPositions");
    auto pmtPosX = (TH1*)fd->Get("pmtPosX");
    auto pmtPosY = (TH1*)fd->Get("pmtPosY");
    auto pmtPosZ = (TH1*)fd->Get("pmtPosZ");
    
    for(int ibin = 1; ibin<=pmtPosX->GetNbinsX(); ibin++)
    {
        std::cout<<ibin<<" ("<<pmtPosX->GetBinContent(ibin)<<","
        <<pmtPosY->GetBinContent(ibin)<<","
        <<pmtPosZ->GetBinContent(ibin)<<")"<<std::endl;
    }
}
