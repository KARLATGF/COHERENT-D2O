
void openLastFile(){

    ifstream infile("scripts/lastFileName.txt");
    TString fname;
    infile >> fname;
    infile.close();
    
    TFile *f = TFile::Open(fname.Data(),"read");
    TBrowser *b = new TBrowser();    

}
