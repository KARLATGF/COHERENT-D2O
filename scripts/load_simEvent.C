void load_simEvent(){
  gSystem->AddIncludePath("-IsimEvent");
  gSystem->Load("G4d2o-build/simEvent/libsimEvent.so");
}
