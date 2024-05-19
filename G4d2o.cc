
#include <stdlib.h>
#include <fstream>
#include <sys/time.h>

//geant controls
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "G4d2oPhysicsList.hh"
#include "G4d2oRunAction.hh"
#include "G4d2oNeutrinoAlley.hh"

#include "G4d2oEventAction.hh"
#include "G4d2oStackingAction.hh"

#include "G4d2oPhotonGun.hh"
#include "G4d2oElectronGun.hh"
#include "G4d2oMichelElectronGun.hh"
#include "G4d2oCosmicGun.hh"

#include "inputVariables.hh"

#include "TRandom.h"

int main(int argc, char** argv)
{
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    //Read in some input variables
    inputVariables *input = inputVariables::GetIVPointer(argc, argv);
    
    G4int numEvents = input->GetNumberOfEvents();
    G4int ivis = input->GetVisualization();
    G4int irand = input->GetRandomStatus();
    G4int runno = input->GetRunNumber();
    G4int iPhysics = input->GetPhysicsType();
    G4int iNeutronHP = input->GetNeutronHP();
    G4int iPGA = input->GetPGAType();
    G4long userSEED = input->GetUserSEED();
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Construct the default run manager
    G4RunManager *runManager = new G4RunManager;
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    //Set the Random Seed
    timeval theTime;
    gettimeofday(&theTime, NULL);
    G4long iseed = (theTime.tv_sec *1000) + (theTime.tv_usec / 1000);
        
    if(irand==0){ //fixed random seed for debugging
        iseed = 20336446;
        iseed = 1531567817239; //Crash in early event
        iseed = 1533952177170; //Large cosmic ray
        if(userSEED>0){
            iseed = userSEED;
            G4cout<<"Setting seed to "<<iseed<<G4endl;
        }
        
        gRandom->SetSeed(1);
    }
    if(irand==1){
        gRandom->SetSeed(0);
    }

    CLHEP::HepRandom::setTheSeed(iseed);

    //keep a log of the starting random seed in case we want to look at an event again
    FILE *outSeed;
    outSeed = fopen("seedHistory.dat","a");
    fprintf(outSeed,"Run %03d: %ld\n",runno,iseed);
    fclose(outSeed);

    //////////////////////////////////////////////////////////////////////////////////////////////
    //Detector construction
    G4cerr<<"Constructing detector..."<<G4endl;
    G4d2oNeutrinoAlley *detCon = new G4d2oNeutrinoAlley();
    runManager->SetUserInitialization(detCon);
    
    //////////////////////////////////////////////////////////////////////////////////////////////
	//Physics list
    G4cerr<<"Constructing physics..."<<G4endl;
    G4bool bPhysics = false; //true=turn on physics, false=no physics (for geometry debugging)
    if(iPhysics>0) bPhysics = true;
    G4bool bNeutronHP = false; //true=HP neutron physics, false=non-HP neutron physics
    if(iNeutronHP==1) bNeutronHP = true;
    
    runManager->SetUserInitialization(new G4d2oPhysicsList(bPhysics,bNeutronHP));

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Primary Generator Action 
    G4cerr<<"Constructing primary generator action..."<<G4endl;
    if(iPGA==0) runManager->SetUserAction(new G4d2oPhotonGun);
    else if(iPGA==1) runManager->SetUserAction(new G4d2oElectronGun);
    else if(iPGA==2) runManager->SetUserAction(new G4d2oCosmicGun);
    else if(iPGA==3) runManager->SetUserAction(new G4d2oMichelElectronGun);

    std::cout<<"TopOfCeiling "<<detCon->GetTopOfConcreteCeiling()<<std::endl;
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Run and Event Action 
    G4cerr<<"Constructing run action..."<<G4endl;
    runManager->SetUserAction(new G4d2oRunAction);
    G4cerr<<"Constructing event action..."<<G4endl;
    runManager->SetUserAction(new G4d2oEventAction);
    G4cout<<"Constructing tracking action..."<<G4endl;
    runManager->SetUserAction(new G4d2oStackingAction(detCon->GetDetectorPtr() ) );

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Initialize G4 kernel
    G4cerr<<"Initializing run..."<<G4endl;
    runManager->Initialize();
    
    //Get visualization ready
    G4VisManager *visManager = 0; 
//#ifdef G4VIS_USE
    if(ivis>0){
        visManager = new G4VisExecutive;
        visManager->Initialize();
    }
//#endif

    // get the pointer to the UI manager 
    G4UImanager *UI = G4UImanager::GetUIpointer();
   
    // execute visualisation macro
    G4String visName;
    if(ivis==0) visName = "/control/execute mac/vis-novis.mac";
    if(ivis==1) visName = "/control/execute mac/vis-opengl.mac";
    if(ivis==2) visName = "/control/execute mac/vis-wired.mac";
    if(ivis==3) visName = "/control/execute mac/vis-dawn.mac";
    
    UI->ApplyCommand(visName);
    //UI->ApplyCommand("/tracking/verbose 1");
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // turn beam on
    runManager->BeamOn(numEvents);

    // job termination
//#ifdef G4VIS_USE
    delete visManager;
//#endif
    
    delete runManager;
    
    return 0;
    
}//END of main()
