
#include "G4d2oElectronGun.hh"
#include "inputVariables.hh"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TSystem.h"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"

G4d2oElectronGun::G4d2oElectronGun()
{
    
    G4cout << "\tConstructing G4d2oElectronGun..." ;
    
    //Set random seed variables
    input = inputVariables::GetIVPointer();
    G4int irand = input->GetRandomStatus();
    if(irand==0) gRandom->SetSeed(1);
    if(irand==1) gRandom->SetSeed(0);
    totalEvents = input->GetNumberOfEvents();

    //Initialize a few things
    G4int n_particle = 1;
    particleGun = new G4ParticleGun(n_particle);
//    sourceEnergy = 4000.0*MeV;
    sourceEnergy = 10.0*MeV;
//    sourceEnergy = 140.0*keV;
    
    //Set the particle name
    G4String particleName = "e-";
    //G4String particleName = "mu-";

    //set up the gun
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    particleGun->SetParticleDefinition(particleTable->FindParticle(particleName));
    
    SourcePosition.setX(0);
    SourcePosition.setY(0);
    SourcePosition.setZ(0);
    
    G4double thetaAngle1 = 0.0; //in degrees
    G4double thetaAngle2 = 180.0; //in degrees
    
    cthetaRange1 = cos(thetaAngle1*deg);
    cthetaRange2 = cos(thetaAngle2*deg);
    
    G4d2oNeutrinoAlley *detCon = (G4d2oNeutrinoAlley*)G4RunManager::GetRunManager()->GetUserDetectorConstruction();
    G4d2oDetector *theDet = (G4d2oDetector*)detCon->GetDetectorPtr();
    // Inner Tank of D20
    tankSize.set(theDet->TankX(),theDet->TankY(),theDet->TankZ());
    // Outer Tank of H20
    //tankSize.set(theDet->OuterTankX(),theDet->OuterTankY(),theDet->OuterTankZ());

    etfTimer = new ReplayTools();
    etfTimer->PrepareETFTimer(5000, input->GetNumberOfEvents()); //time in ms
    
    G4cout << "done." << G4endl;

}//END of constructor

G4d2oElectronGun::~G4d2oElectronGun()
{
    
	G4cout<<"Deleting G4d2oElectronGun...";
	
//    delete particleGun;
	
	G4cout<<"done."<<G4endl;
	
}//END of destructor

void G4d2oElectronGun::GeneratePrimaries(G4Event* anEvent)
{
    
    theEventNum = anEvent->GetEventID();
    etfTimer->SetCurrentEvent(theEventNum);
    gSystem->ProcessEvents();
    
    if(anEvent->GetEventID()==0 && input->GetPrintStatus()!=3) etfTimer->StartUpdateTimer();
        
    G4double px, py, pz;
    G4double ctheta, stheta, phi;

    sourceEnergy = (30.0000)*MeV;
    //sourceEnergy = (0.0+G4UniformRand()*53.0)*MeV;
    //sourceEnergy = (29.5+G4UniformRand())*MeV; //Energies between 29.5 and 30.5 for first study.
  
    particleGun->SetParticleEnergy( sourceEnergy );
  
    G4Navigator* Navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

    int ntries = 0;
    while(true){
      SourcePosition.set(tankSize.x()*(G4UniformRand()-0.5),
                         tankSize.y()*(G4UniformRand()-0.5),
                         tankSize.z()*(G4UniformRand()-0.5));
      
      G4VPhysicalVolume* volume = Navigator->LocateGlobalPointAndSetup(SourcePosition);
//       if(volume->GetName() == "h2oPhysV"
//          || volume->GetName() == "d2oPhysV"
//          || volume->GetName() == "acrylicPhysV"
//          )
            if(volume->GetName() == "d2oPhysV")
//        if(volume->GetName() == "d2oPhysV")
        break;
      ntries++;
    }
  
    particleGun->SetParticlePosition( SourcePosition );

    ctheta = G4UniformRand()*(cthetaRange1-cthetaRange2) + cthetaRange2;
    phi = 2.0*TMath::Pi()*G4UniformRand();
    
    stheta = sqrt( 1 - pow(ctheta,2) );
    
    pz = ctheta;
    py = stheta*cos(phi);
    px = stheta*sin(phi);
    
    
////    ///// temp for fixed initial direction /////
//    px = 0.0;
//    py = -1.0;
//    pz = 0.0;
    //particleGun->SetParticlePosition(G4ThreeVector(-5.0*cm,0.0,92.5*cm+3.8*mm));
    //particleGun->SetParticlePosition(G4ThreeVector(0.0*cm,200.0*cm,0.0*cm));
////    ///// end temp /////
    
    initDir.set(px,py,pz);
    
    particleGun->SetParticleMomentumDirection(initDir);
    
    particleGun->GeneratePrimaryVertex(anEvent);
    
}//END of GeneratePrimaries()

