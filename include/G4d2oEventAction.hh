#ifndef G4d2oEventAction_hh
#define G4d2oEventAction_hh 1

#include "G4d2oDetectorHit.hh"
#include "G4d2oPrimaryGeneratorAction.hh"
#include "inputVariables.hh"
#include "simEvent.h"

#include "TH2D.h"
#include "TVector3.h"

#include "G4UserEventAction.hh"
#include "G4Run.hh"
#include "G4Event.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "TGraph.h"

class G4d2oEventAction : public G4UserEventAction
{
protected:
    
    G4HCofThisEvent *HCoE;
    G4d2oDetectorHitsCollection *detHC[MAX_SEN_DET];
    G4d2oDetectorHitsCollection *sourceHC[MAX_SEN_DET];
    G4d2oDetectorHitsCollection *targHC[MAX_SEN_DET];
    
    G4d2oPrimaryGeneratorAction *thePGA;
    G4d2oDetector *theDet;
    
    G4int numEvents;
    
    TFile *fout;
    TTree *outTree;
    TH2D *hPMTArray;
    TH2D *hPMTNumVsTime;
    TH1D *hPhotonEnergy;
    TH2D *hSidePMT[2];
    TH1D *hTotalPhotons;
  
    TGraph* gBialkali_ev = 0;
  
    inputVariables *input;
    G4int numPMTs, numRows, numInEachRow;
    
    simEvent *theEventData;
    
    Int_t eventNumber;
  
    double qe_ev(double energy_ev);
  
    void GetHitsCollection( void );
    virtual void ProcessEvent( void );
    virtual void ZeroEventVariables( void );
        
public:
    
    G4d2oEventAction();
    ~G4d2oEventAction();
    void BeginOfEventAction( const G4Event *thisEvent);
    void EndOfEventAction( const G4Event *thisEvent);
    
    TFile * OutFilePtr() {return fout;}
    
}; //END of class G4d2oEventAction

#endif
