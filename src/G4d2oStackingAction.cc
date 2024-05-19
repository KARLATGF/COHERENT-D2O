#include "G4d2oDetector.hh"
#include "G4d2oStackingAction.hh"

#include "G4Step.hh"
#include "G4OpticalPhoton.hh"
#include "G4VProcess.hh"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Constructor
G4d2oStackingAction::G4d2oStackingAction(G4d2oDetector * det) :
  fDetector(det)
{
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Destructor
G4d2oStackingAction::~G4d2oStackingAction()
{
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Prepare New event
void G4d2oStackingAction::PrepareNewEvent()
{
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// New Stage
void G4d2oStackingAction::NewStage()
{
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Public Functions
G4ClassificationOfNewTrack G4d2oStackingAction::ClassifyNewTrack(const G4Track * the_track)
{
  G4ClassificationOfNewTrack classification    = fWaiting;
  const G4VProcess * kprocess = the_track->GetCreatorProcess();

  if (the_track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() ) {
    G4double energy_eV = the_track->GetTotalEnergy() / eV;

    // Checks in WCSim StackngAction regardng process type etc...
    // Are those necessary? Or can we just kill based on energy?
    // Checks that the process type !=3 ...which I think just means not created from an optical process
    //  e.g. reflection?
    //  So probably needed so don't throw the dice multiple time for a single photon
    if (kprocess && kprocess->GetProcessType() !=3) {
      if (G4UniformRand() > fDetector->EvalPMTQE(energy_eV) ) {
        classification = fKill;
      }
    }
  }

  return classification;
} // END OF USER STEPPING ACTION
