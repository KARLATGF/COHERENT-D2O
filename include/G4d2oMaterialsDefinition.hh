#ifndef G4d2oMaterialsDefinition_H
#define G4d2oMaterialsDefinition_H 1

#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VPhysicalVolume.hh"

enum materialName
{
    AIR, ALUMINUM, LEAD,  
    POLY, STEEL,
    PMMA, MUMETAL,
    COPPER,
    PLASTIC,
    VACUUM, H2O, D2O,
    FUSEDSILICA, BOROSILICATE,
    PHOTOCATHODE,
    TEFLON,
    CONCRETE
};//END of enum materialName 

class G4d2oMaterialsDefinition
{
public:
    
    G4d2oMaterialsDefinition();
    ~G4d2oMaterialsDefinition();
    
    G4Material * GetMaterial( materialName matName );
    void SetReflector(G4VPhysicalVolume *theExitingVolume, G4VPhysicalVolume *theEnteringVolume,
                            G4double theReflectivity, G4double theSigmaAlpha=0.0);
    
protected:
    
private:
    
    G4NistManager *manager;
    
    G4Material *matAir, *matAl, *matSteel;
    G4Material *matPb, *matPoly;
	G4Material *matPMMA;
	G4Material *matMuMetal;
    G4Material *matCopper;
    G4Material *matPlastic;
    G4Material *matVacuum;
    G4Material *matH2O, *matD2O;
    G4Material *matFusedSilica;
    G4Material *matBorosilicate;
    G4Material *matPhotoCathode;
    G4Material *matTeflon;
    G4Material *matConcrete;

//    G4double temperature;
//    G4double pressure;
    
    void SetUniformOpticalProperties(G4Material *theMat, G4double theIndex, G4double theAbsLength);
    void SetOpticalProperties(G4Material *theMat, G4double theIndex, G4String absLengthFile);
    
};//END of class G4d2oMaterialsDefinition

#endif
