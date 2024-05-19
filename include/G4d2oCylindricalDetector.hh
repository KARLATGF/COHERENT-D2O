#ifndef G4d2oCylindricalDetector_H
#define G4d2oCylindricalDetector_H 1

#include "G4d2oDetector.hh"

class G4d2oCylindricalDetector : public G4d2oDetector
{
public:
    
    G4d2oCylindricalDetector();
    ~G4d2oCylindricalDetector();
    
    virtual G4LogicalVolume * GetDetector();
    
protected:
    virtual void Initialize();
    virtual void DetermineSpacing();
    
    G4LogicalVolume *GetTotalDetectorLogV();
    G4LogicalVolume *GetShieldingLogV();
    G4LogicalVolume *GetOuterVesselLogV();
    G4LogicalVolume *GetTeflonLiningLogV();
    G4LogicalVolume *GetTeflonLiningCapLogV(const char * name, bool boundPMTs);
    G4LogicalVolume *GetH2OLogV();
    G4LogicalVolume *GetAcrylicLogV();
    G4LogicalVolume *GetD2OLogV();
    G4LogicalVolume *GetbubbleLogV();
    G4LogicalVolume *GetVetoLogV();
    G4LogicalVolume *GetInnerVetoLogV();
    
    //-----------------------------------------------------------
    
    G4LogicalVolume *GetTeflonLiningLLogV();
    G4LogicalVolume *GetTeflonLiningFLogV();
    G4LogicalVolume *GetTeflonLiningTLogV();
    
    //-----------------------------------------------------------

    void PlacePMTs(G4LogicalVolume *thePMTLogV, G4LogicalVolume *theMotherLogV);

    G4int iUseBottomPMTs;
    
    G4double acrylicEndCapThickness;
    
    //-----------------------------------------------------------
    
    G4double h2oInnerRadius;
    G4double TCRadius; 
    
    G4double teflonHeight; 
    G4double teflonLength; 
    G4double teflonOuterR; 
    G4double teflonInnerR; 
    G4int numInCircle;
    G4double spacingInCircle;
    G4double totalAngle; 
  
    //-----------------------------------------------------------
    
};//END of class G4d2oCylindricalDetector

#endif
