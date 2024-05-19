#ifndef G4d2oDetector_H
#define G4d2oDetector_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4RotationMatrix.hh"

#include "TH2F.h"

#include "G4d2oMaterialsDefinition.hh"

#include "inputVariables.hh"

class G4d2oDetector
{
public:
    
  G4d2oDetector();
  virtual ~G4d2oDetector();
    
  virtual G4LogicalVolume * GetDetector();
  //virtual G4ThreeVector GetPMTPosition(Int_t iPMT);
//   virtual G4ThreeVector GetTranslation();
//   virtual G4RotationMatrix *GetRotation();

  inline TH1D * GetPMTPositionHistogram(Int_t iDim) {return hPanel[iDim];}
  inline G4int GetTotalPMTS() { return totPMT; }
  inline G4int GetPMTRows() { return numRows; }
  inline G4int GetPMTsInRow() { return numInRows; }
  
  inline G4double TankLength() {return d2oLength;} //x
  inline G4double TankWidth() {return d2oWidth;} //y
  inline G4double TankHeight() {return d2oHeight;} //z
  inline G4double TankX() {return d2oLength;} //x
  inline G4double TankY() {return d2oWidth;} //y
  inline G4double TankZ() {return d2oHeight;} //z
  inline G4double TankRadius() {return d2oRadius;} //z

  inline G4double OuterTankLength() {return h2oInnerLength;} //x
  inline G4double OuterTankWidth() {return h2oInnerWidth;} //y
  inline G4double OuterTankHeight() {return h2oInnerHeight;} //z
  
  inline G4double OuterTankX() {return h2oInnerLength;} //x
  inline G4double OuterTankY() {return h2oInnerWidth;} //y
  inline G4double OuterTankZ() {return h2oInnerHeight;} //z

  inline G4double OuterReflectorX() {return h2oInnerLength;} //x
  inline G4double OuterReflectorY() {return h2oInnerWidth;} //y
  inline G4double OuterReflectorZ() {return h2oInnerHeight;} //z
    
  inline G4double OuterContX() {return contOuterLength;} //x
  inline G4double OuterContY() {return contOuterWidth;} //y
  inline G4double OuterContZ() {return contOuterHeight;} //z
  
  inline G4double OuterVetoX() {return vetoOuterLength;} //x
  inline G4double OuterVetoY() {return vetoOuterWidth;} //y
  inline G4double OuterVetoZ() {return vetoOuterHeight;} //z

  inline G4double SidePMTX() {return sideLiningX;}
  inline G4double SidePMTZ() {return sideLiningZ;}
  
private:
 //virtual G4LogicalVolume *GetTeflon();
// virtual G4LogicalVolume *GetSphericalPMT(); //New
  virtual G4LogicalVolume *GetEllipsoidPMT(); //*New
  //virtual G4LogicalVolume *GetAreaPMT(); //For spherical PMTs
  virtual void DetermineSpacing();
  virtual void Initialize();
  virtual G4LogicalVolume *GetTotalDetectorLogV();
  virtual G4LogicalVolume *GetOuterVesselLogV();
  //virtual G4LogicalVolume *GetTeflonLiningLogV();
  virtual G4LogicalVolume *GetH2OLogV();
  virtual G4LogicalVolume *GetAcrylicLogV();
  virtual G4LogicalVolume *GetD2OLogV();
  virtual void PlacePMTs(G4LogicalVolume *thePMTLogV, G4LogicalVolume *theMotherLogV);
  

  
  virtual G4LogicalVolume *GetTeflonLiningLLogV();
  virtual G4LogicalVolume *GetTeflonLiningFLogV();
  virtual G4LogicalVolume *GetTeflonLiningTLogV();
  


  //Dimensions
  const G4double in;
  
  inputVariables *input;
  G4d2oMaterialsDefinition *matPtr;
  
  TH1D *hPanel[3];
  G4int totPMT;
  G4int numRows;
  G4int numInRows;
  
  G4double d2oRadius;
  G4double d2oLength; //x
  G4double d2oWidth; //y
  G4double d2oHeight; //z
  
  G4double acrylicThickness;
  G4double pmtDepth;
  G4double teflonThickness;
  G4double teflonReflectivity;
  G4double teflonSigmaAlpha; //roughness parameter. 0=polished
  G4double pmtDiameter; //controlled from beamOn.dat or command line now
  
  G4double pmtMinorAxis; // Minor axis for ellipsoidal pmts *New
  G4double pmtMajorAxis; // Major axis for ellipsoidal pmts *New
  G4double pmtLegLength; //*New

  G4double muMetalThickness;
  G4double pmtWindowThickness;  
  G4double minGapBetweenPMTs;
  G4double outerContainerThickness;
  
  G4double h2oLength; //x
  G4double h2oWidth; //y
  G4double h2oHeight; //z
  
  G4double pmtReflectorLength; //x
  G4double pmtReflectorWidth; //y
  G4double pmtReflectorHeight; //z
  G4double pmtReflectorRadius;
  
  G4double h2oInnerLength; //x
  G4double h2oInnerWidth; //y
  G4double h2oInnerHeight; //z
  
  G4double h2oInnerRadius;
  G4double TCRadius; //*New!!
  
  G4double sideLiningX, sideLiningZ;
  G4double theSideLiningThickness;
  
  G4double contOuterLength; //x
  G4double contOuterWidth; //y
  G4double contOuterHeight; //z
  
  G4double muonVetoThickness;
  G4int muonVetoLayers;
  
  G4double vetoOuterLength; //x
  G4double vetoOuterWidth; //y
  G4double vetoOuterHeight; //z
 
  G4int numInX, numInY, numInZ;
  G4double spacingInX, spacingInY, spacingInZ;
  
  G4int numInCircle;
  G4double spacingInCircle;
  G4double totalAngle; //* New!!
  
  G4double teflonHeight; //!!new
  G4double teflonLength; //!!new
  G4double teflonOuterR; //!!new
  G4double teflonInnerR; //!!new
  
  G4double PMTcoatingReflectivity; //!!new
  G4double PMTcoatingSigmaAlpha; //!!new

  G4bool bUseBottomPMTs;
  
protected:

};//END of class G4d2oDetector

#endif
