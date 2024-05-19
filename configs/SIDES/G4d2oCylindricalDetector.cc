//Water Heater (Tank) configuration with PMTs on the SIDES only
// and Water (10 cm) Tail Catcher.
// added ELLIPSOIDAL PMTs from Matthew branch (G4d2o-addCylindricalDet)
// upgraded with branch from June 30th and July 7th.
//This accomodates 48 PMTs

//NEW: Teflon surounds acrylic vessel + holes for PMTs

//******************************************************************************
// Tyvek here is built by 6 different pannels. Sim runs but doesn't produce the 
// correct amount of light.
// NEW: Added top and bottom tyvek panels.
// NEW: PMT Housings added THEN REMOVED
// NEW: NOT USING G4SUBTRACTION FOR LEFT/RIGHT PANELS ANYMORE
// NEW: NEW PMTs WITH SILVER COATING TO PREVENT LIGHT FROM ESCAPING
//******************************************************************************


#include "G4d2oCylindricalDetector.hh"
#include "G4d2oRunAction.hh"
#include "G4d2oSensitiveDetector.hh"

#include "TMath.h"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4PVParameterised.hh"
#include "G4SubtractionSolid.hh"

#include "G4UnionSolid.hh" 
#include "G4IntersectionSolid.hh" 
#include "G4MultiUnion.hh" 
#include "G4OpticalSurface.hh" 
#include "G4LogicalBorderSurface.hh" 
#include "G4Ellipsoid.hh" 

using namespace std;

G4d2oCylindricalDetector::G4d2oCylindricalDetector()
{
    G4cerr << "\n\tUsing vertical cylindrical geometry with PMTs on the sides" << G4endl ;
    
    Initialize();
}

void G4d2oCylindricalDetector::Initialize(){
    
    //relevant parameters
    
    //D2O
    d2oLength = 70.0*cm;//radius
    d2oWidth = d2oLength; //radius
    //d2oHeight = 140.0*cm;//z
    d2oHeight = 161.544*cm;//z
    
    //acrylic
    acrylicThickness = 0.25*in;
    acrylicEndCapThickness = 1.0*in;
    
    //pmt
    //pmtDepth = 10.0*in;
    pmtDepth = 21.55*cm; //**New
    muMetalThickness = 0.03*in;
    pmtWindowThickness = 3.0*mm;  //making this up
    minGapBetweenPMTs = 3.5*cm;
    
    pmtDiameter = input->GetPMTDiameter()*in; //pmtDiameter controlled from beamOn.dat or command line
    pmtMinorAxis = pmtMajorAxis = pmtDiameter;
    
    pmtLegLength = 62.*mm;
    if (pmtDiameter < 0.0) {
      pmtMinorAxis = 76.72*mm * 2.0;
      pmtMajorAxis = 103.6*mm * 2.0;
    }

    //teflon
    //teflonThickness = 0.25*in;
    teflonThickness = 0.2*cm;
    teflonReflectivity = input->GetReflectivity(); //need to verify
    teflonSigmaAlpha = 0.1; //roughness parameter. 0.1 for teflon, use 0 for polished
    
    //H2O
    h2oTCthickness = input->GetTailCatcherThick()*cm;
    h2oLength = d2oLength + 2*(acrylicThickness + h2oTCthickness + pmtMinorAxis +pmtLegLength + 1.5*cm); //x
    h2oWidth = d2oLength + 2*(acrylicThickness + h2oTCthickness); //y
    
    if (pmtDiameter > 0.0) {
      h2oHeight = d2oHeight + 2*(acrylicThickness + h2oTCthickness); //z
    } else {
      h2oHeight = d2oHeight + 2*(acrylicThickness + h2oTCthickness); //z
    }
    
    //not using these for now, but set them anyway so that other parts of code don't have problems
    //These are used to throw random positions for the electron gun to include the outer h2o.
    //The vol0 variable in the event output selects 1 for d2o and 2 for h2o
    h2oInnerLength = h2oLength;
    h2oInnerWidth = h2oWidth;
    h2oInnerHeight = h2oHeight;
    
    //-----------------------------------------------------------
    h2oInnerRadius = d2oWidth/2.0 + acrylicThickness + h2oTCthickness + pmtMinorAxis/2.0;
    TCRadius = d2oWidth/2.0 + acrylicThickness + h2oTCthickness;
    //-----------------------------------------------------------
    
    //outer vessel
    outerContainerThickness = 0.25*in;
    
    ///// Full detector volume /////
    contOuterLength = h2oWidth + 2*outerContainerThickness; //x
    contOuterWidth = h2oLength + 2*outerContainerThickness; //y
    contOuterHeight = h2oHeight + 2*outerContainerThickness; //z

    ///// Shielding volume /////
    shieldThickness = input->GetShieldThickness()*in;
    shieldLength = contOuterLength + 2.0*shieldThickness; // x
    shieldWidth = shieldLength;
    shieldHeight = contOuterHeight + 2.0*shieldThickness;
    if (!input->GetBottomShielding() ) {
      shieldHeight = contOuterHeight + shieldThickness;
    }

    //Veto Layers
    muonVetoThickness = 2.54*cm;
    muonVetoLayers = 2;
    
    vetoOuterLength = shieldLength + 2.0*muonVetoLayers*muonVetoThickness; //x
    vetoOuterWidth = shieldWidth + 2.0*muonVetoLayers*muonVetoThickness; //y
    vetoOuterHeight = shieldHeight + 2.0*muonVetoLayers*muonVetoThickness; //z
    if (!input->GetBottomVeto() ) {
      vetoOuterHeight = shieldHeight + muonVetoLayers*muonVetoThickness;
    }

    //Teflon Shell parameters:
    teflonInnerR = (h2oLength-pmtLegLength - pmtMajorAxis - 7.5*cm)/2.0;
    teflonOuterR = (h2oLength-pmtLegLength - pmtMajorAxis - 7.5*cm   + 2.0*teflonThickness)/2.0;
    teflonHeight = (h2oHeight - 1.0*in)/2.0;
    teflonLength = (TMath::Sqrt(pow( (h2oLength-pmtMinorAxis-pmtLegLength - 7.5*cm),2) - pow((h2oLength/2.0),2) ))/3.415; //for front/back panels
  
  
    iUseBottomPMTs = input->GetBottomPMTs(); //PMTs on the bottom or just teflon reflector

    //DetermineSpacing must be called after all the parameters above are defined
    DetermineSpacing();
    std::cout << "PMT spacing determined\n";
    
    G4cerr << "done." << G4endl;

}//END of constructor

G4d2oCylindricalDetector::~G4d2oCylindricalDetector()
{
    G4cout << "Instance of G4d2oCylindricalDetector Destructed!" << G4endl;
    
}//END of destructor

void G4d2oCylindricalDetector::DetermineSpacing(){
    
    //* From Cylindrical Geometry:
      
      //minimum gap between PMTs and between PMTs and edge of volume
    G4double nominalSpacing = pmtMajorAxis + minGapBetweenPMTs;

    //spacing along FULL circumference
    //G4double totalLengthInCircle = 2*TMath::Pi()*h2oInnerRadius - 2*pmtDepth - 0.5*nominalSpacing - 2*minGapBetweenPMTs;
    //numInCircle = TMath::FloorNint(totalLengthInCircle/nominalSpacing);
    //spacingInCircle = totalLengthInCircle/double(numInCircle);
    
    //spacing along SIDES OF circumference: //*New!!
    totalAngle = TMath::ATan(TCRadius/h2oInnerRadius); //*New!!
    
    G4double totalLengthInCircle = 2*totalAngle*h2oInnerRadius; //*New!!
    numInCircle = TMath::FloorNint(totalLengthInCircle/nominalSpacing);
    spacingInCircle = totalLengthInCircle/double(numInCircle);
    
    G4double nominalRowSpacing = TMath::Sqrt( pow(pmtMajorAxis,2.0) - pow(0.5*spacingInCircle,2.0) ) + minGapBetweenPMTs;
    
    
    G4double totalWidthInZ = h2oHeight  - 2*minGapBetweenPMTs - 2*(pmtMajorAxis-nominalRowSpacing);
    numInZ = TMath::FloorNint(totalWidthInZ/nominalRowSpacing);
    spacingInZ = totalWidthInZ/double(numInZ);
    
    //some totals that get passed to other classes
    numRows = numInZ;
    numInRows = numInCircle;
    totPMT = numRows*numInRows;
    
    sideLiningX = h2oInnerRadius; //not using for now
    sideLiningZ = h2oInnerRadius; //not using for now
      
  
}

G4LogicalVolume * G4d2oCylindricalDetector::GetDetector(){
    
    //Materials pointer
    matPtr = G4d2oRunAction::GetMaterialsPointer();
    
    ///// create total detector volume /////
    G4LogicalVolume *totalDetLogV = GetTotalDetectorLogV();

    ///// create lead shielding volume /////
    G4LogicalVolume *shieldingLogV = GetShieldingLogV();
    
    ///// create steel vessel /////
    G4LogicalVolume *outerVesselLogV = GetOuterVesselLogV();
    
    ///// Create Teflon Lining /////
    //G4LogicalVolume *teflonLining = GetTeflonLiningLogV();
    
     //Left/Right Panel
        G4LogicalVolume *teflonLiningL = GetTeflonLiningLLogV();

     // Front/Back Panel
        G4LogicalVolume *teflonLiningF = GetTeflonLiningFLogV();

     //Top/Bottom Panel
        G4LogicalVolume *teflonLiningT = GetTeflonLiningTLogV();
    

    ///// H2O Inner /////
    G4LogicalVolume *h2oLogV = GetH2OLogV();
    
    ///// PMMA /////
    G4LogicalVolume *acrylicLogV = GetAcrylicLogV();
    
    ///// D2O /////
    G4LogicalVolume *d2oLogV = GetD2OLogV();
    
    ///// PMTs with mu-metal shield /////
    G4LogicalVolume *pmtLogV = 0;
    if (pmtDiameter < 1.0) {
      pmtLogV = GetEllipsoidPMT();
    } else {
      pmtLogV = GetSphericalPMT();
    }

    ///// InnerVeto Layer
    G4LogicalVolume *innerVetoLogV = GetInnerVetoLogV();

    /////OuterVeto Layer
    G4LogicalVolume *vetoLogV = GetVetoLogV();
    
    //     Geometry hierarchy
    //     totalDetLogV (mother volume)
    //       -> veto layers
    //        -> outer steel
    //           -> teflon lining
    //              -> H2O tank
    //                 -(1)-> PMTs
    //                 -(2)-> acrylic tank
    //                        -> D2O
    
    ///// Place the D2O in acrylic tank /////
    new G4PVPlacement(0,G4ThreeVector(0,0,0),d2oLogV,"d2oPhysV",acrylicLogV,false,0,true);
    
    ///// Place the acrylic tank in H2O /////
    new G4PVPlacement(0,G4ThreeVector(0,0,0),acrylicLogV,"acrylicPhysV",h2oLogV,false,0,true);
    
    ///// Place PMTs into H2O /////
    PlacePMTs(pmtLogV,h2oLogV);

    //// Create teflon end caps ////
    G4LogicalVolume *teflonCapTop = GetTeflonLiningCapLogV("topCap", true);
    G4LogicalVolume *teflonCapBot = GetTeflonLiningCapLogV("botCap", iUseBottomPMTs);
    
    ///// Place the H2O in teflon lining /////
    //G4VPhysicalVolume *h2oPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),h2oLogV,"h2oPhysV",teflonLining,false,0,true);
    G4VPhysicalVolume *h2oPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),h2oLogV,"h2oPhysV",outerVesselLogV,false,0,true);
    
    ///// Place the teflon lining in outer vessel/////

    // 4a: Left Curved panel
    
        G4VPhysicalVolume *teflonPhysLV = new G4PVPlacement(0,G4ThreeVector(0,0,0),teflonLiningL,"teflonPhysLV",h2oLogV,false,0,true);
        
        matPtr->SetReflector(h2oPhysV, teflonPhysLV, teflonReflectivity, teflonSigmaAlpha);
    
        // 4b: Right Curved panel
    
        G4RotationMatrix *RotR = new G4RotationMatrix();
        RotR -> rotateX(180.0*deg);
        RotR -> rotateY(180.0*deg);
        
        G4VPhysicalVolume *teflonPhysRV = new G4PVPlacement(RotR,G4ThreeVector(0,0,0),teflonLiningL,"teflonPhysRV",h2oLogV,false,0,true);
        
        matPtr->SetReflector(h2oPhysV, teflonPhysRV, teflonReflectivity, teflonSigmaAlpha);
    
        // 4c: Front panel
    
        G4VPhysicalVolume *teflonPhysFV = new G4PVPlacement(0,G4ThreeVector(-(h2oWidth-3.0*teflonThickness)/2.0,0.,0.),teflonLiningF,"teflonPhysFV",h2oLogV,false,0,true);
        
        matPtr->SetReflector(h2oPhysV, teflonPhysFV, teflonReflectivity, teflonSigmaAlpha);
    
        // 4d: Back  panel
        
        G4VPhysicalVolume *teflonPhysBV = new G4PVPlacement(0,G4ThreeVector((h2oWidth-3.0*teflonThickness)/2.0,0.,0.),teflonLiningF,"teflonPhysBV",h2oLogV,false,0,true);
        
        matPtr->SetReflector(h2oPhysV, teflonPhysBV, teflonReflectivity, teflonSigmaAlpha);
        
        // 4e: Top  panel
    
        G4VPhysicalVolume *teflonPhysTV = new G4PVPlacement(0,G4ThreeVector(0,0,(h2oHeight-2.34*cm)/2.0),teflonLiningT,"teflonPhysTV",h2oLogV,false,0,true);

        matPtr->SetReflector(h2oPhysV, teflonPhysTV, teflonReflectivity, teflonSigmaAlpha);
        
        // 4f: Bottom  panel
    
        G4VPhysicalVolume *teflonPhysBotV = new G4PVPlacement(0,G4ThreeVector(0,0,-(h2oHeight-2.34*cm)/2.0),teflonLiningT,"teflonPhysBotV",h2oLogV,false,0,true);
        
        matPtr->SetReflector(h2oPhysV, teflonPhysBotV, teflonReflectivity, teflonSigmaAlpha);
        
        
        
    /*G4VPhysicalVolume *teflonPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,-(pmtLegLength+pmtMinorAxis-h2oTCthickness+2.0*teflonThickness)/2.0),teflonLining,"teflonPhysV",h2oLogV,false,0,true); 
    
    matPtr->SetReflector(h2oPhysV, teflonPhysV, teflonReflectivity, teflonSigmaAlpha);*/

    /*// Place the teflon caps
    G4VPhysicalVolume *teflonCapTopPhysV = 0;
    G4VPhysicalVolume *teflonCapBotPhysV = 0;
    if (teflonCapBot) {
      if (pmtDiameter < 0.0)
        teflonCapBotPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -(h2oHeight-teflonThickness)/2.0 ), teflonCapBot, "teflonCapBotPhysV", h2oLogV, false, 0, true);//!!!!
      else 
        teflonCapBotPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -(h2oHeight-teflonThickness)/2.0), teflonCapBot, "teflonCapBotPhysV", h2oLogV, false, 0, true);
      matPtr->SetReflector(h2oPhysV, teflonCapBotPhysV, teflonReflectivity, teflonSigmaAlpha);
    }
    if (teflonCapTop) {
      if (pmtDiameter < 0.0)
        teflonCapTopPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, (h2oHeight-pmtMinorAxis)/2.0 - pmtLegLength), teflonCapTop, "teflonCapTopPhysV", h2oLogV, false, 0, true);
      else 
        teflonCapTopPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, (h2oHeight-teflonThickness)/2.0), teflonCapTop, "teflonCapTopPhysV", h2oLogV, false, 0, true);
      matPtr->SetReflector(h2oPhysV, teflonCapTopPhysV, teflonReflectivity, teflonSigmaAlpha);
    } */
    
    ///// Place the outer vessel inside the shielding
   /* G4double offset = 0.0;
    if (!input->GetBottomShielding() )
      offset = -0.5*shieldThickness;
    G4VPhysicalVolume *outerVesselPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,offset),outerVesselLogV,"outerVesselPhysV",shieldingLogV,false,0,true);

    ///// Place the Pb shielding inside the veto
    offset = 0.0;
    if (!input->GetBottomVeto() ) 
      offset = -0.5*muonVetoThickness;
    G4VPhysicalVolume *shieldingPhysV = new G4PVPlacement(0, G4ThreeVector(0,0,offset), shieldingLogV, "shieldingPhysV", innerVetoLogV, false, 0, true);
    
    /// // Place the inner veto inside the outer veto layer
    G4VPhysicalVolume *outerVetoPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,offset),innerVetoLogV,"innerVetoPhysV",vetoLogV,false,0,true);*/

    ///// Place the outer vessel in total detector/////

    new G4PVPlacement(0,G4ThreeVector(0,0,0),outerVesselLogV,"outerVesselLogV",totalDetLogV,false,0,true);
    
    // Add Event Handlers for Muon Vetos
    // Muon Vetos Inner
   /* G4d2oSensitiveDetector *senDetmuVetoInner = new G4d2oSensitiveDetector("muVetoInner", innerVetoLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoInner);
    innerVetoLogV->SetSensitiveDetector(senDetmuVetoInner);
    G4d2oSensitiveDetector *senDetmuVetoOuter = new G4d2oSensitiveDetector("muVetoOuter", vetoLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoOuter);
    vetoLogV->SetSensitiveDetector(senDetmuVetoOuter); */

    return totalDetLogV;
    
    
}//END of GetDetector()

G4LogicalVolume *G4d2oCylindricalDetector::GetTotalDetectorLogV(){
    
    //define the thickness of the side lining
    G4Box* outerSolid = new G4Box("outerSolid",contOuterLength/2.0,contOuterWidth/2.0,contOuterHeight/2.0);
    G4LogicalVolume *totalDetLogV = new G4LogicalVolume(outerSolid,
                                                        matPtr->GetMaterial( AIR ),
                                                        "totalDetLogV");

  G4VisAttributes *visTotal = new G4VisAttributes();
  visTotal->SetColour(G4Color::Green());
  visTotal->SetForceSolid(true);
  totalDetLogV->SetVisAttributes( visTotal );

    return totalDetLogV;
    
}

G4LogicalVolume *G4d2oCylindricalDetector::GetShieldingLogV() {
  G4LogicalVolume * shieldingLogV = 0;

  // Lead sheidling around the steel vessel
  // Do we need to do anything special in the event that the shielding thickness is set to 0? To prevent any sort of geant issues with teh geo?
  G4Tubs * shieldVessel = new G4Tubs("shieldVessel", 0.0, shieldLength/2.0, shieldHeight/2.0, 0.0, 360.0*deg);
  shieldingLogV = new G4LogicalVolume(shieldVessel,
                                    matPtr->GetMaterial( LEAD ),
                                    "shieldingLogV");

  G4VisAttributes *visShielding = new G4VisAttributes();
  visShielding->SetColour(G4Color::Gray());
  visShielding->SetForceSolid(true);
  shieldingLogV->SetVisAttributes( visShielding );

  return shieldingLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetOuterVesselLogV(){
    
    //Full detector volume will be steel for now. Steel will be displaced
    //as we add materials inside it. Steel has no optical parameters, so any photons that pass from
    //the H2O to the steel will be killed
    G4Box* outerVessel = new G4Box("outerVessel",contOuterLength/2.0,contOuterWidth/2.0,contOuterHeight/2.0);
    G4LogicalVolume *outerVesselLogV = new G4LogicalVolume(outerVessel,
                                                           matPtr->GetMaterial( STEEL ),
                                                           "outerVesselLogV");
    //    outerVesselLogV->SetVisAttributes(G4VisAttributes::Invisible);

    return outerVesselLogV;
    
}


G4LogicalVolume *G4d2oCylindricalDetector::GetInnerVetoLogV(){
    G4double halfheight = shieldHeight/2.0+1.0*muonVetoThickness;
    if (!input->GetBottomVeto() )
      halfheight = (shieldHeight + muonVetoThickness) / 2.0;
    G4Tubs* innerVeto = new G4Tubs("v",0.0,shieldLength/2.0+1.0*muonVetoThickness,
                                   halfheight,0.0,360.0*deg);
    G4LogicalVolume *innerVetoLogV = new G4LogicalVolume(innerVeto,
                                                           matPtr->GetMaterial( PLASTIC ),
                                                           "innerVetoLogV");
    return innerVetoLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetVetoLogV(){
    
    //Full veto volume will be acrylic for now.
    G4double halfheight = shieldHeight/2.0+2.0*muonVetoThickness;
    if (!input->GetBottomVeto() )
      halfheight = (shieldHeight + 2.0*muonVetoThickness)/2.0;
    G4Tubs* outerVeto = new G4Tubs("v",0.0,shieldLength/2.0+2.0*muonVetoThickness,
                                   halfheight,0.0,360.0*deg);
    G4LogicalVolume *outerVetoLogV = new G4LogicalVolume(outerVeto,
                                                           matPtr->GetMaterial( PLASTIC ),
                                                           "outerVetoLogV");
    return outerVetoLogV;
    
}

G4LogicalVolume *G4d2oCylindricalDetector::GetH2OLogV(){
    
    G4Box* h2oSolid = new G4Box("h2oSolid",(h2oWidth-teflonThickness)/2.0,(h2oLength-teflonThickness)/2.0,(h2oHeight-teflonThickness)/2.0);
    
        
    G4LogicalVolume *h2oLogV = new G4LogicalVolume(h2oSolid,
                                                   matPtr->GetMaterial( H2O ),
                                                   "h2oLogV");
    
    G4VisAttributes *visLightWater = new G4VisAttributes();
    visLightWater->SetColour(G4Color::Cyan());
    visLightWater->SetForceSolid(true);
    h2oLogV->SetVisAttributes( visLightWater );

    return h2oLogV;
    
}

G4LogicalVolume *G4d2oCylindricalDetector::GetAcrylicLogV(){
    
    G4double acrylicLength = d2oLength + 2.0*acrylicThickness; //x
    G4double acrylicHeight = d2oHeight + 2.0*acrylicEndCapThickness; //z
    
    //make it solid acrylic - D2O gets placed inside and displaces central region
    G4Tubs* acrylicSolid = new G4Tubs("acrylicSolid",0.0,acrylicLength/2.0,acrylicHeight/2.0,0.0,360.0*deg);
    G4LogicalVolume *acrylicLogV = new G4LogicalVolume(acrylicSolid,
                                                       matPtr->GetMaterial( PMMA ),
                                                       "acrylicLogV");

    G4VisAttributes *visAcrylic = new G4VisAttributes();
    visAcrylic->SetColour(G4Color::Magenta());
    visAcrylic->SetForceSolid(true);
    acrylicLogV->SetVisAttributes( visAcrylic );

    return acrylicLogV;
    
}

G4LogicalVolume *G4d2oCylindricalDetector::GetD2OLogV(){
    
    G4Tubs* d2oSolid = new G4Tubs("d2oSolid",0.0,d2oLength/2.0,d2oHeight/2.0,0.0,360.0*deg);
    G4LogicalVolume *d2oLogV = new G4LogicalVolume(d2oSolid,
                                                   matPtr->GetMaterial( D2O ),
                                                   "d2oLogV");
    G4VisAttributes *visAttD2O = new G4VisAttributes();
    visAttD2O->SetColour(G4Color::Blue());
    visAttD2O->SetForceSolid(true);
    d2oLogV->SetVisAttributes( visAttD2O );

    return d2oLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetTeflonLiningCapLogV(const char * name, bool boundPMTs=true){
  G4LogicalVolume *capLogV = 0;

  G4double width = (h2oLength-teflonThickness)/2.0 - 1.0*mm;
  G4double height = (teflonThickness)/2.0;

  G4VSolid * capSolid = new G4Tubs(Form("%sSolid", name), 0.0, width, height, 0.0, 360.0*deg);

  if (boundPMTs) {
    G4Tubs *capSub = new G4Tubs("subCap", 0.0, pmtMajorAxis/2.0 + 1.0*mm, height+1.0*cm, 0.0 ,360.0);
    
    G4Transform3D transform;
    int counter = 0;
    for (int bin=1; bin<=hPanel[2]->GetNbinsX(); ++bin) {
      // For now we can assume that the bottom pmts are located directly below the top PMTS (ie no x/y offset)
      if (hPanel[2]->GetBinContent(bin) <= 0.0) continue;

      // For now we can assume that the bottom pmts are located directly below the top PMTS (ie no x/y offset)
      transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(hPanel[0]->GetBinContent(bin), hPanel[1]->GetBinContent(bin), -1.0*mm));

      capSolid = new G4SubtractionSolid(Form("%sSolid", name), capSolid, capSub, transform);
      ++counter;
    }
    capLogV = new G4LogicalVolume(capSolid,
                                  matPtr->GetMaterial( TEFLON ),
                                  Form("%sLogV", name) );

    G4VisAttributes *visCap = new G4VisAttributes();
    visCap->SetColour(G4Color::Red());
    visCap->SetForceSolid(true);
    capLogV->SetVisAttributes( visCap );
  } else {
    capLogV = new G4LogicalVolume(capSolid,
                                  matPtr->GetMaterial( TEFLON ),
                                  Form("%sLogV", name) );

    G4VisAttributes *visCap = new G4VisAttributes();
    visCap->SetColour(G4Color::Red());
    visCap->SetForceSolid(true);
    capLogV->SetVisAttributes( visCap );
  }


  return capLogV;
}

G4LogicalVolume *G4d2oCylindricalDetector::GetTeflonLiningLogV(){

    
    G4double teflonheight = (d2oHeight/2.0) + acrylicEndCapThickness + h2oTCthickness + (pmtMinorAxis/2.0) - (pmtLegLength+3.0*teflonThickness)/2.0; 
    
    G4Tubs * teflonSolid = new G4Tubs("teflonSolid",(h2oLength-teflonThickness)/2.0, h2oLength/2.0, teflonheight, 0.0, 360.0*deg);
    G4LogicalVolume *theSideLiningLogV = new G4LogicalVolume(teflonSolid,
                                                             matPtr->GetMaterial( TEFLON ),
                                                             "theSideLiningLogV");
    G4VisAttributes *visTeflonLining = new G4VisAttributes();
    visTeflonLining->SetColour(G4Color::Red());
    visTeflonLining->SetForceSolid(true);
    theSideLiningLogV->SetVisAttributes( visTeflonLining );
    
    return theSideLiningLogV;
    
}


//-----------Teflon Panels: needed to separate to avoid overlap issues---------------
//------------------------------------------------------------------------------------

//-----------Left/Right panel (with holes)

G4LogicalVolume *G4d2oCylindricalDetector::GetTeflonLiningLLogV(){
    
  G4double OPMTTopESemiX      = 103.6*mm;
  G4double OPMTTopESemiY      = OPMTTopESemiX;
  G4double OPMTTopESemiZ      = 76.72*mm;
  G4double OPMTTopEBottomCut  = 0;
  G4double OPMTTopETopCut     = OPMTTopESemiZ;
  G4double PMTFrameThickness  = pmtWindowThickness;
    
    // Left curved pannel, using G4SubtractionSolid (
    //G4Tubs* TefInnerL = new G4Tubs("TefInnerL",0.0,teflonInnerR,teflonHeight,28.51*deg,123.0*deg);
    
    G4Tubs* TefOuterL = new G4Tubs("TefOuterL",teflonInnerR,teflonOuterR,teflonHeight,28.51*deg,123.0*deg); 
    
    G4SubtractionSolid* teflonVesselL;// = new G4SubtractionSolid("teflonVesselL", TefOuterL, TefInnerL); 
  
    
    //----- Placing holes for PMTs +x side LEFT CURVED PANNEL ( 
    
    G4int iCopy = 0;
    G4double angularSpacing = 2.*totalAngle/numInCircle; //*New!!
    
    
    
    for(G4int iCopyC = 0; iCopyC<numInCircle; iCopyC++){

        for(G4int iCopyZ = 0; iCopyZ<numInZ; iCopyZ++){
	  
            G4double theZZ = -spacingInZ*numInZ/2.0 + double(iCopyZ+0.5)*spacingInZ;
 
            G4double angleInCircle = iCopyC*angularSpacing - angularSpacing;
            G4RotationMatrix *rotH = new G4RotationMatrix();
            
            
            rotH -> rotateZ(angleInCircle + 180.0*deg); // PMT surrounding the cylinder NEW
            rotH->rotateY(90.0*deg); // PMT laying down
            rotH->rotateX(90.0*deg); //PMT facing the fiducial volume
            
            
            G4double theX = h2oInnerRadius*TMath::Sin(angleInCircle);
            G4double theYY = h2oInnerRadius*TMath::Cos(angleInCircle);
            
            G4ThreeVector thisPMTPlacement(theX, theYY, theZZ);
                
            if(iCopyC == 0 && iCopyZ == 0) {
                
                //teflonVesselL = new G4SubtractionSolid(Form("teflonVesselL_%d",iCopy), TefOuterL, new G4Tubs(Form("teflonHoleL_%d",iCopy),0.0,(pmtMajorAxis + 1.0*mm)/2.0,pmtDepth/2.0,0.0,360*deg),rotH,thisPMTPlacement);
                teflonVesselL = new G4SubtractionSolid(Form("teflonVesselL_%d",iCopy), TefOuterL, new G4Ellipsoid(Form("teflonHoleL_%d",iCopy),OPMTTopESemiX, OPMTTopESemiY, OPMTTopESemiZ, OPMTTopEBottomCut, OPMTTopETopCut),rotH,thisPMTPlacement);
                
                
            } 
            else{
                
                teflonVesselL = new G4SubtractionSolid(Form("teflonVesselL_%d",iCopy), teflonVesselL, new G4Ellipsoid(Form("teflonHoleL_%d",iCopy),OPMTTopESemiX, OPMTTopESemiY, OPMTTopESemiZ, OPMTTopEBottomCut, OPMTTopETopCut),rotH,thisPMTPlacement);
                
            }
            
            

            
            iCopy++;
                            
        }
        
    } 
    
    
    
    G4LogicalVolume *theSideLiningLogLV = new G4LogicalVolume(teflonVesselL,
                                                             matPtr->GetMaterial( TEFLON ),
                                                             "theSideLiningLogLV");
    
    G4VisAttributes *visAttTeflon = new G4VisAttributes();
    visAttTeflon->SetColour(G4Color::Red());
    visAttTeflon->SetForceSolid(true);
    theSideLiningLogLV->SetVisAttributes( visAttTeflon );
    
    return theSideLiningLogLV;
    
}


//-------------------Front/Back panel

G4LogicalVolume *G4d2oCylindricalDetector::GetTeflonLiningFLogV(){
    

    G4Box* teflonSolidF = new G4Box("teflonSolidF",teflonThickness/2.0, teflonLength, teflonHeight );
    
    
    G4LogicalVolume *theSideLiningFLogV = new G4LogicalVolume(teflonSolidF,
                                                             matPtr->GetMaterial( TEFLON ),
                                                             "theSideLiningFLogV");
    
    G4VisAttributes *visAttTeflon = new G4VisAttributes();
    visAttTeflon->SetColour(G4Color::Red());
    visAttTeflon->SetForceSolid(true);
    theSideLiningFLogV->SetVisAttributes( visAttTeflon );
    
    return theSideLiningFLogV;
    
    
                }
                
                

//-------------------Top/Bottom panel

G4LogicalVolume *G4d2oCylindricalDetector::GetTeflonLiningTLogV(){
    
    G4Tubs* CylTopSolid = new G4Tubs("CylTopSolid",0.0,teflonOuterR, teflonThickness/2.0,0.0*deg,360.0*deg); 

    G4Box* BoxTopSolid = new G4Box("BoxTopSolid", (h2oWidth - 2.0*teflonThickness)/2.0, teflonInnerR, teflonThickness/2.0 ); 

    G4IntersectionSolid* teflonSolidT =
    new G4IntersectionSolid("teflonSolidT", BoxTopSolid, CylTopSolid);
  
    G4LogicalVolume *theSideLiningTLogV = new G4LogicalVolume(teflonSolidT,
                                                             matPtr->GetMaterial( TEFLON ),
                                                             "theSideLiningTLogV");
    
    G4VisAttributes *visAttTeflon = new G4VisAttributes();
    visAttTeflon->SetColour(G4Color::Red());
    visAttTeflon->SetForceSolid(true);
    theSideLiningTLogV->SetVisAttributes( visAttTeflon );
    
    return theSideLiningTLogV;
    
}


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------




void G4d2oCylindricalDetector::PlacePMTs(G4LogicalVolume *thePMTLogV, G4LogicalVolume *theMotherLogV){

    
    G4int iCopy = 0;
    
    //****************************** x-y plane, PMTs on +x side *******************************
    
    //angularSpacing = 2.*TMath::Pi()/numInCircle;
    G4double angularSpacing = 2.*totalAngle/numInCircle; //*New!!

    for(G4int iCopyC = 0; iCopyC<numInCircle; iCopyC++){

        for(G4int iCopyZ = 0; iCopyZ<numInZ; iCopyZ++){
	  
            G4double theZZ = -spacingInZ*numInZ/2.0 + double(iCopyZ+0.5)*spacingInZ;
            //G4double theShift = 0.25*angularSpacing;
            G4double theShift = 0.0; //*New
            //if(iCopyZ%2==0) theShift = -0.25*angularSpacing;

            
            //G4double angleInCircle = iCopyC*angularSpacing + theShift;
            G4double angleInCircle = iCopyC*angularSpacing - angularSpacing;
            G4RotationMatrix *rotA = new G4RotationMatrix();
            //G4RotationMatrix *rotTop = new G4RotationMatrix();
             
            //rotA->rotateY(180.0*deg); //PMT head pointing downwards!
            //rotA -> rotateZ(TMath::Pi() - angleInCircle); // PMT surrounding the cylinder
            rotA -> rotateZ(angleInCircle + 180.0*deg); // PMT surrounding the cylinder NEW
            
            rotA->rotateY(90.0*deg); // PMT laying down
            rotA->rotateX(90.0*deg); //PMT facing the fiducial volume
            

            G4double theX = h2oInnerRadius*TMath::Sin(angleInCircle);
            G4double theYY = h2oInnerRadius*TMath::Cos(angleInCircle);
            

            G4ThreeVector thisPMTPlacement(theX, theYY, theZZ);
                
            new G4PVPlacement(rotA,thisPMTPlacement,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);
            
            
            //--- Placing housings:
            
            /*G4double housingShift = pmtHousingDepth/2.;
            G4ThreeVector thisPMTHousingPlacement(theX + housingShift*TMath::Sin(angleInCircle), theYY + housingShift*TMath::Cos(angleInCircle) , theZZ);

            G4RotationMatrix *rotAA = new G4RotationMatrix();
            
            rotAA -> rotateZ(angleInCircle + 180.0*deg); 
            rotAA -> rotateY(90.0*deg); 
            rotAA -> rotateX(90.0*deg); 
            
            new G4PVPlacement(rotAA,thisPMTHousingPlacement,pmtHousingLogV,Form("pmtHousingPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);*/

            iCopy++;
                            

        }
        
    }
    
    // Now store PMT locations:
    
    
    iCopy = 0;
    
    //histograms that hold coordinates of the front faces of the PMTs
    for(G4int iDim=0; iDim<3; iDim++)
        hPanel[iDim] = new TH1D(Form("hPanel[%d]",iDim),Form("Dimension %d",iDim),totPMT,-0.5,totPMT-0.5);
    
    //angularSpacing = 2.*TMath::Pi()/numInCircle;
     angularSpacing = 2.*totalAngle/numInCircle; //*New!!

    for(G4int iCopyC = 0; iCopyC<numInCircle; iCopyC++){

        for(G4int iCopyZ = 0; iCopyZ<numInZ; iCopyZ++){
	  
            G4double theZZ = -spacingInZ*numInZ/2.0 + double(iCopyZ+0.5)*spacingInZ;
            //G4double theShift = 0.25*angularSpacing;
            G4double theShift = 0.0; //*New
            //if(iCopyZ%2==0) theShift = -0.25*angularSpacing;

	    
            //G4double angleInCircle = iCopyC*angularSpacing + theShift;
            G4double angleInCircle = iCopyC*angularSpacing - angularSpacing;
            G4RotationMatrix *rotA = new G4RotationMatrix();
            
            //rotA->rotateY(180.0*deg); //PMT head pointing downwards!
            //rotA -> rotateZ(TMath::Pi() - angleInCircle); // PMT surrounding the cylinder
            rotA -> rotateZ(angleInCircle + 180.0*deg); // PMT surrounding the cylinder NEW
            
            rotA->rotateY(90.0*deg); // PMT laying down
            rotA->rotateX(90.0*deg); //PMT facing the fiducial volume
            

            G4double theX = h2oInnerRadius*TMath::Sin(angleInCircle);
            G4double theYY = h2oInnerRadius*TMath::Cos(angleInCircle);
            

            G4ThreeVector thisPMTPlacement(theX, theYY, theZZ);
            
            
            iCopy++;
            
            //save coordinates of the front faces of the PMTs
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x()-pmtMinorAxis/2.0 - pmtLegLength);
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y()-pmtMinorAxis/2.0 - pmtLegLength);
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z());
            
            
            
            }
        
    } 
    
    
//*New!!:

    
    iCopy = 0;
    
    //****************************** x-y plane, PMTs on -x side *******************************

        angularSpacing = 2.*totalAngle/numInCircle; 

    for(G4int iCopyC = 0; iCopyC<numInCircle; iCopyC++){

        for(G4int iCopyZ = 0; iCopyZ<numInZ; iCopyZ++){
	  
            G4double theZZ = -spacingInZ*numInZ/2.0 + double(iCopyZ+0.5)*spacingInZ;
            
            G4double angleInCircle = iCopyC*angularSpacing - angularSpacing;
            G4RotationMatrix *rotB = new G4RotationMatrix();
             
            rotB -> rotateZ(angleInCircle + 180.0*deg); // PMT surrounding the cylinder NEW
            
            rotB->rotateY(90.0*deg); // PMT laying down
            rotB->rotateX(-90.0*deg); //PMT facing the fiducial volume
            

            G4double theX = h2oInnerRadius*TMath::Sin(angleInCircle);
            G4double theYY = h2oInnerRadius*TMath::Cos(angleInCircle);
            

            G4ThreeVector thisPMTPlacement2(-theX, -theYY, theZZ);
                
            new G4PVPlacement(rotB,thisPMTPlacement2,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);
            
            //--- Placing housings:
            
            /*G4double housingShift = pmtHousingDepth/2.;
            G4ThreeVector thisPMTHousingPlacement2(-theX - housingShift*TMath::Sin(angleInCircle), -theYY - housingShift*TMath::Cos(angleInCircle) , theZZ);

            G4RotationMatrix *rotBB = new G4RotationMatrix();
            
            rotBB -> rotateZ(angleInCircle + 180.0*deg); 
            rotBB -> rotateY(90.0*deg); 
            rotBB -> rotateX(-90.0*deg); 
            
            new G4PVPlacement(rotBB,thisPMTHousingPlacement2,pmtHousingLogV,Form("pmtHousingPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);*/
            
            

            iCopy++;
                            

        }
        
    }
    
    // Now store PMT locations:
    
    
    iCopy = 0;
    
    //histograms that hold coordinates of the front faces of the PMTs
    for(G4int iDim=0; iDim<3; iDim++)
        hPanel[iDim] = new TH1D(Form("hPanel[%d]",iDim),Form("Dimension %d",iDim),totPMT,-0.5,totPMT-0.5);
    
     angularSpacing = 2.*totalAngle/numInCircle; //*New!!

    for(G4int iCopyC = 0; iCopyC<numInCircle; iCopyC++){

        for(G4int iCopyZ = 0; iCopyZ<numInZ; iCopyZ++){
	  
            G4double theZZ = -spacingInZ*numInZ/2.0 + double(iCopyZ+0.5)*spacingInZ;
	    
            G4double angleInCircle = iCopyC*angularSpacing - angularSpacing;
            G4RotationMatrix *rotB = new G4RotationMatrix();

            rotB -> rotateZ(angleInCircle + 180.0*deg); // PMT surrounding the cylinder NEW
            
            rotB->rotateY(90.0*deg); // PMT laying down
            rotB->rotateX(-90.0*deg); //PMT facing the fiducial volume
            

            G4double theX = h2oInnerRadius*TMath::Sin(angleInCircle);
            G4double theYY = h2oInnerRadius*TMath::Cos(angleInCircle);
            

            G4ThreeVector thisPMTPlacement2(-theX, -theYY, theZZ);
            
            
            iCopy++;
            
            //save coordinates of the front faces of the PMTs
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement2.x()-pmtMinorAxis/2.0 - pmtLegLength);
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement2.y()-pmtMinorAxis/2.0 - pmtLegLength);
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement2.z());
            
            
            
            }
        
    }
    
    
}




