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

//NEED TO COORDINATE DIMENSIONS WITH D2O GROUP 


#include "G4d2oDetector.hh"
#include "G4d2oRunAction.hh"
#include "G4d2oSensitiveDetector.hh"

#include "TMath.h"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Ellipsoid.hh" //*New
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4SDManager.hh"
#include "G4PVParameterised.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh" //*New
#include "G4IntersectionSolid.hh" //**New
#include "G4MultiUnion.hh" //**New
#include "G4OpticalSurface.hh" //NEW
#include "G4LogicalBorderSurface.hh" //NEW

using namespace std;

G4d2oDetector::G4d2oDetector(): in(2.54*cm){
    G4cerr << "\n\tBuilding vertical cylindrical geometry" << G4endl ;
    
    Initialize();
}

void G4d2oDetector::Initialize(){
    
  //relevant parameters
  
  //D2O
  d2oLength = 69.99*cm;
  d2oWidth = d2oLength; //diameter
  d2oHeight = 161.544*cm;//z //Using bottom PMT space for d2o 
  
  //acrylic
  acrylicThickness = 0.25*in;
  acrylicEndCapThickness = 1.0*in;
  
  //pmt
  pmtDepth = 10.0*in;
  //pmtDepth = 21.55*cm; //**New
  muMetalThickness = 0.03*in;
  pmtWindowThickness = 3.0*mm;  //making this up
  minGapBetweenPMTs = 3.0*mm;
  //input = inputVariables::GetIVPointer();
  pmtDiameter = input->GetPMTDiameter()*in; //pmtDiameter controlled from beamOn.dat or command line. 
  pmtMinorAxis = pmtMajorAxis = pmtDiameter;
  
     pmtLegLength = 62.*mm;
     if (pmtDiameter < 0.0) {
       pmtMinorAxis = 76.72*mm * 2.0;
       pmtMajorAxis = 103.6*mm * 2.0;
     }

  //teflon
  teflonThickness = 0.25*in;
  teflonReflectivity = 0.97; //need to verify
  teflonSigmaAlpha = 0.1; //roughness parameter. use 0 for polished
  
  
  //H2O
  //h2oTCthickness = input->GetTailCatcherThick()*cm;
  h2oTCthickness = 10.0*cm;
  h2oLength = d2oLength + 2*(acrylicThickness + h2oTCthickness); //x
  h2oWidth = h2oLength; //y
    if (pmtDiameter > 0.0) {
      h2oHeight = d2oHeight + 2*(acrylicEndCapThickness + h2oTCthickness + pmtMinorAxis/2.0); //z
    } else {
      h2oHeight = d2oHeight + 2*(acrylicEndCapThickness + h2oTCthickness + pmtMinorAxis + pmtLegLength); //z
    }
  
  
  //not using these for now, but set them anyway so that other parts of code don't have problems
  h2oInnerLength = h2oLength;
  h2oInnerWidth = h2oWidth;
  h2oInnerHeight = h2oHeight;
  
  //**New:
  h2oInnerRadius = d2oWidth/2.0 + acrylicThickness + h2oTCthickness + pmtMinorAxis/2.0;
  TCRadius = d2oWidth/2.0 + acrylicThickness + h2oTCthickness;
  
  //outer vessel
  outerContainerThickness = 0.25*in;
  
  ///// Full detector volume /////
  contOuterLength = h2oLength + 2*outerContainerThickness; //x
  contOuterWidth = contOuterLength; //y
  contOuterHeight = h2oHeight + 2*outerContainerThickness; //z
  
    ///// Shielding volume /////
    //shieldThickness = input->GetShieldThickness()*in;
    shieldThickness = 2.0*in;
    shieldLength = contOuterLength + 2.0*shieldThickness; // x
    shieldWidth = shieldLength;
    shieldHeight = contOuterHeight + 2.0*shieldThickness; 
    
   //Veto Layers
    muonVetoThickness = 2.54*cm;
    muonVetoLayers = 2;
    
    iUseBottomVeto = false;
    vetoOuterLength = shieldLength + 2.0*muonVetoLayers*muonVetoThickness; //x
    vetoOuterWidth = shieldWidth + 2.0*muonVetoLayers*muonVetoThickness; //y
    vetoOuterHeight = shieldHeight + 2.0*muonVetoLayers*muonVetoThickness; //z
    if (iUseBottomVeto )
      vetoOuterHeight = shieldHeight + muonVetoLayers*muonVetoThickness; 
  
  iUseBottomPMTs = false;
  //iUseBottomPMTs = input->GetBottomPMTs(); //PMTs on the bottom or just teflon reflector

  std::cout << "All parameters set!\n";

  //DetermineSpacing must be called after all the parameters above are defined
  DetermineSpacing();
  std::cout << "PMT spacing determined\n";

  G4cerr << "done." << G4endl;
  
}//END of constructor



G4d2oDetector::~G4d2oDetector(){
    G4cout << "Instance of G4d2oDetector Destructed!" << G4endl;
    
}//END of destructor



void G4d2oDetector::DetermineSpacing(){

spacingInX = pmtMajorAxis + (2.0*minGapBetweenPMTs);
  numInX = TMath::FloorNint( (h2oLength / spacingInX) );

  spacingInY = TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + 2.0*minGapBetweenPMTs);
  numInY = TMath::FloorNint( (h2oLength / spacingInY) );
        
}



G4LogicalVolume * G4d2oDetector::GetDetector(){
  std::cout << "Getting Detector\n";

    //Materials pointer
    matPtr = G4d2oRunAction::GetMaterialsPointer();
    
    ///// create total detector volume /////
    G4LogicalVolume *totalDetLogV = GetTotalDetectorLogV();
    
    ///// create lead shielding volume /////
    G4LogicalVolume *shieldingLogV = GetShieldingLogV();
    
    ///// create steel vessel /////
    G4LogicalVolume *outerVesselLogV = GetOuterVesselLogV();
    
    ///// Create Teflon Lining /////
    G4LogicalVolume *teflonLining = GetTeflonLiningLogV();
    
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
      //pmtLogV = GetSphericalPMT();
      pmtLogV = GetEllipsoidPMT();
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
    
    ///// 1.- Place the D2O in acrylic tank ///// 
    new G4PVPlacement(0,G4ThreeVector(0,0,0),d2oLogV,"d2oPhysV",acrylicLogV,false,0,true);
    
    ///// 2.- Place the acrylic tank in H2O /////
    new G4PVPlacement(0,G4ThreeVector(0,0,0),acrylicLogV,"acrylicPhysV",h2oLogV,false,0,true);
    
    ///// 3.- Place PMTs into H2O /////
    PlacePMTs(pmtLogV,h2oLogV);
    
    //// Create teflon end caps ////
    G4LogicalVolume *teflonCapTop = GetTeflonLiningCapLogV("topCap", true);
    G4LogicalVolume *teflonCapBot = GetTeflonLiningCapLogV("botCap", iUseBottomPMTs);
    
    ///// 4.- Place the H2O in teflon lining /////
    G4VPhysicalVolume *h2oPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),h2oLogV,"h2oPhysV",teflonLining,false,0,true);
    //G4VPhysicalVolume *h2oPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),h2oLogV,"h2oPhysV",outerVesselLogV,false,0,true);
    
    ///// 5.- Place the teflon lining in outer vessel/////
    G4VPhysicalVolume *teflonPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),teflonLining,"teflonPhysV",outerVesselLogV,false,0,true);
    //G4VPhysicalVolume *teflonPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),teflonLining,"teflonPhysV",h2oLogV,false,0,true);
    //set reflectivity of teflon
    matPtr->SetReflector(h2oPhysV, teflonPhysV, teflonReflectivity, teflonSigmaAlpha);
    
    // Place the teflon caps
    G4VPhysicalVolume *teflonCapTopPhysV = 0;
    G4VPhysicalVolume *teflonCapBotPhysV = 0;
    if (teflonCapBot) {
      if (pmtDiameter < 0.0)
        teflonCapBotPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -((h2oHeight-pmtMinorAxis)/2.0 - pmtLegLength)), teflonCapBot, "teflonCapBotPhysV", h2oLogV, false, 0, true);
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
    } 
    
    ///// Place the outer vessel inside the shielding
    G4VPhysicalVolume *outerVesselPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),outerVesselLogV,"outerVesselPhysV",shieldingLogV,false,0,true);
    
    ///// Place the Pb shielding inside the veto
    G4double offset = 0.0;
    if (iUseBottomVeto) 
      offset = -0.5*muonVetoThickness;
    G4VPhysicalVolume *shieldingPhysV = new G4PVPlacement(0, G4ThreeVector(0,0,offset), shieldingLogV, "shieldingPhysV", innerVetoLogV, false, 0, true);
    
    /// // Place the inner veto inside the outer veto layer
    G4VPhysicalVolume *outerVetoPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,offset),innerVetoLogV,"innerVetoPhysV",vetoLogV,false,0,true);
    
    ///// Place the outer vessel in total detector/////
    new G4PVPlacement(0,G4ThreeVector(0,0,0),vetoLogV,"vetoPhysV",totalDetLogV,false,0,true);
        
    // Add Event Handlers for Muon Vetos
    // Muon Vetos Inner
    G4SDManager *sdManager = G4SDManager::GetSDMpointer(); //New
    
    G4d2oSensitiveDetector *senDetmuVetoInner = new G4d2oSensitiveDetector("muVetoInner", innerVetoLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoInner);
    innerVetoLogV->SetSensitiveDetector(senDetmuVetoInner);
    G4d2oSensitiveDetector *senDetmuVetoOuter = new G4d2oSensitiveDetector("muVetoOuter", vetoLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoOuter);
    vetoLogV->SetSensitiveDetector(senDetmuVetoOuter);
    
    ///// 7.- Place the outer vessel in total detector/////
    //new G4PVPlacement(0,G4ThreeVector(0,0,0),outerVesselLogV,"outerVesselPhysV",totalDetLogV,false,0,true);
    
    return totalDetLogV;
    
    
}//END of GetDetector()


G4LogicalVolume *G4d2oDetector::GetTotalDetectorLogV(){
    
    //define the thickness of the side lining
    G4Tubs* outerSolid = new G4Tubs("outerSolid",0.0,vetoOuterLength/2.0,vetoOuterHeight/2.0,0.0,360.0*deg);
    G4LogicalVolume *totalDetLogV = new G4LogicalVolume(outerSolid,
                                                        matPtr->GetMaterial( AIR ),
                                                        "totalDetLogV");

  G4VisAttributes *visTotal = new G4VisAttributes();
  visTotal->SetColour(G4Color::Green());
  visTotal->SetForceSolid(true);
  totalDetLogV->SetVisAttributes( visTotal );

    return totalDetLogV;
    
}



G4LogicalVolume *G4d2oDetector::GetShieldingLogV() {
  G4LogicalVolume * shieldingLogV = 0;

  // Lead sheidling around the steel vessel
  // Do we need to do anything special in the event that the shielding thickness is set to 0? To prevent any sort of geant issues with the geo?
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

G4LogicalVolume *G4d2oDetector::GetOuterVesselLogV(){
    
    //Full detector volume will be steel for now. Steel will be displaced
    //as we add materials inside it. Steel has no optical parameters, so any photons that pass from
    //the H2O to the steel will be killed
    G4Tubs* outerVessel = new G4Tubs("outerVessel",0.0,contOuterLength/2.0,contOuterHeight/2.0,0.0,360.0*deg);
    G4LogicalVolume *outerVesselLogV = new G4LogicalVolume(outerVessel,
                                                           matPtr->GetMaterial( STEEL ),
                                                           "outerVesselLogV");
    //    outerVesselLogV->SetVisAttributes(G4VisAttributes::Invisible);

    return outerVesselLogV;
    
}


G4LogicalVolume *G4d2oDetector::GetInnerVetoLogV(){
    
    G4double halfheight = shieldHeight/2.0 + 1.0*muonVetoThickness;
    if (iUseBottomVeto )
      halfheight = (shieldHeight + muonVetoThickness) / 2.0;
    G4Tubs* innerVeto = new G4Tubs("v",0.0,shieldLength/2.0 + 1.0*muonVetoThickness,
                                   halfheight,0.0,360.0*deg);
    G4LogicalVolume *innerVetoLogV = new G4LogicalVolume(innerVeto,
                                                           matPtr->GetMaterial( PLASTIC ),
                                                           "innerVetoLogV");
    return innerVetoLogV;
}


G4LogicalVolume *G4d2oDetector::GetVetoLogV(){
    
    //Full veto volume will be acrylic for now.
    G4double halfheight = shieldHeight/2.0 + 2.0*muonVetoThickness;
    if (iUseBottomVeto )
      halfheight = (shieldHeight + 2.0*muonVetoThickness)/2.0;
    G4Tubs* outerVeto = new G4Tubs("v",0.0,shieldLength/2.0 + 2.0*muonVetoThickness,
                                   halfheight,0.0,360.0*deg);
    G4LogicalVolume *outerVetoLogV = new G4LogicalVolume(outerVeto,
                                                           matPtr->GetMaterial( PLASTIC ),
                                                           "outerVetoLogV");
    return outerVetoLogV;
    
}



G4LogicalVolume *G4d2oDetector::GetH2OLogV(){
    
    //G4Tubs* h2oSolid = new G4Tubs("h2oSolid",0.0,(h2oLength-teflonThickness)/2.0,
    //                              (h2oHeight-teflonThickness)/2.0,0.0,360.0*deg);
    G4Tubs* h2oSolid = new G4Tubs("h2oSolid",0.0,(h2oLength)/2.0,
                                  (h2oHeight)/2.0,0.0,360.0*deg);
    G4LogicalVolume *h2oLogV = new G4LogicalVolume(h2oSolid,
                                                   matPtr->GetMaterial( H2O ),
                                                   "h2oLogV");
    
    G4VisAttributes *visLightWater = new G4VisAttributes();
    visLightWater->SetColour(G4Color::Cyan());
    visLightWater->SetForceSolid(true);
    h2oLogV->SetVisAttributes( visLightWater );

    return h2oLogV;
    
}



G4LogicalVolume *G4d2oDetector::GetAcrylicLogV(){
    
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

G4LogicalVolume *G4d2oDetector::GetD2OLogV(){
    
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



G4LogicalVolume *G4d2oDetector::GetTeflonLiningCapLogV(const char * name, bool boundPMTs=true){
    
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

      capSolid = new G4SubtractionSolid(Form("%sSolid", name), capSolid, capSub, transform); //HOLES?
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




G4LogicalVolume *G4d2oDetector::GetTeflonLiningLogV(){
    
  //d2oHeight + 2*(acrylicThickness + h2oTCthickness + pmtMinorAxis + pmtLegLength);
    //teflonThickness = 2.55*in;
    G4double teflonheight = (d2oHeight/2.0) + acrylicEndCapThickness + h2oTCthickness + (pmtMinorAxis/2.0);
    //d2oHeight + 2.0*(acrylicThickness + h2oTCthickness) + pmtMinorAxis;
    
    //G4Tubs* teflonSolid = new G4Tubs("teflonSolid",0.0,h2oLength/2.0,h2oHeight/2.0,0.0,360.0*deg);
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



void G4d2oDetector::PlacePMTs(G4LogicalVolume *thePMTLogV, G4LogicalVolume *theMotherLogV){



    
    G4int iCopy = 0;
    
    //histograms that hold coordinates of the front faces of the PMTs
    //x-y plane, +z
    G4RotationMatrix *rotTop = new G4RotationMatrix();
    rotTop->rotateY(180.0*deg);

    if (!(numInY%2) && !(numInX%2)) {
      G4double yshift = 0.0;
      if (!(numInY%2)) {
        yshift = (pmtMajorAxis / 2.0) / TMath::Sqrt(3.0);
      }
      for (int iY=0; iY<numInY; ++iY) {
        G4int ycount = iY - (numInY/2);
        G4double ycoord = static_cast<G4double>(ycount);
        ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
        ycoord += yshift;

        G4int nx = numInX - TMath::Abs(ycount);
        G4double xshift = 0.0;
        if (!(nx%2)) {
          xshift = pmtMajorAxis/2.0;
        }

        for (int iX=0; iX<nx; ++iX) {
          G4int xcount = iX - (nx/2);
          G4double xcoord = static_cast<G4double>(xcount);
          xcoord *= pmtMajorAxis + minGapBetweenPMTs;
          xcoord += xshift;

          G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
          radius += (pmtMajorAxis / 2.0);
          if (radius >= (h2oLength-teflonThickness)/2.0)
            continue;

          G4ThreeVector thisPMTPlacement;
          if (pmtDiameter > 0.0) {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness);
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,(h2oHeight-teflonThickness)/2.0 - pmtMinorAxis/2.0 - pmtLegLength);
          }
          new G4PVPlacement(rotTop,thisPMTPlacement,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);

          ++iCopy;
        }

      }
    } else {
      G4double yshift = 0.0;
      if (!(numInY%2)) {
        yshift = pmtMajorAxis / 2.0;
      }
      for (int iY=0; iY<numInY; ++iY) {
        G4int ycount = iY - (numInY/2);
        G4double ycoord = static_cast<G4double>(ycount);
        ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
        ycoord += yshift;

        G4int nx = numInX - TMath::Abs(ycount);
        G4double xshift = 0.0;
        if (!(nx%2)) {
          xshift = pmtMajorAxis/2.0;
        }

        for (int iX=0; iX<nx; ++iX) {
          G4int xcount = iX - (nx/2);
          G4double xcoord = static_cast<G4double>(xcount);
          xcoord *= pmtMajorAxis + minGapBetweenPMTs;
          xcoord += xshift;

          G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
          radius += (pmtMajorAxis / 2.0);
          if (radius >= (h2oLength+teflonThickness)/2.0)
            continue;

          G4ThreeVector thisPMTPlacement;
          if (pmtDiameter > 0.0) {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness);
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength);
          }
          new G4PVPlacement(rotTop,thisPMTPlacement,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);

          ++iCopy;
        }
      }
    }

    //x-y plane, -z
    if(iUseBottomPMTs) {
      G4RotationMatrix *rotBottom = new G4RotationMatrix();

      if (!(numInY%2) && !(numInX%2)) {
        G4double yshift = 0.0;
        if (!(numInY%2)) {
          yshift = (pmtMajorAxis / 2.0) / TMath::Sqrt(3.0);
        }
        for (int iY=0; iY<numInY; ++iY) {
          G4int ycount = iY - (numInY/2);
          G4double ycoord = static_cast<G4double>(ycount);
          ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
          ycoord += yshift;

          G4int nx = numInX - TMath::Abs(ycount);
          G4double xshift = 0.0;
          if (!(nx%2)) {
            xshift = pmtMajorAxis/2.0;
          }

          for (int iX=0; iX<nx; ++iX) {
            G4int xcount = iX - (nx/2);
            G4double xcoord = static_cast<G4double>(xcount);
            xcoord *= pmtMajorAxis + minGapBetweenPMTs;
            xcoord += xshift;

            G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
            radius += (pmtMajorAxis / 2.0);
            if (radius >= h2oLength/2.0)
              continue;

            G4ThreeVector thisPMTPlacement;
            if (pmtDiameter > 0.0) {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-teflonThickness));
            } else {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength));
            }
            new G4PVPlacement(rotBottom,thisPMTPlacement,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);

            ++iCopy;
          }

        }
      } else {
        G4double yshift = 0.0;
        if (!(numInY%2)) {
          yshift = pmtMajorAxis / 2.0;
        }
        for (int iY=0; iY<numInY; ++iY) {
          G4int ycount = iY - (numInY/2);
          G4double ycoord = static_cast<G4double>(ycount);
          ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
          ycoord += yshift;

          G4int nx = numInX - TMath::Abs(ycount);
          G4double xshift = 0.0;
          if (!(nx%2)) {
            xshift = pmtMajorAxis/2.0;
          }

          for (int iX=0; iX<nx; ++iX) {
            G4int xcount = iX - (nx/2);
            G4double xcoord = static_cast<G4double>(xcount);
            xcoord *= pmtMajorAxis + minGapBetweenPMTs;
            xcoord += xshift;

            G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
            radius += (pmtMajorAxis / 2.0);
            if (radius >= h2oLength/2.0)
              continue;

            G4ThreeVector thisPMTPlacement;
            if (pmtDiameter > 0.0) {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-teflonThickness));
            } else {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength));
            }
            new G4PVPlacement(rotBottom,thisPMTPlacement,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);

            ++iCopy;
          }
        }
      }
    }
  
    totPMT = iCopy;
    numInRows = totPMT;
    numRows = totPMT;

    // Now store PMT locations...
    iCopy = 1;
    for(G4int iDim=0; iDim<3; iDim++)
        hPanel[iDim] = new TH1D(Form("hPanel[%d]",iDim),Form("Dimension %d",iDim),totPMT,-0.5,totPMT-0.5);
    
    if (!(numInY%2) && !(numInX%2)) {
      G4double yshift = 0.0;
      if (!(numInY%2)) {
        yshift = (pmtMajorAxis / 2.0) / TMath::Sqrt(3.0);
      }
      for (int iY=0; iY<numInY; ++iY) {
        G4int ycount = iY - (numInY/2);
        G4double ycoord = static_cast<G4double>(ycount);
        ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
        ycoord += yshift;

        G4int nx = numInX - TMath::Abs(ycount);
        G4double xshift = 0.0;
        if (!(nx%2)) {
          xshift = pmtMajorAxis/2.0;
        }

        for (int iX=0; iX<nx; ++iX) {
          G4int xcount = iX - (nx/2);
          G4double xcoord = static_cast<G4double>(xcount);
          xcoord *= pmtMajorAxis + minGapBetweenPMTs;
          xcoord += xshift;

          G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
          radius += (pmtMajorAxis / 2.0);
          if (radius >= h2oLength/2.0)
            continue;

          G4ThreeVector thisPMTPlacement;
          if (pmtDiameter > 0.0) {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness);
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength);
          }
          hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
          hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
          if (pmtDiameter > 0.0) {
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0);
          } else {
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0 - pmtLegLength);
          }

          ++iCopy;
        }

      }
    } else {
      G4double yshift = 0.0;
      if (!(numInY%2)) {
        yshift = pmtMajorAxis / 2.0;
      }
      for (int iY=0; iY<numInY; ++iY) {
        G4int ycount = iY - (numInY/2);
        G4double ycoord = static_cast<G4double>(ycount);
        ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
        ycoord += yshift;

        G4int nx = numInX - TMath::Abs(ycount);
        G4double xshift = 0.0;
        if (!(nx%2)) {
          xshift = pmtMajorAxis/2.0;
        }

        for (int iX=0; iX<nx; ++iX) {
          G4int xcount = iX - (nx/2);
          G4double xcoord = static_cast<G4double>(xcount);
          xcoord *= pmtMajorAxis + minGapBetweenPMTs;
          xcoord += xshift;

          G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
          radius += (pmtMajorAxis / 2.0);
          if (radius >= h2oLength/2.0)
            continue;

          G4ThreeVector thisPMTPlacement;
          if (pmtDiameter > 0.0) {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness);
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength);
          }
          hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
          hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
          if (pmtDiameter > 0.0) {
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0);
          } else {
            hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0 - pmtLegLength);
          }

          ++iCopy;
        }
      }
    }

    //x-y plane, -z
    if(iUseBottomPMTs) {
      G4RotationMatrix *rotBottom = new G4RotationMatrix();

      if (!(numInY%2) && !(numInX%2)) {
        G4double yshift = 0.0;
        if (!(numInY%2)) {
          yshift = (pmtMajorAxis / 2.0) / TMath::Sqrt(3.0);
        }
        for (int iY=0; iY<numInY; ++iY) {
          G4int ycount = iY - (numInY/2);
          G4double ycoord = static_cast<G4double>(ycount);
          ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
          ycoord += yshift;

          G4int nx = numInX - TMath::Abs(ycount);
          G4double xshift = 0.0;
          if (!(nx%2)) {
            xshift = pmtMajorAxis/2.0;
          }

          for (int iX=0; iX<nx; ++iX) {
            G4int xcount = iX - (nx/2);
            G4double xcoord = static_cast<G4double>(xcount);
            xcoord *= pmtMajorAxis + minGapBetweenPMTs;
            xcoord += xshift;

            G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
            radius += (pmtMajorAxis / 2.0);
            if (radius >= h2oLength/2.0)
              continue;

            G4ThreeVector thisPMTPlacement;
            if (pmtDiameter > 0.0) {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-teflonThickness));
            } else {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength));
            }
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
            if (pmtDiameter > 0.0) {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0);
            } else {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0-pmtLegLength);
            }

            ++iCopy;
          }

        }
      } else {
        G4double yshift = 0.0;
        if (!(numInY%2)) {
          yshift = pmtMajorAxis / 2.0;
        }
        for (int iY=0; iY<numInY; ++iY) {
          G4int ycount = iY - (numInY/2);
          G4double ycoord = static_cast<G4double>(ycount);
          ycoord *= TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + minGapBetweenPMTs);
          ycoord += yshift;

          G4int nx = numInX - TMath::Abs(ycount);
          G4double xshift = 0.0;
          if (!(nx%2)) {
            xshift = pmtMajorAxis/2.0;
          }

          for (int iX=0; iX<nx; ++iX) {
            G4int xcount = iX - (nx/2);
            G4double xcoord = static_cast<G4double>(xcount);
            xcoord *= pmtMajorAxis + minGapBetweenPMTs;
            xcoord += xshift;

            G4double radius = TMath::Sqrt(TMath::Power(ycoord, 2.0) + TMath::Power(xcoord, 2.0) );
            radius += (pmtMajorAxis / 2.0);
            if (radius >= h2oLength/2.0)
              continue;

            G4ThreeVector thisPMTPlacement;
            if (pmtDiameter > 0.0) {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-teflonThickness));
            } else {
              thisPMTPlacement = G4ThreeVector(xcoord,ycoord,-(h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength));
            }
            //save coordinates of the front faces of the PMTs
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
            if (pmtDiameter > 0.0) {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0);
            } else {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0-pmtLegLength);
            }

            ++iCopy;
          }
        }
      }
    }
    
} //END OF PLACEPMTS()

 
 

  G4LogicalVolume *G4d2oDetector::GetEllipsoidPMT(){
     ///// PMTs /////
    //--------Ellipsoids--------
    
    // Taken from cenns10geant4 code
  G4cout<<"Creating Ellipsoidal PMT\n";
  G4LogicalVolume *pmtLogV = 0;

  // Top ellipsoid of the PMT
  // 76.72mm long, 103.6mm wide
  // Nominal PMT axes are 103.6mm, 76.72mm. Keep that ratio even for alternative pmt diameters
  G4double OPMTTopESemiX      = 103.6*mm;
  G4double OPMTTopESemiY      = OPMTTopESemiX;
  G4double OPMTTopESemiZ      = 76.72*mm;
  G4double OPMTTopEBottomCut  = 0;
  G4double OPMTTopETopCut     = OPMTTopESemiZ;
  G4double PMTFrameThickness  = pmtWindowThickness;

  G4String OPMTTopEVolName    = "PMTTopEVol";
  G4Ellipsoid *OPMTTopESolid = new G4Ellipsoid(OPMTTopEVolName,OPMTTopESemiX, OPMTTopESemiY, OPMTTopESemiZ, OPMTTopEBottomCut, OPMTTopETopCut);

  // Bottom ellipsoid of the PMT
  G4String OPMTBottomEName      = "OPMTBottomEVol";
  G4double OPMTBottomEBottomCut = -OPMTTopETopCut;
  G4double OPMTBottomETopCut    = 0;
  
 
  G4Ellipsoid *OPMTBottomESolid = new G4Ellipsoid(OPMTBottomEName, OPMTTopESemiX, OPMTTopESemiY, OPMTTopESemiZ, OPMTBottomEBottomCut, OPMTBottomETopCut);

  // Outer Spherical PMT Part
  G4String OPMTSPartName = "OuterPMTSphericalVol";
  G4Transform3D Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,0));
  
  G4UnionSolid *OPMTSPartSolid  = new G4UnionSolid(OPMTSPartName, OPMTTopESolid, OPMTBottomESolid, Transform); //joining top with bottom OUTER ellipsoids.

  // Outer PMT leg
  // Should probably have a check here that the leg doesn't extend out both sides ofthe PMT...
  G4double PMTLegPenetration = OPMTTopESemiZ - 60.84*mm;  // 60.84 - radius value in leg's begining poing
  G4double PMTLegProtruding  = 62*mm;
  G4double PMTLegDiameter    = 84.5*mm;

  G4String OPMTLegName        = "OPMTLegVol";
  G4double OPMTLegInnerRad    = 0;
  G4double OPMTLegOuterRad    = PMTLegDiameter/2;
  G4double OPMTLegZsize       = PMTLegPenetration + PMTLegProtruding;
  G4double OPMTLegSphi        = 0;
  G4double OPMTLegPphi        = 2*M_PI*rad;
  G4double OPMTZposition      = OPMTLegZsize/2 + 60.84*mm;             // Center Position
  G4Tubs *OPMTLegSolid = new G4Tubs(OPMTLegName, OPMTLegInnerRad, OPMTLegOuterRad, OPMTLegZsize/2, OPMTLegSphi, OPMTLegPphi);

  
  // Outer PMT Full Shape
  
  Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-OPMTZposition));
  
  G4String OPMTVolName = "pmtLogV";
  G4UnionSolid      *OPMTSolid  = new G4UnionSolid(OPMTVolName, OPMTSPartSolid, OPMTLegSolid, Transform); //Joining OUTER PMT head with leg
  pmtLogV    = new G4LogicalVolume(OPMTSolid, matPtr->GetMaterial( BOROSILICATE ), OPMTVolName); //Creating full PMT logical volume, made of Borosilicate.  

  // Inner PMT shape

  // FIrst the TOP ellipsoid
  G4String IPMTTopEName         = "IPMTtopEVol";
  G4double IPMTTopESemiX        = OPMTTopESemiX - PMTFrameThickness;
  G4double IPMTTopESemiY        = OPMTTopESemiY - PMTFrameThickness;
  G4double IPMTTopESemiZ        = OPMTTopESemiZ - PMTFrameThickness;
  G4double IPMTTopEBottomCut    = 0;
  G4double IPMTTopECut          = IPMTTopESemiZ;
  G4Ellipsoid *IPMTTopESolid = new G4Ellipsoid(IPMTTopEName,IPMTTopESemiX, IPMTTopESemiY, IPMTTopESemiZ, IPMTTopEBottomCut, IPMTTopECut); //Top Inner ellipsoid (VACUUM)
  
  G4LogicalVolume *ITopPMTLog= new G4LogicalVolume(IPMTTopESolid, matPtr->GetMaterial( VACUUM ), "ITopPMTLog"); //This is where the photocathode will be placed into.
  
//----------
  
  
  //---Inner BOTTOM PMT Spherical part solid (silver coating to prevent photons to escape)
  
  // Inner PMT Bottom Ellipsoid
  G4String IPMTBottomEName      = "IPMTbottomSilverEVol";
  G4double IPMTBottomESemiX     = IPMTTopESemiX;
  G4double IPMTBottomESemiY     = IPMTTopESemiY;
  G4double IPMTBottomESemiZ     = IPMTTopESemiZ;
  G4double IPMTBottomEBottomCut = -IPMTBottomESemiZ;
  G4double IPMTBottomETopCut    = 0;
  G4Ellipsoid *IPMTBottomESolid = new G4Ellipsoid(IPMTBottomEName,IPMTBottomESemiX, IPMTBottomESemiY, IPMTBottomESemiZ, IPMTBottomEBottomCut, IPMTBottomETopCut); //Bottom Inner ellipsoid
  
  // Inner PMT Leg
  G4String IPMTLegName      = "IPMTSilverLegVol";
  G4double IPMTLegInnerRad  = 0;
  G4double IPMTLegOuterRad  = OPMTLegOuterRad - PMTFrameThickness;
  G4double IPMTLegZsize     = OPMTLegZsize - PMTFrameThickness; //****** 2*
  G4double IPMTLegSphi      = 0;
  G4double IPMTLegPphi      = 2*M_PI*rad;
  G4double IPMTZposition    = 60.84*mm  + IPMTLegZsize/2;
  G4Tubs *IPMTLegSolid = new G4Tubs(IPMTLegName, IPMTLegInnerRad, IPMTLegOuterRad, IPMTLegZsize/2, IPMTLegSphi, IPMTLegPphi); //Inner leg

  
  
  //Joining Inner Bottom and Leg 
  G4String  IPMTSilverPartName = "InnerPMTSilverVol";
  Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-IPMTZposition));
  
  G4UnionSolid *IPMTSPartSolid = new G4UnionSolid(IPMTSilverPartName, IPMTBottomESolid, IPMTLegSolid, Transform); 


  
  // Inner Bottom SILVER shape:
  
  G4String IPMTSilverVolName = "InnerSilverPMTVol";
  G4LogicalVolume *IPMTSilverLog = new G4LogicalVolume(IPMTSPartSolid, matPtr->GetMaterial(ALUMINUM), IPMTSilverVolName); //Creating logical volume of INNER HEAD and make it vacuum
  
  
  //IPMTSilverLog->SetVisAttributes(G4VisAttributes::Invisible);
  
  //!!new:
  G4VisAttributes *visAttIPMT = new G4VisAttributes();
  visAttIPMT->SetColour(G4Color::Green());
  visAttIPMT->SetForceSolid(true);
  IPMTSilverLog->SetVisAttributes(visAttIPMT);
  
  
  //---BOTTOM Inner made of VACUUM:
  
    // Inner PMT Bottom Ellipsoid
  G4String IPMTBottomEVName      = "IPMTbottomVacuumEVol";
  G4double IPMTBottomEVSemiX     = IPMTTopESemiX - PMTFrameThickness/2.0;
  G4double IPMTBottomEVSemiY     = IPMTTopESemiY - PMTFrameThickness/2.0;
  G4double IPMTBottomEVSemiZ     = IPMTTopESemiZ - PMTFrameThickness/2.0;
  G4double IPMTBottomEVBottomCut = -IPMTBottomEVSemiZ;
  G4double IPMTBottomEVTopCut    = 0;
  G4Ellipsoid *IPMTBottomVESolid = new G4Ellipsoid(IPMTBottomEVName,IPMTBottomEVSemiX, IPMTBottomEVSemiY, IPMTBottomEVSemiZ, IPMTBottomEVBottomCut, IPMTBottomEVTopCut); //Bottom Inner ellipsoid
  
  // Inner PMT Leg
  G4String IPMTLegVName      = "IPMTVacuumLegVol";
  G4double IPMTLegVInnerRad  = 0;
  G4double IPMTLegVOuterRad  = OPMTLegOuterRad - 2.0*PMTFrameThickness;
  G4double IPMTLegVZsize     = OPMTLegZsize - 2.0*PMTFrameThickness;
  G4double IPMTLegSVphi      = 0;
  G4double IPMTLegPVphi      = 2*M_PI*rad;
  G4double IPMTVZposition    = 60.84*mm  + IPMTLegVZsize/2;
  G4Tubs *IPMTVLegSolid = new G4Tubs(IPMTLegVName, IPMTLegVInnerRad, IPMTLegVOuterRad, IPMTLegVZsize/2, IPMTLegSVphi, IPMTLegPVphi); //Inner leg

  
  //Joining Inner Bottom and Leg 
  G4String  IPMTSVPartName = "InnerPMTVacuumVol";
  Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,-IPMTVZposition));
  //Transform = G4Transform3D(*PMTRotation, G4ThreeVector(0,0,0));//*New
  G4UnionSolid *IPMTSVPartSolid = new G4UnionSolid(IPMTSVPartName, IPMTBottomVESolid, IPMTVLegSolid, Transform); 

  
  // Inner BOTTOM VACUUM Shape
  
  G4String IPMTVolName = "InnerVacuumPMTVol";
  
  G4LogicalVolume   *IPMTVacuumLog    = new G4LogicalVolume(IPMTSVPartSolid, matPtr->GetMaterial( VACUUM ), IPMTVolName); //Creating logical volume of INNER HEAD and make it vacuum
  
  //IPMTVacuumLog->SetVisAttributes(G4VisAttributes::Invisible);
  
  //!!new:
  G4VisAttributes *visAttIVPMT = new G4VisAttributes();
  visAttIVPMT->SetColour(G4Color::Magenta());
  visAttIVPMT->SetForceSolid(true);
  IPMTVacuumLog->SetVisAttributes(visAttIVPMT);
  
  
  
  
  //--- PLACE VACUUM BOTTOM PART INTO SILVER BOTTOM PART 
  
  G4String  BottomPMTPartName = "BottomPart";
  
  G4VPhysicalVolume *BottomPMT  = 
  new G4PVPlacement(0, G4ThreeVector(0,0,0), IPMTVacuumLog, BottomPMTPartName, IPMTSilverLog, false, 0, true); // "Vacuum in Silver"
  
  // SILVER IN OUTER pmtLogV
  G4VPhysicalVolume *AgCoating = 
  new G4PVPlacement(0, G4ThreeVector(0,0,0), IPMTSilverLog, "BottomTotal", pmtLogV, false, 0, true);
  
  
//////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////
//   Optical properties of the interface between the Air and Reflective Surface
//   For Mirror, reflectivity is set at 100% and specular reflection is assumed.
  
     
    G4OpticalSurface* OpticalAirMirror = new G4OpticalSurface("OpticalAirMirror");
    
    new G4LogicalBorderSurface("OpticalAirMirror", AgCoating , BottomPMT, OpticalAirMirror );
    
    OpticalAirMirror->SetModel(unified);
    OpticalAirMirror->SetType(dielectric_dielectric);
    OpticalAirMirror->SetFinish(polishedfrontpainted);
    
    const G4int NUM = 2;
    G4double XX[NUM] = {0.1*MeV, 55.0*MeV} ; 
    G4double ICEREFLECTIVITY[NUM]  = { 1.0, 1.0 };
    
     G4MaterialPropertiesTable *AirMirrorMPT = new G4MaterialPropertiesTable();
     AirMirrorMPT->AddProperty("REFLECTIVITY", XX, ICEREFLECTIVITY,NUM);
     OpticalAirMirror->SetMaterialPropertiesTable(AirMirrorMPT);
  
    //matPtr->SetReflector(BottomPMT, AgCoating, PMTcoatingReflectivity, PMTcoatingSigmaAlpha);
     
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////     
     
  
  // TOP PART IN OUTER pmtLogV
  
  //G4VPhysicalVolume *TopPMT = 
  new G4PVPlacement(0, G4ThreeVector(0,0,0), ITopPMTLog, "TopTotal", pmtLogV, false, 0, true);
  
  
  // Then place inner PMT volume inside outer PMT volume
  //G4VPhysicalVolume *IPMTPhys  = new G4PVPlacement(0, G4ThreeVector(0,0,0), IPMTLog, IPMTVolName, pmtLogV, false, 0, true); //top
 

  // Now place the photocathode
  
  G4String PMTPhtVolName    = "PhotocathodeVol";
  G4double PMTPhtThickness  =  2.6E-5*mm;
  G4double PMTPhtSRad       =  131*mm;                                                // Photocathode Cons Surface Side
  G4double PMTPhtSurDiam    =  190*mm;                                                // Photocathode surface diameter
  G4double PMTPhtAngle      =  acos(PMTPhtSurDiam/(2*PMTPhtSRad))*rad;
  G4double PMTPhtOuterRad   =  PMTPhtSRad-PMTFrameThickness;
  G4double PMTPhtInnerRad   =  PMTPhtSRad - PMTFrameThickness - PMTPhtThickness;
//  G4double PMTPhtZpozition  =  PMTPhtSRad - OPMTTopESemiZ + PMTFrameThickness + PMTPhtThickness/2;
  G4double PMTPhtZpozition  =  - PMTPhtSRad + OPMTTopESemiZ - PMTFrameThickness;
  
  G4Sphere          *PMTPhtSolid  = new G4Sphere(PMTPhtVolName, PMTPhtInnerRad, PMTPhtOuterRad, 0, 2*M_PI*rad, 0, PMTPhtAngle);
  G4LogicalVolume   *PMTPhtLog    = new G4LogicalVolume(PMTPhtSolid, matPtr->GetMaterial( PHOTOCATHODE ), PMTPhtVolName);
  
  G4String PMTPhtTopVolName = "PhotocathodeVol";
  G4VPhysicalVolume *PMTPhtPhys  = new G4PVPlacement(0, G4ThreeVector(0, 0, PMTPhtZpozition), PMTPhtLog, PMTPhtTopVolName, ITopPMTLog, false, 0,  true);  //top
  

  ////// Make it sensitive ////////
  G4cout<<G4endl;
  G4d2oSensitiveDetector *senDet = new G4d2oSensitiveDetector("pmt", PMTPhtLog->GetMaterial()->GetName(), 0);
  G4SDManager *sdManager = G4SDManager::GetSDMpointer(); //New
  sdManager->AddNewDetector(senDet);
  PMTPhtLog->SetSensitiveDetector(senDet);
  G4cout<<G4endl;

  G4VisAttributes *visAttPhotoC = new G4VisAttributes();
  visAttPhotoC->SetColour(G4Color::Yellow());
  visAttPhotoC->SetForceSolid(true);
  
  G4VisAttributes *visAttPhotoM = new G4VisAttributes();
  visAttPhotoM->SetColour(G4Color::White());
  visAttPhotoM->SetForceSolid(true);
  
  PMTPhtLog->SetVisAttributes( visAttPhotoC );
  pmtLogV->SetVisAttributes( visAttPhotoM );

  return pmtLogV;
}
//*


