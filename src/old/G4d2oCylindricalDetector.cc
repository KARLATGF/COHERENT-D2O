// PMTs ON TOP GEOMETRY!!

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

using namespace std;

G4d2oCylindricalDetector::G4d2oCylindricalDetector()
{
    G4cerr << "\n\tUsing cylindrical geometry" << G4endl ;
    
    Initialize();
}

void G4d2oCylindricalDetector::Initialize(){
    
    //relevant parameters
    
    //D2O
    d2oLength = 70.0*cm;//radius
    d2oWidth = d2oLength; //radius
    //d2oHeight = 140.0*cm;//z
    d2oHeight = 161.0*cm;//z
    
    //acrylic
    acrylicThickness = 0.25*in;
    acrylicEndCapThickness = 1.0*in;
    
    //pmt
    pmtDepth = 10.0*in;
    muMetalThickness = 0.03*in;
    pmtWindowThickness = 3.0*mm;  //making this up
    minGapBetweenPMTs = 3.0*mm;
    pmtDiameter = input->GetPMTDiameter()*in; //pmtDiameter controlled from beamOn.dat or command line
    pmtMinorAxis = pmtMajorAxis = pmtDiameter;
    pmtLegLength = 62.*mm;
    if (pmtDiameter < 0.0) {
      pmtMinorAxis = 76.72*mm * 2.0;
      pmtMajorAxis = 103.6*mm * 2.0;
    }

    //teflon
    teflonThickness = 0.25*in;
    //teflonReflectivity = 0.99; //need to verify
    teflonReflectivity = input->GetReflectivity(); //need to verify
    teflonSigmaAlpha = 0.1; //roughness parameter. 0.1 for teflon, use 0 for polished
    
    //H2O
    h2oTCthickness = input->GetTailCatcherThick()*cm;
    h2oLength = d2oLength + 2*(acrylicThickness + h2oTCthickness); //x
    h2oWidth = h2oLength; //y
    if (pmtDiameter > 0.0) {
      h2oHeight = d2oHeight + 2*(acrylicEndCapThickness + h2oTCthickness + pmtMinorAxis/2.0); //z
    } else {
      h2oHeight = d2oHeight + 2*(acrylicEndCapThickness + h2oTCthickness + pmtMinorAxis/2.0 + pmtLegLength/2.0);
     
    }
    
    //not using these for now, but set them anyway so that other parts of code don't have problems
    //These are used to throw random positions for the electron gun to include the outer h2o.
    //The vol0 variable in the event output selects 1 for d2o and 2 for h2o
    h2oInnerLength = h2oLength;
    h2oInnerWidth = h2oWidth;
    h2oInnerHeight = h2oHeight;
    
    //outer vessel
    outerContainerThickness = 0.25*in;
    
    ///// Full detector volume /////
    contOuterLength = h2oLength + 2*outerContainerThickness; //x
    contOuterWidth = contOuterLength; //y
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

    //bUseBottomPMTs = true;
    iUseBottomPMTs = input->GetBottomPMTs(); //PMTs on the bottom or just teflon reflector

    //DetermineSpacing must be called after all the parameters above are defined
    DetermineSpacing();
    
    G4cerr << "done." << G4endl;

}//END of constructor

G4d2oCylindricalDetector::~G4d2oCylindricalDetector()
{
    G4cout << "Instance of G4d2oCylindricalDetector Destructed!" << G4endl;
    
}//END of destructor

void G4d2oCylindricalDetector::DetermineSpacing(){
  spacingInX = pmtMajorAxis + (2.0*minGapBetweenPMTs);
  numInX = TMath::FloorNint( (h2oLength / spacingInX) );

  spacingInY = TMath::Sqrt(3.0) * ( (pmtMajorAxis/2.0) + 2.0*minGapBetweenPMTs);
  numInY = TMath::FloorNint( (h2oLength / spacingInY) );
  //numInY = numInX;

/*    
    //minimum gap between PMTs and between PMTs and edge of volume
    G4double nominalSpacing = pmtMajorAxis + minGapBetweenPMTs;
    
    G4double totalRadius = (h2oLength-nominalSpacing)/2.0;
    
    numInX = TMath::FloorNint(totalRadius/nominalSpacing)+1;
    G4cout<<h2oLength<<"  "<<nominalSpacing<<"  "<<totalRadius<<"  "<<numInX<<"  "<<pmtMajorAxis<<G4endl;
    spacingInX = nominalSpacing;

    spacingInY = 2.0*TMath::ASin(0.5);
    numInY = TMath::FloorNint(TMath::TwoPi()/spacingInY);
    
    numInZ = 0;
    spacingInZ = 0.0;
    
    numRows = numInX;
    if(iUseBottomPMTs) numRows *= 2;
    numInRows = numInY;
    totPMT = numRows*numInRows;
*/        
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
    new G4PVPlacement(0,G4ThreeVector(0,0,0),d2oLogV,"d2oPhysV",acrylicLogV,false,0,true); // FOR H2O ONLY, COMMENT THIS OUT
    
    ///// Place the acrylic tank in H2O /////
    new G4PVPlacement(0,G4ThreeVector(0,0,-(pmtLegLength+pmtMinorAxis)/2.0),acrylicLogV,"acrylicPhysV",h2oLogV,false,0,true); // added -(pmtLegLength+pmtMinorAxis)/2.0!!!! // FOR H2O ONLY, COMMENT THIS OUT
    
    ///// Place PMTs into H2O /////
    PlacePMTs(pmtLogV,h2oLogV); 

    //// Create teflon end caps ////
    G4LogicalVolume *teflonCapTop = GetTeflonLiningCapLogV("topCap", true);
    G4LogicalVolume *teflonCapBot = GetTeflonLiningCapLogV("botCap", iUseBottomPMTs);
    
    ///// Place the H2O in teflon lining /////
    //G4VPhysicalVolume *h2oPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),h2oLogV,"h2oPhysV",teflonLining,false,0,true);
    
    G4VPhysicalVolume *h2oPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),h2oLogV,"h2oPhysV",outerVesselLogV,false,0,true); //bring this back when pertinent!!!
    //G4VPhysicalVolume *h2oPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,-(0.5*pmtLegLength+0.5*pmtMinorAxis)/2.0),h2oLogV,"h2oPhysV",outerVesselLogV,false,0,true); // FOR H2O ONLY!
    
    ///// Place the teflon lining in outer vessel/////
    //G4VPhysicalVolume *teflonPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),teflonLining,"teflonPhysV",outerVesselLogV,false,0,true); // FOR H2O ONLY, COMMENT THIS OUT
    
    G4VPhysicalVolume *teflonPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,-(pmtLegLength+pmtMinorAxis-h2oTCthickness+2.0*teflonThickness)/2.0),teflonLining,"teflonPhysV",h2oLogV,false,0,true); // added -(pmtLegLength+pmtMinorAxis-h2oTCthickness+teflonThickness)/2.0!!!! //bring this back when pertinent!!!
    //G4VPhysicalVolume *teflonPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,-(pmtLegLength-0.5*pmtMinorAxis-h2oTCthickness+2.0*teflonThickness)/2.0),teflonLining,"teflonPhysV",h2oLogV,false,0,true); // FOR H2O ONLY, COMMENT THIS IN
    
    //set reflectivity of teflon
    matPtr->SetReflector(h2oPhysV, teflonPhysV, teflonReflectivity, teflonSigmaAlpha);

    // Place the teflon caps
    G4VPhysicalVolume *teflonCapTopPhysV = 0;
    G4VPhysicalVolume *teflonCapBotPhysV = 0;
    if (teflonCapBot) {
      if (pmtDiameter < 0.0) 
        teflonCapBotPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -(h2oHeight-teflonThickness)/2.0 ), teflonCapBot, "teflonCapBotPhysV", h2oLogV, false, 0, true);//Bring back when prudent!!!!!
          //teflonCapBotPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -(h2oHeight-teflonThickness-0.5*pmtLegLength-0.5*pmtMinorAxis)/2.0 ), teflonCapBot, "teflonCapBotPhysV", h2oLogV, false, 0, true); // FOR H2O ONLY, COMMENT THIS IN
      else 
        teflonCapBotPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -(h2oHeight-teflonThickness)/2.0), teflonCapBot, "teflonCapBotPhysV", h2oLogV, false, 0, true);//Bring back when prudent!!!!!
          //teflonCapBotPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -(h2oHeight-teflonThickness-0.5*pmtLegLength-0.5*pmtMinorAxis)/2.0), teflonCapBot, "teflonCapBotPhysV", h2oLogV, false, 0, true); // FOR H2O ONLY, COMMENT THIS IN
      matPtr->SetReflector(h2oPhysV, teflonCapBotPhysV, teflonReflectivity, teflonSigmaAlpha);
    }
    
    
    if (teflonCapTop) { // Tyvek Top Cap, bring back when prudent!!!
      if (pmtDiameter < 0.0)
        teflonCapTopPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, (h2oHeight-pmtMinorAxis)/2.0 - pmtLegLength), teflonCapTop, "teflonCapTopPhysV", h2oLogV, false, 0, true);
      else 
        teflonCapTopPhysV = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, (h2oHeight-teflonThickness)/2.0), teflonCapTop, "teflonCapTopPhysV", h2oLogV, false, 0, true);
      matPtr->SetReflector(h2oPhysV, teflonCapTopPhysV, teflonReflectivity, teflonSigmaAlpha);
    } // FOR H2O ONLY, COMMENT THIS OUT
    
    ///// Place the outer vessel inside the shielding
   /* // Muon Veto, bring back when prudent!!!
    G4double offset = 0.0;
    if (!input->GetBottomShielding() )
      offset = -0.5*shieldThickness;
    G4VPhysicalVolume *outerVesselPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,offset),outerVesselLogV,"outerVesselPhysV",shieldingLogV,false,0,true);

    ///// Place the Pb shielding inside the veto
    offset = 0.0;
    if (!input->GetBottomVeto() ) 
      offset = -0.5*muonVetoThickness;
    //G4VPhysicalVolume *shieldingPhysV = new G4PVPlacement(0, G4ThreeVector(0,0,offset), shieldingLogV, "shieldingPhysV", innerVetoLogV, false, 0, true);
    
    /// // Place the inner veto inside the outer veto layer
    G4VPhysicalVolume *outerVetoPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,offset),innerVetoLogV,"innerVetoPhysV",vetoLogV,false,0,true);

    ///// Place the outer vessel in total detector/////
    new G4PVPlacement(0,G4ThreeVector(0,0,0),vetoLogV,"vetoPhysV",totalDetLogV,false,0,true); */
    new G4PVPlacement(0,G4ThreeVector(0,0,0),outerVesselLogV,"outerVesselLogV",totalDetLogV,false,0,true);
    
    // Add Event Handlers for Muon Vetos
    // Muon Vetos Inner
   /* // Muon Veto Handlers, bring back when prudent!!!
    G4d2oSensitiveDetector *senDetmuVetoInner = new G4d2oSensitiveDetector("muVetoInner", innerVetoLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoInner);
    innerVetoLogV->SetSensitiveDetector(senDetmuVetoInner);
    G4d2oSensitiveDetector *senDetmuVetoOuter = new G4d2oSensitiveDetector("muVetoOuter", vetoLogV->GetMaterial()->GetName(), 0);
    sdManager->AddNewDetector(senDetmuVetoOuter);
    vetoLogV->SetSensitiveDetector(senDetmuVetoOuter); */

    return totalDetLogV;
    
    
}//END of GetDetector()

G4LogicalVolume *G4d2oCylindricalDetector::GetTotalDetectorLogV(){
    
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
    G4Tubs* outerVessel = new G4Tubs("outerVessel",0.0,contOuterLength/2.0,contOuterHeight/2.0,0.0,360.0*deg);
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
    
    //G4Tubs* h2oSolid = new G4Tubs("h2oSolid",0.0,(h2oLength-teflonThickness)/2.0,
    //                              (h2oHeight-teflonThickness)/2.0,0.0,360.0*deg);
    
    G4Tubs* h2oSolid = new G4Tubs("h2oSolid",0.0,(h2oLength)/2.0,(h2oHeight)/2.0,0.0,360.0*deg); // Bring this back this when appropiate!!!
    
    //G4Tubs* h2oSolid = new G4Tubs("h2oSolid",0.0,(h2oLength)/2.0,(h2oHeight-0.5*pmtLegLength-0.5*pmtMinorAxis)/2.0,0.0,360.0*deg); // FOR H2O ONLY, COMMENT THIS BACK IN
    
    
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
  //d2oHeight + 2*(acrylicThickness + h2oTCthickness + pmtMinorAxis + pmtLegLength);
    //teflonThickness = 2.55*in;
    
    G4double teflonheight = (d2oHeight/2.0) + acrylicEndCapThickness + h2oTCthickness + (pmtMinorAxis/2.0) - (pmtLegLength+3.0*teflonThickness)/2.0; // Bring back when prudent!!! //added - (pmtLegLength+2.0*teflonThickness)/2.0 !!!!
    //G4double teflonheight = (d2oHeight/2.0) + acrylicEndCapThickness + h2oTCthickness + (pmtMinorAxis/2.0) - (pmtLegLength+3.0*teflonThickness)/2.0 + pmtLegLength; // FOR H2O ONLY, COMMENT THIS BACK IN!
    
    ////G4double teflonheight = (d2oHeight/2.0) + acrylicEndCapThickness + h2oTCthickness + (pmtMinorAxis+pmtLegLength)/2.0;
    
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

void G4d2oCylindricalDetector::PlacePMTs(G4LogicalVolume *thePMTLogV, G4LogicalVolume *theMotherLogV){
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
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness); //Bring back when prudent!!!
            //thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness+pmtLegLength); // FOR H2O ONLY, COMMENT THIS BACK IN
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,(h2oHeight-teflonThickness)/2.0 - pmtMinorAxis/2.0 - pmtLegLength); //Bring back when prudent!!!
            //thisPMTPlacement = G4ThreeVector(xcoord,ycoord,(h2oHeight-teflonThickness)/2.0 - pmtMinorAxis/2.0 - pmtLegLength+pmtLegLength); // FOR H2O ONLY, COMMENT THIS BACK IN
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
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness); //Bring back when prudent!!!
            //thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness+pmtLegLength); // FOR H2O ONLY, COMMENT THIS BACK IN
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength); //Bring back when prudent!!!
            //thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength+pmtLegLength); // FOR H2O ONLY, COMMENT THIS BACK IN
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
/*
    for(G4int iR=0; iR<numInX; iR++){
        
        G4double numInThisRow = TMath::FloorNint(numInY*iR);
        G4double angularSpacingThisRow = TMath::TwoPi()/double(numInThisRow);
        
        if(iR==0) numInThisRow = 1;
        
        for(G4int iCircumference=0; iCircumference<numInThisRow; iCircumference++){
            
            G4ThreeVector thisPMTPlacement;
            G4ThreeVector thisPMTHousingPlacement;
            G4double theX, theY;
            if(iR==0){
                theX = 0.0;
                theY = 0.0;
            }
            else{
                TVector2 theLoc(iR*spacingInX,0.0);
                TVector2 theLocP = theLoc.Rotate(iCircumference*angularSpacingThisRow);
                theX = theLocP.X();
                theY = theLocP.Y();
            }
            
            if (pmtDiameter > 0.0) {
              thisPMTPlacement = G4ThreeVector(theX,theY,h2oHeight/2.0-teflonThickness);
            } else {
              thisPMTPlacement = G4ThreeVector(theX,theY,h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength);
            }
            new G4PVPlacement(rotTop,thisPMTPlacement,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);
            
            iCopy++;
            
            //save coordinates of the front faces of the PMTs
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
            if (pmtDiameter > 0.0) {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0);
            } else {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0 - pmtLegLength);
            }
        }
    }
    
    //x-y plane, -z
    G4RotationMatrix *rotBottom = new G4RotationMatrix();
    
    for(G4int iR=0; iR<numInX && iUseBottomPMTs; iR++){
        
        G4double numInThisRow = TMath::FloorNint(numInY*iR);
        G4double angularSpacingThisRow = TMath::TwoPi()/double(numInThisRow);
        
        if(iR==0) numInThisRow = 1;
        
        for(G4int iCircumference=0; iCircumference<numInThisRow; iCircumference++){
            
            G4ThreeVector thisPMTPlacement;
            G4ThreeVector thisPMTHousingPlacement;
            G4double theX, theY;
            if(iR==0){
                theX = 0.0;
                theY = 0.0;
            }
            else{
                TVector2 theLoc(iR*spacingInX,0.0);
                TVector2 theLocP = theLoc.Rotate(iCircumference*angularSpacingThisRow);
                theX = theLocP.X();
                theY = theLocP.Y();
            }
            
            if (pmtDiameter > 0.0) {
              thisPMTPlacement = G4ThreeVector(theX,theY,-(h2oHeight/2.0-teflonThickness));
            } else {
              thisPMTPlacement = G4ThreeVector(theX,theY,-(h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength));
            }
            new G4PVPlacement(rotBottom,thisPMTPlacement,thePMTLogV,Form("pmtPhysV_%d",iCopy),theMotherLogV,true,iCopy,true);
            
            iCopy++;
            
            //save coordinates of the front faces of the PMTs
            hPanel[0]->SetBinContent(iCopy,thisPMTPlacement.x());
            hPanel[1]->SetBinContent(iCopy,thisPMTPlacement.y());
            if (pmtDiameter > 0.0) {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0);
            } else {
              hPanel[2]->SetBinContent(iCopy,thisPMTPlacement.z()-pmtMinorAxis/2.0-pmtLegLength);
            }
        }
    }
*/    
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
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness); //Bring this back when prudent!!!
            //thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness+pmtLegLength); // FOR H2O ONLY, COMMENT THIS BACK IN
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength); //Bring this back when prudent!!!
            //thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength+pmtLegLength); // FOR H2O ONLY, COMMENT THIS BACK IN
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
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness);//Bring this back when prudent!!!
            //thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness+pmtLegLength); // FOR H2O ONLY, COMMENT THIS BACK IN
          } else {
            thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength);//Bring this back when prudent!!!
            //thisPMTPlacement = G4ThreeVector(xcoord,ycoord,h2oHeight/2.0-teflonThickness - pmtMinorAxis/2.0 - pmtLegLength+pmtLegLength); // FOR H2O ONLY, COMMENT THIS BACK IN
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
    
} //END OF PlacePMTs




