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

G4d2oDetector::~G4d2oDetector(){
    G4cout << "Instance of G4d2oDetector Destructed!" << G4endl;
    
}//END of destructor

void G4d2oDetector::Initialize(){
  //relevant parameters
  
  //D2O

  d2oLength = 69.99*cm;
  d2oWidth = d2oLength; //diameter
//   d2oHeight = 140.0*cm;//z //Considers PMTs on bottom
  d2oHeight = 161.544*cm;//z //Using bottom PMT space for d2o 
  
  //acrylic
  acrylicThickness = 0.25*in;
  
  //pmt
  //pmtDepth = 10.0*in;
  pmtDepth = 21.55*cm; //**New
  muMetalThickness = 0.03*in;
  pmtWindowThickness = 3.0*mm;  //making this up
  //minGapBetweenPMTs = 1.0*mm;
  minGapBetweenPMTs = 3.5*cm; //*New

  input = inputVariables::GetIVPointer();
  pmtDiameter = input->GetPMTDiameter()*in; //pmtDiameter controlled from beamOn.dat or command line. 
  
     //*New:
   pmtMinorAxis = pmtMajorAxis = pmtDiameter;
     pmtLegLength = 62.*mm;
     if (pmtDiameter < 0.0) {
       pmtMinorAxis = 76.72*mm * 2.0;
       pmtMajorAxis = 103.6*mm * 2.0;
     }
   //*

  
  //teflon
  teflonThickness = 0.2*cm;
  teflonReflectivity = 0.97; //need to verify
  teflonSigmaAlpha = 0.1; //roughness parameter. use 0 for polished
  
  //aluminum/silver PMT coating:
  PMTcoatingReflectivity = 1.0;
  PMTcoatingSigmaAlpha = 0.0;

  
  //H2O
   double h2oTCthickness = 10*cm;
//   double h2oTCthickness = input->GetTailCatcherThick()*cm; //New. TC thickness controlled in BeamOn
  h2oWidth = d2oLength + 2*(acrylicThickness + h2oTCthickness); //y 
  h2oLength = d2oLength + 2*(acrylicThickness + h2oTCthickness + pmtMinorAxis +pmtLegLength + 1.5*cm); //x
  //h2oLength = d2oLength + 2*(acrylicThickness + h2oTCthickness + 8.0*cm);
  
    //*New:
  if (pmtDiameter > 0.0) {
      h2oHeight = d2oHeight + 2*(acrylicThickness + h2oTCthickness); //z
    } else {
      h2oHeight = d2oHeight + 2*(acrylicThickness + h2oTCthickness); //z
    }
  //*
  
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
   contOuterWidth = h2oLength + 2*outerContainerThickness; //y?contOuterWidth = h2oLength + 2*(outerContainerThickness) + pmtMinorAxis + 2.0*pmtLegLength; //y?
  contOuterLength= h2oWidth + 2*outerContainerThickness; //x?
  contOuterHeight = h2oHeight + 2*outerContainerThickness; //z
  
  
  //Teflon Shell parameters for panels BEHIND PMTs:
  /*teflonInnerR = (h2oLength - 2.54*cm)/2.0;
  teflonOuterR = (h2oLength - 2.54*cm + 2.0*teflonThickness)/2.0;
  teflonHeight = (h2oHeight - 2.54*cm)/2.0;*/
  
  //Teflon Shell parameters:
  teflonInnerR = (h2oLength-pmtLegLength - pmtMajorAxis - 7.5*cm)/2.0;
  teflonOuterR = (h2oLength-pmtLegLength - pmtMajorAxis - 7.5*cm   + 2.0*teflonThickness)/2.0;
  teflonHeight = (h2oHeight - 1.0*in)/2.0;
  teflonLength = (TMath::Sqrt(pow( (h2oLength-pmtMinorAxis-pmtLegLength - 7.5*cm),2) - pow((h2oLength/2.0),2) ))/3.415; //for front/back panels
  
  
//   bUseBottomPMTs = false;

  std::cout << "All parameters set!\n";

  //DetermineSpacing must be called after all the parameters above are defined
  
  DetermineSpacing();
  std::cout << "PMT spacing determined\n";

  G4cerr << "done." << G4endl;
  
}//END of constructor

void G4d2oDetector::DetermineSpacing(){

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
    
    //here we fix total width to use in calculation b/c we don't want to overlap the edge,
    //thus -2*(pmtMajorAxis-nominalRowSpacing)
    
    G4double totalWidthInZ = h2oHeight  - 2*minGapBetweenPMTs - 2*(pmtMajorAxis-nominalRowSpacing);
    numInZ = TMath::FloorNint(totalWidthInZ/nominalRowSpacing);
    spacingInZ = totalWidthInZ/double(numInZ);
    
    //some totals that get passed to other classes
    numRows = numInZ;
    numInRows = numInCircle;
    totPMT = numRows*numInRows;
    
    sideLiningX = h2oInnerRadius; //not using for now
    sideLiningZ = h2oInnerRadius; //not using for now
      
      //*
        
}

G4LogicalVolume * G4d2oDetector::GetDetector(){
  std::cout << "Getting Detector\n";

    //Materials pointer
    matPtr = G4d2oRunAction::GetMaterialsPointer();
    
    ///// create total detector volume /////
    G4LogicalVolume *totalDetLogV = GetTotalDetectorLogV();
    
    ///// create steel vessel /////
    G4LogicalVolume *outerVesselLogV = GetOuterVesselLogV();
    
    ///// H2O Inner /////
    G4LogicalVolume *h2oLogV = GetH2OLogV();
    
    ///// Create Teflon Lining /////
    
     // Outer teflon:
        //G4LogicalVolume *teflonLining = GetTeflonLiningLogV();
    
     // Inner teflon:
        
    //Left/Right Panel
        G4LogicalVolume *teflonLiningL = GetTeflonLiningLLogV();

     // Front/Back Panel
        G4LogicalVolume *teflonLiningF = GetTeflonLiningFLogV();

     //Top/Bottom Panel
        G4LogicalVolume *teflonLiningT = GetTeflonLiningTLogV();

        
    
    ///// PMMA /////
    G4LogicalVolume *acrylicLogV = GetAcrylicLogV();
    
    ///// D2O /////
    G4LogicalVolume *d2oLogV = GetD2OLogV();
    
    ///// PMTs with mu-metal shield /////

    //*New:
    G4LogicalVolume *pmtLogV = 0;
    if (pmtDiameter < 1.0) {
      pmtLogV = GetEllipsoidPMT();
    } else {
        //pmtLogV = GetAreaPMT(); //For spherical PMTs
        pmtLogV = GetEllipsoidPMT();
    }
    //*
    
    //     OLD Geometry hierarchy
    //     totalDetLogV (mother volume)
    //        -> outer steel
    //           -> teflon lining
    //              -> H2O tank
    //                 -(1)-> PMTs
    //                 -(2)-> acrylic tank
    //                        -> D2O
    
    //     NEW Geometry hierarchy
    //     totalDetLogV (mother volume)
    //        -> outer steel
    //           
    //           -> H2O tank
    //               -(0)-> teflon lining
    //               -(1)-> PMTs
    //               -(2)-> acrylic tank
    //                        -> D2O
    
    ///// 1.- Place the D2O in acrylic tank ///// 
    new G4PVPlacement(0,G4ThreeVector(0,0,0),d2oLogV,"d2oPhysV",acrylicLogV,false,0,true);
    
    ///// 2.- Place the acrylic tank in H2O /////
    new G4PVPlacement(0,G4ThreeVector(0,0,0),acrylicLogV,"acrylicPhysV",h2oLogV,false,0,true);
    
    ///// 3.- Place PMTs into H2O /////
    PlacePMTs(pmtLogV,h2oLogV);
    
    
    /////-- 5.- Place the H2O in outer outer teflon/////
    G4VPhysicalVolume *h2oPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),h2oLogV,"h2oPhysV",outerVesselLogV,false,0,true);
    //G4VPhysicalVolume *h2oPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),h2oLogV,"h2oPhysV",teflonLining,false,0,true);
    
    /////-- 4.- Place the teflon linings in (outer vessel) H2O/////
    
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
        
    ///// Place the teflon lining in outer vessel/////
    //G4VPhysicalVolume *teflonPhysV = new G4PVPlacement(0,G4ThreeVector(0,0,0),teflonLining,"teflonPhysV",outerVesselLogV,false,0,true);
    //set reflectivity of teflon
    //matPtr->SetReflector(h2oPhysV, teflonPhysV, teflonReflectivity, teflonSigmaAlpha);
    
    
    ///// 7.- Place the outer vessel in total detector/////
    new G4PVPlacement(0,G4ThreeVector(0,0,0),outerVesselLogV,"outerVesselPhysV",totalDetLogV,false,0,true);
    
    return totalDetLogV;
    
    
}//END of GetDetector()

G4LogicalVolume *G4d2oDetector::GetTotalDetectorLogV(){
    
    //define the thickness of the side lining
    G4Box* outerSolid = new G4Box("outerSolid",contOuterLength/2.0,contOuterWidth/2.0,contOuterHeight/2.0);
    G4LogicalVolume *totalDetLogV = new G4LogicalVolume(outerSolid,
                                                        matPtr->GetMaterial( AIR ),
                                                        "totalDetLogV");

    return totalDetLogV;
    
}

G4LogicalVolume *G4d2oDetector::GetOuterVesselLogV(){
    
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

G4LogicalVolume *G4d2oDetector::GetH2OLogV(){
    
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

G4LogicalVolume *G4d2oDetector::GetAcrylicLogV(){
    
    G4double acrylicLength = d2oLength + 2.0*acrylicThickness; //x
    G4double acrylicHeight = d2oHeight + 2.0*acrylicThickness; //z
    
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


//-----------Teflon Panels: needed to separate to avoid overlap issues---------------
//------------------------------------------------------------------------------------

/*G4LogicalVolume *G4d2oDetector::GetTeflonLiningLogV(){
    
    //G4Box* teflonSolid = new G4Box("teflonSolid",h2oLength/2.0,h2oWidth/2.0,h2oHeight/2.0);
    G4Box* teflonSolid = new G4Box("teflonSolid",h2oWidth/2.0,h2oLength/2.0,h2oHeight/2.0);
    G4LogicalVolume *theSideLiningLogV = new G4LogicalVolume(teflonSolid,
                                                             matPtr->GetMaterial( TEFLON ),
                                                             "theSideLiningLogV");
    
    G4VisAttributes *visAttTeflon = new G4VisAttributes();
    visAttTeflon->SetColour(G4Color::Red());
    visAttTeflon->SetForceSolid(true);
    theSideLiningLogV->SetVisAttributes( visAttTeflon );
    
    
    return theSideLiningLogV;
    
}*/

//-----------Left/Right panel (with holes)

G4LogicalVolume *G4d2oDetector::GetTeflonLiningLLogV(){
    
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

G4LogicalVolume *G4d2oDetector::GetTeflonLiningFLogV(){
    

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

G4LogicalVolume *G4d2oDetector::GetTeflonLiningTLogV(){
    
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


//*NEW:

void G4d2oDetector::PlacePMTs(G4LogicalVolume *thePMTLogV, G4LogicalVolume *theMotherLogV){


// ****************************** From Cylindrical Geometry: ********************************//

    
    //--------------------------------------------------------------------------------------
  
      //Mu-metal housing
    //at present, this is included to ensure that photons don't enter PMT volume from behind

   /* G4double pmtHousingDepth = pmtDepth - 76.0*mm;
   
    G4Tubs *pmtHousingSolid = new G4Tubs("pmtHousingSolid",0.0,103.6*mm,
                                         pmtHousingDepth/2.0,0.0,359.999*deg);
    G4LogicalVolume *pmtHousingLogV = new G4LogicalVolume(pmtHousingSolid,
                                                          matPtr->GetMaterial( VACUUM ),
                                                          "pmtHousingLogV");
    G4VisAttributes *visAttPhotoHousing = new G4VisAttributes();
    visAttPhotoHousing->SetColour(G4Color::Gray());
    pmtHousingLogV->SetVisAttributes( visAttPhotoHousing );
    
    

    G4Tubs *pmtMuMetalSolid = new G4Tubs("pmtMuMetalSolid",103.6*mm - muMetalThickness,103.6*mm,
                                         pmtHousingDepth/2.0,0.0,359.999*deg);
    G4LogicalVolume *pmtMuMetalLogV = new G4LogicalVolume(pmtMuMetalSolid,
                                                          matPtr->GetMaterial( MUMETAL ),
                                                          "pmtMuMetalLogV");
    G4VisAttributes *visAttMuMetal = new G4VisAttributes();
    visAttMuMetal ->SetColour(G4Color::White());
    pmtMuMetalLogV -> SetVisAttributes( visAttMuMetal );
    new G4PVPlacement(0,G4ThreeVector(0,0,0),pmtMuMetalLogV,"pmtMuMetalPhysV",pmtHousingLogV,false,0,true);*/
    
    
    //-------------------------------------------------------------------------------------
    
    
    
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
    
 // ENDof *New!!   
    
 // ********************************************************************************************   


 
} //*NEW

/* G4LogicalVolume *G4d2oDetector::GetTeflon(){
    
    ///// Teflon reflector /////
    G4Box *teflonSolid = new G4Box("teflonSolid", h2oLength/2.0,teflonThickness/2.0, h2oHeight/2.0);
    G4LogicalVolume *teflonLogV = new G4LogicalVolume(teflonSolid,
                                                      matPtr->GetMaterial( TEFLON ),
                                                      "teflonLogV");
    
    return teflonLogV;
    
}

G4ThreeVector G4d2oDetector::GetPMTPosition(Int_t iPMT){
    
    if(!hPanel[0] || !hPanel[1] || !hPanel[2]){
        printf("Error in G4d2oDetector::GetPMTPosition() - pixel centers not defined.\n\n");
        exit(0);
    }
    
    G4ThreeVector thePosition(hPanel[0]->GetBinContent(iPMT+1),
                              hPanel[1]->GetBinContent(iPMT+1),
                              hPanel[2]->GetBinContent(iPMT+1));
    
    return thePosition;
  
} */


/*G4LogicalVolume *G4d2oDetector::GetEllipsoidPMT()
{
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
  G4UnionSolid *OPMTSPartSolid  = new G4UnionSolid(OPMTSPartName, OPMTTopESolid, OPMTBottomESolid, Transform);

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
  G4UnionSolid      *OPMTSolid  = new G4UnionSolid(OPMTVolName, OPMTSPartSolid, OPMTLegSolid, Transform);
  pmtLogV    = new G4LogicalVolume(OPMTSolid, matPtr->GetMaterial( BOROSILICATE ), OPMTVolName);
  

  // Inner PMT shape

  // FIrst the top ellipsoid
  G4String IPMTTopEName         = "IPMTtopEVol";
  G4double IPMTTopESemiX        = OPMTTopESemiX - PMTFrameThickness;
  G4double IPMTTopESemiY        = OPMTTopESemiY - PMTFrameThickness;
  G4double IPMTTopESemiZ        = OPMTTopESemiZ - PMTFrameThickness;
  G4double IPMTTopEBottomCut    = 0;
  G4double IPMTTopECut          = IPMTTopESemiZ;
  G4Ellipsoid *IPMTTopESolid = new G4Ellipsoid(IPMTTopEName,IPMTTopESemiX, IPMTTopESemiY, IPMTTopESemiZ, IPMTTopEBottomCut, IPMTTopECut);

  // THen Inner PMT Bottom Ellipsoid
  G4String IPMTBottomEName      = "IPMTbottomEVol";
  G4double IPMTBottomESemiX     = IPMTTopESemiX;
  G4double IPMTBottomESemiY     = IPMTTopESemiY;
  G4double IPMTBottomESemiZ     = IPMTTopESemiZ;
  G4double IPMTBottomEBottomCut = -IPMTBottomESemiZ;
  G4double IPMTBottomETopCut    = 0;
  G4Ellipsoid *IPMTBottomESolid = new G4Ellipsoid(IPMTBottomEName,IPMTBottomESemiX, IPMTBottomESemiY, IPMTBottomESemiZ, IPMTBottomEBottomCut, IPMTBottomETopCut);

  // Inner PMT Spherical part solid
  G4String  IPMTSPartName = "InnerPMTSphericalVol";
  Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0,0,0));
  G4UnionSolid *IPMTSPartSolid = new G4UnionSolid(IPMTSPartName, IPMTTopESolid,IPMTBottomESolid, Transform);

  // Inner PMT Leg
  G4String IPMTLegName      = "IPMTLegVol";
  G4double IPMTLegInnerRad  = 0;
  G4double IPMTLegOuterRad  = OPMTLegOuterRad - PMTFrameThickness;
  G4double IPMTLegZsize     = OPMTLegZsize - 2*PMTFrameThickness;
  G4double IPMTLegSphi      = 0;
  G4double IPMTLegPphi      = 2*M_PI*rad;
  G4double IPMTZposition    = 60.84*mm  + IPMTLegZsize/2;
  G4Tubs *IPMTLegSolid = new G4Tubs(IPMTLegName, IPMTLegInnerRad, IPMTLegOuterRad, IPMTLegZsize/2, IPMTLegSphi, IPMTLegPphi);

  // Inner PMT Full Shape
  G4String IPMTVolName = "InnerPMTVol";
  Transform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, -IPMTZposition));
  G4LogicalVolume   *IPMTLog    = new G4LogicalVolume(IPMTSPartSolid, matPtr->GetMaterial( VACUUM ), IPMTVolName);
  IPMTLog->SetVisAttributes(G4VisAttributes::Invisible);

  // Then place inner PMT volume inside outer PMT volume
  G4RotationMatrix *IPMTRotation = new G4RotationMatrix();
  IPMTRotation->rotateX(M_PI*rad);
  G4VPhysicalVolume *IPMTPhys      = new G4PVPlacement(0, G4ThreeVector(0,0,0), IPMTLog, IPMTVolName, pmtLogV, false, 0, true); //top

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
  G4VPhysicalVolume *PMTPhtPhys      = new G4PVPlacement(0, G4ThreeVector(0, 0, PMTPhtZpozition), PMTPhtLog, PMTPhtTopVolName, IPMTLog, false, 0,  true);  //top
  

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
}*/

//*


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


