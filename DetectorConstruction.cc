
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "TrackerSD.hh"
#include "G4OpticalSurface.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4ElementTable.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include <G4SubtractionSolid.hh>
#include "G4PhysicalConstants.hh"
#include "G4PVReplica.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4VHit.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4bool checkOverlaps = true;
  // ------------- Materials -------------
  G4double a, z, density;
 
  G4int nelements;

  G4NistManager *nistManager = G4NistManager::Instance();
    G4bool isotopes = false;

    G4Material *fAir = nistManager->FindOrBuildMaterial("G4_AIR");
    
    
    G4Element*  O = nistManager->FindOrBuildElement("O" , isotopes);
    G4Element* Si = nistManager->FindOrBuildElement("Si", isotopes);
    G4Element* Lu = nistManager->FindOrBuildElement("Lu", isotopes);
    G4Element* Y = nistManager->FindOrBuildElement("Y", isotopes);
    G4Element* N = nistManager->FindOrBuildElement("N", isotopes);
    G4Element* C = nistManager->FindOrBuildElement("C", isotopes);
    G4Element* H = nistManager->FindOrBuildElement("H", isotopes);
    G4Element* F = nistManager->FindOrBuildElement("F", isotopes);
    
    
  

    G4Material* LYSO = new G4Material("LYSO", 7.4*g/cm3, 4);
    LYSO->AddElement(Lu, 6);
    LYSO->AddElement(Si, 10);
    LYSO->AddElement(O , 50);
    LYSO->AddElement(Y, 14);
    
    G4Material* BC422Q = new G4Material("BC422Q", 1.032*g/cm3, 2);
    BC422Q->AddElement(C, 10);
    BC422Q->AddElement(H, 11);
    
    G4Material *fTeflon = new G4Material("Teflon", 2.2*g/cm3, 2);
  fTeflon->AddElement(C, 0.240183);
  fTeflon->AddElement(F, 0.759817);

  G4Material* Al = new G4Material("Aluminium", 13., 26.98*g/mole, 2.700*g/cm3);
  G4Material* Ni = new G4Material("Nickel", 28., 58.69*g/mole, 8.910*g/cm3);
  G4Material* Au = new G4Material("Gold", 79., 196.97*g/mole, 19.32*g/cm3);
  G4Material* Ag = new G4Material("Silber", 47., 107.87*g/mole, 10.49*g/cm3);
  G4Material* Pb = new G4Material("Blei", 82., 207.20*g/mole, 11.35*g/cm3);

  std::vector<G4double> pdTeflonPhotonMomentum  = { 1.0*eV, 6.05*eV};
  std::vector<G4double> pdTeflonRefractiveIndex = {1.3, 1.3};
  std::vector<G4double> pdTeflonReflectivity = {1., 1.};
  //std::vector<G4double> pdTeflonBackscatter = {0.8, 0.8};
  std::vector<G4double> pdTeflonAbs = {0.001 *um, 0.001 *um};
  
  G4MaterialPropertiesTable *pTeflonPropertiesTable = new G4MaterialPropertiesTable();

  // Optical Properties
  pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, pdTeflonRefractiveIndex);
  pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, pdTeflonReflectivity);
  //pTeflonPropertiesTable->AddProperty("BACKSCATTERCONSTANT", pdTeflonPhotonMomentum, pdTeflonBackscatter);
  pTeflonPropertiesTable->AddProperty("ABSLENGTH", pdTeflonPhotonMomentum, pdTeflonAbs);
  fTeflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);
  


	
  std::vector<G4double> energy = {3.54*eV, 3.44*eV,3.35*eV, 3.29*eV, 3.26*eV, 3.18*eV, 3.10*eV, 3.02*eV, 2.95*eV, 2.88*eV, 2.82*eV, 2.76*eV, 2.70*eV};
	
  std::vector<G4double> rAir = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  
	std::vector<G4double> fast = {0.07, 0.43, 0.78, 1.00, 0.97, 0.80, 0.65, 0.51, 0.38, 0.29, 0.17, 0.09, 0.01};

	std::vector<G4double> rBC422Q = {1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58};
  
	
	std::vector<G4double> abs = {8.0*cm,8.0*cm,8.0*cm,8.0*cm,8.0*cm,8.0*cm,8.0*cm,8.0*cm,8.0*cm,8.0*cm,8.0*cm,8.0*cm,8.0*cm};


	std::vector<G4double> fraction = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    std::vector<G4double> reflectivity = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    
    
    G4MaterialPropertiesTable *mptBC422Q = new G4MaterialPropertiesTable();
    
    mptBC422Q->AddProperty("SCINTILLATIONCOMPONENT1", energy, fast);
    //mptBC422Q->AddProperty("SCINTILLATIONCOMPONENT2", energy, fast);
    mptBC422Q->AddProperty("RINDEX",energy, rBC422Q);
    mptBC422Q->AddProperty("ABSLENGTH", energy, abs);
    mptBC422Q->AddConstProperty("SCINTILLATIONYIELD", 3306./MeV,true);
    mptBC422Q->AddConstProperty("RESOLUTIONSCALE", 1.0);
    mptBC422Q->AddConstProperty("SCINTILLATIONRISETIME1", 110.*ps,true);
    mptBC422Q->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 0.7*ns,true);
    mptBC422Q->AddConstProperty("SCINTILLATIONYIELD1", 1.);
    BC422Q->SetMaterialPropertiesTable(mptBC422Q);
  
 
/// for LYSO as Scintillator
  /*std::vector<G4double> energy   =  {1.79*eV, 1.85*eV, 1.91*eV, 1.97*eV, 2.04*eV, 2.11*eV, 2.19*eV, 2.27*eV, 2.36*eV, 2.45*eV, 2.56*eV, 2.67*eV, 2.80*eV, 2.94*eV, 3.09*eV, 3.25*eV,3.44*eV, 3.65*eV, 3.89*eV, 4.16*eV};



  std::vector<G4double> fastLYSO  =  {0.01, 0.10, 0.20, 0.50, 0.90, 1.70, 2.90, 5.00, 8.30, 12.5, 17.0, 22.9,26.4, 25.6, 16.8, 4.20, 0.30, 0.20, 0.10, 0.01};



  std::vector<G4double> rLyso =  {1.81, 1.81, 1.81, 1.81, 1.81, 1.81, 1.81, 1.81, 1.81, 1.81, 1.81, 1.81,1.81, 1.81, 1.81, 1.81, 1.81, 1.81, 1.81, 1.81};



  std::vector<G4double> absLYSO  =  {1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm,1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm, 1.2*cm};
  
  std::vector<G4double> rAir = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    

    G4MaterialPropertiesTable *mptLYSO = new G4MaterialPropertiesTable();

    mptLYSO->AddProperty("SCINTILLATIONCOMPONENT1",  energy, fastLYSO);
    //mptBC422Q->AddProperty("SCINTILLATIONCOMPONENT2", energy, fast);
    mptLYSO->AddProperty("RINDEX",energy, rLyso);
    mptLYSO->AddProperty("ABSLENGTH", energy, absLYSO);
    mptLYSO->AddConstProperty("SCINTILLATIONYIELD", 33200./MeV,true);
    mptLYSO->AddConstProperty("RESOLUTIONSCALE", 1.0);
    mptLYSO->AddConstProperty("SCINTILLATIONRISETIME1", 70.*ps,true);
    mptLYSO->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 36.*ns,true);
    mptLYSO->AddConstProperty("SCINTILLATIONYIELD1", 1.);
    LYSO->SetMaterialPropertiesTable(mptLYSO);*/

 

 
    G4MaterialPropertiesTable *mptAir = new G4MaterialPropertiesTable();
    mptAir->AddProperty("RINDEX", energy, rAir, true);

    
    
    
    fAir->SetMaterialPropertiesTable(mptAir);
	
  G4MaterialPropertiesTable *OpSurfaceProperty = new G4MaterialPropertiesTable();


  G4OpticalSurface* opTeflonSurface = new G4OpticalSurface("TeflonSurface");
  opTeflonSurface->SetType(dielectric_LUT);
  opTeflonSurface->SetModel(LUT);
  opTeflonSurface->SetFinish(polishedteflonair);
  
  
  opTeflonSurface->SetMaterialPropertiesTable(OpSurfaceProperty);
	
	
	
    //World size
    G4double worldSizeXY = 20. * cm;
    G4double worldSizeZ = 20. * cm;

    //Crystal Parameters
    G4double cryst_dX = 40.* mm, cryst_dY = 40.* mm, cryst_dZ = 10.* mm;
    G4double foilThickness = 1 * mm;
    


    // Get materials
    G4Material* defaultMaterial = G4Material::GetMaterial("G4_AIR");
    
    G4Material* cryst_mat = G4Material::GetMaterial("BC422Q");
    
      //G4Material* cryst_mat = G4Material::GetMaterial("LYSO");

    // World
    //
    G4VSolid *worldS
            = new G4Box("World",           // its name
                        worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2); // its size

    G4LogicalVolume *worldLV
            = new G4LogicalVolume(
                    worldS,           // its solid
                    defaultMaterial,  // its material
                    "World");         // its name

    G4VPhysicalVolume *worldPV
            = new G4PVPlacement(
                    0,                // no rotation
                    G4ThreeVector(),  // at (0,0,0)
                    worldLV,          // its logical volume
                    "World",          // its name
                    0,                // its mother  volume
                    false,            // no boolean operation
                    0,                // copy number
                    fCheckOverlaps);  // checking overlaps
                    
	
	
	G4MaterialPropertiesTable *mptWorld = new G4MaterialPropertiesTable();
	mptWorld->AddProperty("RINDEX", energy, rAir, 2);
	defaultMaterial->SetMaterialPropertiesTable(mptWorld);

    

   
    G4double cryst_xpos = cryst_dX , cryst_ypos = cryst_dY ;

    //Pyrmiden geometry for crystal>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    G4Trd* crystalS = new G4Trd("crystal", cryst_dX/2, 10.* mm, cryst_dY/2, 10.* mm, cryst_dZ/2);
  
    //Box geometry for crystal>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    //G4Box* crystalS = new G4Box("crystal", cryst_dX/2, cryst_dY/2, cryst_dZ/2);

    //Geometry for Cylinder Crystals>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    //G4Tubs *crystalS = new G4Tubs("crystal", 0.,cryst_dX/2, cryst_dZ/2, 0.,2*M_PI );

    //Geometry for Cons>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    //G4Cons *crystalS = new G4Cons("crystal", 0.,5.*mm, 0., cryst_dX/2, cryst_dZ/2, 0.,2*M_PI );

    G4LogicalVolume* crystalLV
            = new G4LogicalVolume(
                    crystalS,     // its solid
                    cryst_mat,  // its material0.75*mm
                    "Crystal");   // its name
    
    //Sample Cration and if you want shielding
    G4Tubs* Sample = new G4Tubs("sample", 0.,5.*mm, 0.75*mm, 0.,2*M_PI );
    G4Tubs* SampleShield = new G4Tubs("sampleShield", 0.,cryst_dX/2, 0.75*mm, 0.,2*M_PI );

    G4VSolid * substractSampleshiled = new G4SubtractionSolid("substractSampleshiled", SampleShield, Sample,0, G4ThreeVector(0.,0.,0.0*mm));

    G4LogicalVolume* sampleLVShiedl
            = new G4LogicalVolume(
                    SampleShield,     // its solid
                    Pb,  // its material
                    "sampleShield");


    G4LogicalVolume* sampleLV
            = new G4LogicalVolume(
                    Sample,     // its solid
                    Al,  // its material
                    "sample");

    

     G4LogicalSkinSurface *skin = new G4LogicalSkinSurface("skin",worldLV, mirrorSurface); 
    

    //Box Geometrie Detektor>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    G4Box * solidDetector = new G4Box("solidDetector", cryst_dX/2, cryst_dY/2, 0.5*mm);
    logicDetector = new G4LogicalVolume(solidDetector, defaultMaterial, "logicDetector");


    //Detector End-wrapping Cone
    /*G4Tubs *DetectorWrapping = new G4Tubs("DetectorWrapping",0.,5.*mm, 0.5*mm, 0.,2*M_PI);

    G4VSolid *DetectorSubstract = new G4SubtractionSolid("DetectorSubstract", DetectorWrapping, solidDetector,0, G4ThreeVector(0.,0.,0.0*mm));

    logicDetectorSubstract = new G4LogicalVolume(DetectorSubstract, fTeflon, "logicDetectorSubstract");*/

    //Detector End-wrapping Box
    /*G4Box *DetectorWrapping = new G4Box("DetectorWrapping",30. *mm,30.*mm, 0.5*mm);

    G4VSolid *DetectorSubstract = new G4SubtractionSolid("DetectorSubstract", DetectorWrapping, solidDetector,0, G4ThreeVector(0.,0.,0.0*mm));

    logicDetectorSubstract = new G4LogicalVolume(DetectorSubstract, fTeflon, "logicDetectorSubstract");*/
    
    //Cylinder Geometrie Detektor>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    /*G4Tubs *solidDetector = new G4Tubs("solidDetector", 0.,cryst_dX/2, 0.5*mm, 0.,2*M_PI );

    logicDetector = new G4LogicalVolume(solidDetector, defaultMaterial, "logicDetector");*/
   
   
    //wrapping for cone geometry>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /*G4VSolid* crystalShape = new G4Cons("crystalshape",  0.,5.*mm, 0., cryst_dX/2, cryst_dZ/2, 0.,2*M_PI);

    G4Cons * BigTrd = new G4Cons("BigTrd",  0.,5.*mm + 0.5*mm, 0., cryst_dX/2+0.5*mm, cryst_dZ/2, 0.,2*M_PI);
    
    G4VSolid * substract = new G4SubtractionSolid("substract", BigTrd, crystalShape,0, G4ThreeVector(0.,0.,0.0*mm));

    logicSubstract = new G4LogicalVolume(substract , fTeflon, "logicSubstract");*/


    //wrapping for cylinder geometry>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /*G4VSolid* crystalShape = new G4Tubs("crystalshape", 0.,cryst_dX/2, cryst_dZ/2, 0.,2*M_PI);

    G4Tubs * BigTrd = new G4Tubs("BigTrd", 0.,cryst_dX/2+0.5*mm, cryst_dZ/2, 0.,2*M_PI);
    
    G4VSolid * substract = new G4SubtractionSolid("substract", BigTrd, crystalShape,0, G4ThreeVector(0.,0.,0.0*mm));

    logicSubstract = new G4LogicalVolume(substract , fTeflon, "logicSubstract");*/


    //wrapping for pyramid geometry>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    G4VSolid* crystalShape = new G4Trd("crystalshape", cryst_dX/2, 10.* mm, cryst_dY/2, 10.* mm, cryst_dZ/2);

    G4Trd * BigTrd = new G4Trd("BigTrd", cryst_dX/2 + 1*mm, 10.* mm+0.5*mm , cryst_dY/2 + 1*mm, 10.* mm +0.5*mm, cryst_dZ/2);
    
    G4VSolid * substract = new G4SubtractionSolid("substract", BigTrd, crystalShape,0, G4ThreeVector(0.,0.,0.0*mm));

    logicSubstract = new G4LogicalVolume(substract , fTeflon, "logicSubstract");
    

    //wrapping Box>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    /*G4VSolid* crystalShape = new G4Box("crystalshape", cryst_dX/2, cryst_dY/2, cryst_dZ/2);

    G4Box * BigBox = new G4Box("BigBox", cryst_dX/2 + 1*mm, cryst_dY/2 + 1*mm, cryst_dZ/2);
    
    G4VSolid * substract = new G4SubtractionSolid("substract", BigBox, crystalShape,0, G4ThreeVector(0.,0.,0.0*mm));

    logicSubstract = new G4LogicalVolume(substract , fTeflon, "logicSubstract");*/

    //Front wrapping Box>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    G4Box * solidFoilD = new G4Box("solidFoilD",  10.*mm+0.5*mm,  10.*mm+0.5*mm,   0.2 * mm);
    
    logicFoilD = new G4LogicalVolume(solidFoilD, fTeflon, "logicFoilD");


    //Front Wrapping cylinder>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    /*G4Tubs * solidFoilD = new G4Tubs("solidFoilD",  0.,cryst_dX/2+0.5*mm, 0.2*mm, 0.,2*M_PI );
    
    logicFoilD = new G4LogicalVolume(solidFoilD, fTeflon, "logicFoilD");*/


  
for(G4int i=0; i<2; i++)
	{
		//90 Grad Arrangment>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		/*
    G4Rotate3D rotZ((i*90)*deg, G4ThreeVector(1,0,0));
		G4Rotate3D rotZsample((315)*deg, G4ThreeVector(1,0,0));
    G4Translate3D transXScint(G4ThreeVector(0,i*25.*mm, 1.5*i*(cryst_dZ+(2*Z*mm))-(cryst_dZ/2+Z*mm)));
		G4Transform3D transformScint =(transXScint)*(rotZ);
		
		G4Translate3D transXDet(G4ThreeVector(0,i*30.5*mm, 1.125*i*(2*cryst_dZ+(2*Z*mm)+1.*mm)-(cryst_dZ+Z*mm+0.5*mm)));
		G4Transform3D transformDet =(transXDet)*(rotZ);
		
		G4Translate3D transXFolO(G4ThreeVector(0,i*25.*mm, 1.5*i*(cryst_dZ+(2*Z*mm))-(cryst_dZ/2+Z*mm)));
		G4Transform3D transformFolO =(transXFolO)*(rotZ);
		
    G4Translate3D transXFolD(G4ThreeVector(0,i*19.75*mm, 2.9*i*(2*Z*mm-0.4*mm)-(Z*mm-0.2*mm)));
		G4Transform3D transformFolD =(transXFolD)*(rotZ);

    G4Translate3D transXSample(G4ThreeVector(0,0, 20.*mm+i*(3*mm)-(0.5*mm)));
		G4Transform3D transformSample =(transXSample)*(rotZsample);*/

    G4Rotate3D rotZ((i*180)*deg, G4ThreeVector(1,1,0));
    //Abstand der Scintilatoren
    G4int Z = 4;

		G4Translate3D transXScint(G4ThreeVector(0,0, i*(cryst_dZ+(2*Z*mm))-(cryst_dZ/2+Z*mm)));
		G4Transform3D transformScint =(transXScint)*(rotZ);
		
		G4Translate3D transXDet(G4ThreeVector(0,0,i*(2*cryst_dZ+(2*Z*mm)+1.*mm)-(cryst_dZ+Z*mm+0.5*mm)));
		G4Transform3D transformDet =(transXDet)*(rotZ);
		
		G4Translate3D transXFolO(G4ThreeVector(0,0, i*(cryst_dZ+(2*Z*mm))-(cryst_dZ/2+Z*mm)));
		G4Transform3D transformFolO =(transXFolO)*(rotZ);
		
    G4Translate3D transXFolD(G4ThreeVector(0,0.,i*(2*Z*mm-0.4*mm)-(Z*mm-0.2*mm)));
		G4Transform3D transformFolD =(transXFolD)*(rotZ);

    /*G4Translate3D transXDetWrap(G4ThreeVector(0,0,i*(2*cryst_dZ+(2*Z*mm)+1.*mm)-(cryst_dZ+Z*mm+0.5*mm)));
		G4Transform3D transformDetWrap =(transXDetWrap)*(rotZ);*/

    G4Translate3D transXSample(G4ThreeVector(0,0.,i*(2.5*mm)-(1.25*mm)));
		G4Transform3D transformSample =(transXSample);
		
		
		G4VPhysicalVolume *physDetector = new G4PVPlacement(transformDet, logicDetector, "physDetector", worldLV, false,i, true);
    //G4VPhysicalVolume *physDetectorWrap = new G4PVPlacement(transformDetWrap, logicDetectorSubstract, "physDetectorWrap", worldLV, false,i, true);
    //G4VPhysicalVolume *physSample = new G4PVPlacement(transformSample, sampleLV, "physDSample", worldLV, false,i, true);
    //G4VPhysicalVolume *physSampleShield = new G4PVPlacement(transformSample, sampleLVShiedl, "physDSampleShield", worldLV, false,i, true);
		
		cAbsorberPV
            = new G4PVPlacement(
            transformScint,  // at (0,0,0)
            crystalLV,          // its logical volume
            "crystalPV",    // its name
            worldLV,          // its mother  volume
            false,            // no boolean operation
            i,                // copy number
            fCheckOverlaps);  // checking overlaps
            
            physFoilO = new G4PVPlacement(transformFolO, logicSubstract, "physFoilO", worldLV, false, i, true);
            
            physFoilD = new G4PVPlacement(transformFolD, logicFoilD, "physFoilD", worldLV, false, i, true);

		
  
  
   
    }

    G4LogicalSkinSurface* teflonSurface1 = new G4LogicalSkinSurface(
    "TeflonSurface1", logicSubstract, opTeflonSurface);
    G4LogicalSkinSurface* teflonSurface5 = new G4LogicalSkinSurface(
    "TeflonSurface5", logicFoilD, opTeflonSurface);
    

    fScoringVolume = crystalLV;
    
    
    
   
   

    

    G4VisAttributes *simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    simpleBoxVisAtt->SetVisibility(true);
    

    //
    // Always return the physical World
    //
    return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerSDname = "/TrackerChamberSD";
  TrackerSD* aTrackerSD = new TrackerSD(trackerSDname);
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  
  
  
 
  SetSensitiveDetector("logicDetector", aTrackerSD, true);
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}
void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((fStepLimit)&&(maxStep>0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


