//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Paraboloid.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fPWorld(0), fLWorld(0), fMaterial(0), 
 fPhysiPool(0), fLogicPool(0), fPoolMater(0), 
 fPhysiFilter3(0), fLogicFilter3(0), fFilter3Mater(0), 
 fPhysiReflector(0), fLogicReflector(0), fReflectorMater(0), 
 fPhysiNGuide(0), fLogicNGuide(0), fNGuideMater(0), 
 fPhysiDDGenerator(0), fLogicDDGenerator(0), fDDGeneratorMater(0),
 fPhysiNShield(0), fLogicNShield(0), fNShieldMater(0), 
 fDetectorMessenger(0)
{
  fWorldSize = 80*m;
  DefineMaterials();
  SetWorldMaterial("G4_AIR");
  
//  // ProtoDUNE Lar pool  
//  fPoolLength      = 8.9*m;
//  fPoolWidth      = 7.8*m;
//  fPoolHeight      = 7.3*m;
//  // ProtoDUNE Gas argon buffer
//  fBufferLength      = 8.9*m; 
//  fBufferWidth      = 7.8*m; 
//  fBufferHeight      = 0.8*m; 


// DUNE Lar pool  
  fPoolLength      = 58.0*m;
  fPoolWidth      = 14.5*m;
  fPoolHeight      = 12.0*m;  
  // DUNE Gas argon buffer
  fBufferLength      = 58.0*m;
  fBufferWidth      = 14.5*m;
  fBufferHeight      = 0.8*m;
  
  fCryostatThickness = 1.0*mm;
  fCryostatLength = fPoolLength + 2*fCryostatThickness;
  fCryostatWidth = fPoolWidth + 2*fCryostatThickness;
  fCryostatHeight= fPoolHeight + fBufferHeight + 2*fCryostatThickness;
  
  // Polyethylene insulator
  fInsulatorThickness = 90*cm;
  fInsulatorLength = fCryostatLength + 2*fInsulatorThickness;
  fInsulatorWidth = fCryostatWidth + 2*fInsulatorThickness;
  fInsulatorHeight= fCryostatHeight + 2*fInsulatorThickness;
  
  // neutron DD generator
  fDDGeneratorHeight = 50.0*cm;
  fDDGeneratorRadius = 5.0*cm;
  
  // Moderator
  fModerator_Thickness = 20*cm;
  fModerator_Height = fDDGeneratorHeight + fModerator_Thickness; // must be larger than the DD generator height
  fModerator_Radius = fDDGeneratorRadius + 8*cm;
  
  // 1st Filter
  fFilter1Height = 20*cm;
  fFilter1Radius = fModerator_Radius + 8.0*cm;
  
  // 2nd Filter
  fFilter2Height = 10*cm;
  fFilter2Radius_top = fFilter1Radius;
  fFilter2Radius_bottom = fFilter1Radius;
  fFilter2Radius_bottom = 12.5*cm;
  
  // 3rd Filter
  fFilter3Height        = 8.0*cm;
  fFilter3Radius_top    = fFilter2Radius_bottom;
  fFilter3Radius_bottom = fFilter2Radius_bottom;
  
  // neutron reflector
  fReflectorThickness = 12.0*cm;
  fReflectorHeight = fFilter1Height + fModerator_Height + fReflectorThickness;
  fReflectorRadius = fFilter1Radius + fReflectorThickness;
  
  // neutron guide
  fNGuideHeight = fFilter2Height;
  fNGuideRadius_top = fFilter2Radius_top + fReflectorThickness;
  fNGuideRadius_bottom = fFilter2Radius_bottom + fReflectorThickness;
  
  // neutron shield
  fNShieldThickness = 12.0*cm;
  fNShieldHeight = fReflectorHeight + fFilter2Height + fFilter3Height + fNShieldThickness;
  fNShieldRadius = fReflectorRadius + fNShieldThickness;
  
  // Feedthrough port
  //fPortHeight = fInsulatorThickness; // Moderator is on top of the cryostat port
  fPortHeight = 5*cm; // Moderator is inside the cryostat port
  fPortOuterRadius = fNShieldRadius;
  
  // Thermal neutron absorber
  fAbsorberHeight = 3*cm;
  fThermalAbsorberRadius = fPortOuterRadius;
  
  // define neutron source location
  
  // corner manhole location
  //fCenterX = -28.0* m;
  //fCenterY = -6.25* m;
  
  // center of cryostat roof
  fCenterX = 0* m;
  fCenterY = 0* m;
  
  // define colors
	blue_		= new G4VisAttributes(G4Colour(0, 0.4, 0.8));
	green_	= new G4VisAttributes(G4Colour(0, 0.6, 0.3));
	red_		= new G4VisAttributes(G4Colour(0.82, 0.1, 0.1));
	grey_		= new G4VisAttributes(G4Colour(0.545, 0.533, 0.470));
	silver_	= new G4VisAttributes(G4Colour(0.625, 0.625, 0.625));
	orange_	= new G4VisAttributes(G4Colour(1,     0.549, 0));
	copper_	= new G4VisAttributes(G4Colour(0.722, 0.451, 0.2));
	cyan_ 	= new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
	yellow_ = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
	pink_ = new G4VisAttributes(G4Colour(255./256,20./256,147./256));

  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::DefineMaterials()
{
  // specific element name for thermal neutronHP
  // (see G4ParticleHPThermalScatteringNames.cc)
  G4int Z, A, a, ncomponents, natoms;
  G4double fractionmass, abundance;
  
  // vacuum
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;  
  fVacuum = new G4Material("interGalactic", atomicNumber, massOfMole, density, kStateGas, temperature, pressure);
  
  // Define elements for all materials
  G4NistManager* man = G4NistManager::Instance();
  G4Element* H = man->FindOrBuildElement("H");
  G4Element* Li = man->FindOrBuildElement("Li");
  G4Element* B = man->FindOrBuildElement("B");
  G4Element* C  = man->FindOrBuildElement("C");  
  G4Element* O  = man->FindOrBuildElement("O");
  G4Element* F = man->FindOrBuildElement("F");
  G4Element* Na = man->FindOrBuildElement("Na");
  G4Element* Mg = man->FindOrBuildElement("Mg");
  G4Element* Al = man->FindOrBuildElement("Al");
  G4Element* Si = man->FindOrBuildElement("Si");
  G4Element* S = man->FindOrBuildElement("S");
  G4Element* K = man->FindOrBuildElement("K");
  G4Element* Ca = man->FindOrBuildElement("Ca");
  G4Element* Ti = man->FindOrBuildElement("Ti");
  G4Element* Cr = man->FindOrBuildElement("Cr");
  G4Element* Mn = man->FindOrBuildElement("Mn");
  G4Element* Fe = man->FindOrBuildElement("Fe");
  G4Element* Ni = man->FindOrBuildElement("Ni");
  G4Element* Sb = man->FindOrBuildElement("Sb");
  G4Element* Xe = man->FindOrBuildElement("Xe");
  G4Element* Cs = man->FindOrBuildElement("Cs");
  G4Element* Bi = man->FindOrBuildElement("Bi");
  
  // stainless steel
  G4Material* StainlessSteel = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
 
  
  // MgF2
  G4Material* MgF2 = new G4Material("MgF2", 3.15*g/cm3, ncomponents=2, kStateSolid);
  MgF2->AddElement(Mg, natoms=1);
  MgF2->AddElement(F, natoms=2);
  
  // TiF3
  G4Material* TiF3 = new G4Material("TiF3", 3.4*g/cm3, ncomponents=2, kStateSolid);
  TiF3->AddElement(Ti, natoms=1);
  TiF3->AddElement(F, natoms=3);
      
  // Fe-56 isotope
  G4Isotope* iso_Fe = new G4Isotope("iso_Fe", Z=26, A=56, a=55.9349363*g/mole);
  G4Element* ele_Fe = new G4Element("ele_Fe", "Fe", ncomponents=1);
  ele_Fe->AddIsotope(iso_Fe,abundance=100.*perCent);
  G4Material* mat_Fe=new G4Material("mat_Fe",7.874*g/cm3, ncomponents = 1);
	mat_Fe->AddElement(ele_Fe, fractionmass = 1 );
	
  // Li-6 isotope
  G4Isotope* iso_Li = new G4Isotope("iso_Li", Z=3, A=6, a=6.015122795*g/mole);
  G4Element* ele_Li = new G4Element("ele_Li", "Li", ncomponents=1);
  ele_Li->AddIsotope(iso_Li,abundance=100.*perCent);
  G4Material* Li6=new G4Material("Li6",0.534*g/cm3, ncomponents = 1);
  Li6->AddElement(ele_Li, fractionmass = 1 );
    
  // B-10 isotope
  G4Isotope* iso_B = new G4Isotope("iso_B", Z=5, A=10, a=10.013*g/mole);
  G4Element* ele_B = new G4Element("ele_B", "B", ncomponents=1);
  ele_B->AddIsotope(iso_B,abundance=100.*perCent);
  G4Material* B10=new G4Material("B10",2.34*g/cm3, ncomponents = 1);
  B10->AddElement(ele_B, fractionmass = 1 );
    
  // SiO2
  G4Isotope* iso_Si30 = new G4Isotope("iso_Si30", Z=14, A=30, a=29.9738*g/mole);
  G4Element* ele_Si30 = new G4Element("ele_Si30", "Si30", ncomponents=1);
  ele_Si30->AddIsotope(iso_Si30,abundance=100.*perCent);
  G4Material* Si30=new G4Material("Si30",2.3296*g/cm3, ncomponents = 1);
  Si30->AddElement(ele_Si30, fractionmass = 1 );
  G4Material* SiO2 = new G4Material("SiO2", 2.65*g/cm3, ncomponents=2, kStateSolid);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O, natoms=2);
  
  // SiS2
  G4Material* SiS2 = new G4Material("SiS2", 1.85*g/cm3, ncomponents=2, kStateSolid);
  SiS2->AddElement(Si, natoms=1);
  SiS2->AddElement(S, natoms=2);
  
  // Fluental
  G4Material* AlF3 = new G4Material("AlF3", 2.88*g/cm3, ncomponents=2, kStateSolid); //ALF3
  AlF3->AddElement(Al, natoms=1);
  AlF3->AddElement(F, natoms=3);
  G4Material* LiF = new G4Material("LiF", 2.64*g/cm3, ncomponents=2, kStateSolid);   //LiF
  LiF->AddElement(Li, natoms=1);
  LiF->AddElement(F, natoms=1);
  G4Material* Fluental = new G4Material ("Fluental", density=2.831, ncomponents=3);
  Fluental->AddElement (Al, fractionmass = 30*perCent);
  Fluental->AddMaterial (AlF3, fractionmass = 69*perCent);
  Fluental->AddMaterial (LiF, fractionmass = 1*perCent);
      
  //Teflon (CF2)
  G4Material* Teflon = new G4Material("Teflon", 2.20*g/cm3, ncomponents=2, kStateSolid);
  Teflon->AddElement(C, natoms=1);
  Teflon->AddElement(F, natoms=2); 
  
  //Al2O3
  G4Material* Al2O3 = new G4Material("Al203", 3.987*g/cm3, ncomponents=2, kStateSolid);
  Al2O3->AddElement(Al, natoms=2);
  Al2O3->AddElement(O, natoms=3);
  //BiF3
  G4Material* BiF3 = new G4Material("BiF3", 5.32*g/cm3, ncomponents=2, kStateSolid);
  BiF3->AddElement(Bi, natoms=1);
  BiF3->AddElement(F, natoms=3);
  //BiF5
  G4Material* BiF5 = new G4Material("BiF5", 5.40*g/cm3, ncomponents=2, kStateSolid);
  BiF5->AddElement(Bi, natoms=1);
  BiF5->AddElement(F, natoms=5);
  //CaF2
  G4Material* CaF2 = new G4Material("CaF2", 3.18*g/cm3, ncomponents=2, kStateSolid);
  CaF2->AddElement(Ca, natoms=1);
  CaF2->AddElement(F, natoms=2);  
  
  //Lithium Polyethylene
  G4Material* polyethylene = man->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material* LiPoly = new G4Material("LiPoly", 1.06*g/cm3, ncomponents=2);
  LiPoly->AddMaterial (Li6, 7.54*perCent);
  LiPoly->AddMaterial (polyethylene, 92.46*perCent);
  
  // world mater
  fWorldMater = man->FindOrBuildMaterial("G4_AIR");
  
//  // insulator material: protodune foam
//  fInsulatorMater = new G4Material("ProtoduneFoam", 0.135*g/cm3, 4);
//  fInsulatorMater->AddElement(C, 17);
//  fInsulatorMater->AddElement(H, 16);
//  fInsulatorMater->AddElement(N, 2);
//  fInsulatorMater->AddElement(O, 4); 
  
  fInsulatorMater = LiPoly; // insulator material: Lithium Polyethylene
  //fInsulatorMater =  man->FindOrBuildMaterial("G4_POLYETHYLENE"); // insulator material: pure Polyethylene
  
  // Feedthrough port
  fPortMater = man->FindOrBuildMaterial("G4_AIR");
  // feedthrough port reflector
  fPortRefMater = man->FindOrBuildMaterial("G4_Pb");
  // cryostat
  fCryostatMater = fVacuum;
  // liquid argon pool
  fPoolMater = man->FindOrBuildMaterial("G4_lAr");
  // gas argon buffer
  fBufferMater = man->FindOrBuildMaterial("G4_Ar");
  // neutron reflectorMater
  fReflectorMater = man->FindOrBuildMaterial("G4_Pb");
  // neutron guide
  fNGuideMater = man->FindOrBuildMaterial("G4_Pb");
  // DD generator
  fDDGeneratorMater = fVacuum; 
  // Moderator 
  fModerator_Mater = man->FindOrBuildMaterial("G4_Fe");
  // Filter 1 
  fFilter1Mater = man->FindOrBuildMaterial("G4_S");
  // Filter 2
  fFilter2Mater = man->FindOrBuildMaterial("G4_S");
  // Filter 3
  fFilter3Mater = man->FindOrBuildMaterial("G4_S");
  // neutron thermal absorber
  fThermalAbsorberMater = Li6;
  // neutron shield
  fNShieldMater = LiPoly; 
 
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;
 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);
 return material;
}          

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  // world box
  G4Box*
  sWorld = new G4Box("Container",                                
                     fWorldSize/2,fWorldSize/2,fWorldSize/2);    
  fLWorld = new G4LogicalVolume(sWorld,                          
                                fWorldMater,                     
                                "World_l");                      
  fPWorld = new G4PVPlacement(0,                                 
                              G4ThreeVector(),                   
                              fLWorld,                           
                              "World_p",                         
                              0,                                 
                              false,                             
                              0);  
  fLWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  	                             
  // insulator
  G4Box*
  sInsulator1 = new G4Box("Insulator1",                         
                   fInsulatorLength/2,fInsulatorWidth/2, fInsulatorHeight/2);     
  G4Tubs*
  sInsulator2 = new G4Tubs("Insulator2",
                   0, fPortOuterRadius, fInsulatorThickness/2, 0.,CLHEP::twopi );
  G4ThreeVector zTransInsulator(fCenterX, fCenterY, fInsulatorHeight/2 - fInsulatorThickness/2); 
  G4SubtractionSolid* sInsulator =
  new G4SubtractionSolid("Insulator1-Insulator2", sInsulator1, sInsulator2, 0, zTransInsulator); 
  fLogicInsulator = new G4LogicalVolume(sInsulator,
                             fInsulatorMater,
                             "Insulator_l"); 
  fPhysiInsulator = new G4PVPlacement(0,     
                            G4ThreeVector(0, 0, 0),
                            fLogicInsulator,       
                            "Insulator_p",
                            fLWorld,               
                            false,                 
                            0);  
 fLogicInsulator->SetVisAttributes(silver_); 
                  
  // Feedthrough port
  G4Tubs* 
  sPort = new G4Tubs("Port", 0, fPortOuterRadius, fPortHeight/2, 0., CLHEP::twopi );                   
  fLogicPort = new G4LogicalVolume(sPort, fPortMater, "Port_l"); 
  fPhysiPort = new G4PVPlacement(0, G4ThreeVector(fCenterX, fCenterY, fInsulatorHeight/2 - fInsulatorThickness + fPortHeight/2), 
                                 fLogicPort, "Port_p", fLWorld, false, 0);   
  fLogicPort->SetVisAttributes(red_);                                                         

  // cryostat
  G4Box*
  sCryostat = new G4Box("Cryostat",                                              
                        fCryostatLength/2, fCryostatWidth/2, fCryostatHeight/2); 
  fLogicCryostat = new G4LogicalVolume(sCryostat,                                
                                       fCryostatMater,                           
                                       "Cryostat_l");                            
  fPhysiCryostat = new G4PVPlacement(0,                                          
                                     G4ThreeVector(0, 0, 0),                     
                                     fLogicCryostat,                             
                                     "Cryostat_p",                               
                                     fLogicInsulator,                            
                                     false,                                      
                                     0);  
  fLogicCryostat-> SetVisAttributes(yellow_);  
                                      
  // liquid argon pool
  G4Box*
  sLarPool = new G4Box("Larpool",                                
                       fPoolLength/2,fPoolWidth/2,fPoolHeight/2);
  fLogicPool = new G4LogicalVolume(sLarPool,                     
                                   fPoolMater,                   
                                   "LarPool_l");                 
  fPhysiPool = new G4PVPlacement(0,                                                                       
                                 G4ThreeVector(0, 0, -fCryostatHeight/2+fCryostatThickness+fPoolHeight/2),
                                 fLogicPool,                                                              
                                 "LarPool_p",                                                             
                                 fLogicCryostat,                                                          
                                 false,                                                                   
                                 0);  
  fLogicPool-> SetVisAttributes(blue_);  
                                                                                                     
  // gas argon buffer
  G4Box*
  sGArBuffer = new G4Box("GArbuffer",                                    
                         fBufferLength/2,fBufferWidth/2,fBufferHeight/2);
  fLogicBuffer = new G4LogicalVolume(sGArBuffer,                         
                                     fBufferMater,                       
                                     "GArbuffer_l");                     
  fPhysiBuffer = new G4PVPlacement(0,                                                                                       
                                   G4ThreeVector(0, 0, -fCryostatHeight/2+fCryostatThickness+fPoolHeight + fBufferHeight/2),
                                   fLogicBuffer,                                                                            
                                   "GArbuffer_p",                                                                           
                                   fLogicCryostat,                                                                          
                                   false,                                                                                   
                                   0);  
  fLogicBuffer->SetVisAttributes(blue_); 
                                                                                           
  // neutron shield
  G4Tubs* 
  sNShield = new G4Tubs("NShield_s",                                            
                        0, fNShieldRadius, 0.5*fNShieldHeight, 0.,CLHEP::twopi);
  fLogicNShield = new G4LogicalVolume(sNShield,                                 
                                      fNShieldMater,                            
                                      "NShield_l");                             
  fPhysiNShield = new G4PVPlacement(0,                                                              
                                    G4ThreeVector(fCenterX, fCenterY,  fInsulatorHeight/2-fInsulatorThickness+fPortHeight  + fNShieldHeight/2),      
                                    fLogicNShield , 
                                    "NShield_p",    
                                    fLWorld,        
                                    false,          
                                    0);  
  fLogicNShield->SetVisAttributes(silver_);  
           
  // neutron thermal absorber                        
  G4Tubs* 
  sThermalAbsorber = new G4Tubs("ThermalAbsorber_s",    
                                0,                      
                                fThermalAbsorberRadius, 
                                fAbsorberHeight/2,      
                                0.,                     
                                CLHEP::twopi);          
  fLogicThermalAbsorber = new G4LogicalVolume(sThermalAbsorber,     
                                              fThermalAbsorberMater,
                                              "ThermalAbsorber_l"); 
  fPhysiThermalAbsorber = new G4PVPlacement(0, 
                                            G4ThreeVector(fCenterX, fCenterY,  -fPortHeight/2 + fAbsorberHeight/2),  
                                            fLogicThermalAbsorber, 
                                            "ThermalAbsorber_p",
                                            fLogicPort,
                                            false,
                                            0);
  fLogicThermalAbsorber->SetVisAttributes(copper_);
  
  // neutron reflector
  G4Tubs* 
  sReflector = new G4Tubs("Reflector_s",                 
                          0,                             
                          fReflectorRadius,              
                          fReflectorHeight/2,            
                          0.,                            
                          CLHEP::twopi);                 
  fLogicReflector = new G4LogicalVolume(sReflector,      
                                        fReflectorMater, 
                                        "Reflector_l");  
  fPhysiReflector = new G4PVPlacement(0,                 
                                      G4ThreeVector(0, 0, -fNShieldHeight/2 + fFilter3Height + fFilter2Height + fReflectorHeight/2),
                                      fLogicReflector ,  
                                      "Reflector_p",     
                                      fLogicNShield,     
                                      false,             
                                      0);                
  // Neutron guide
  G4Cons* 
  sNGuide = new G4Cons("NGuide_s",            
                        fFilter2Radius_bottom,
                        fNGuideRadius_bottom, 
                        fFilter2Radius_top,   
                        fNGuideRadius_top,    
                        fNGuideHeight/2,      
                        0.,                   
                        CLHEP::twopi);         
  fLogicNGuide = new G4LogicalVolume(sNGuide,      
                                      fNGuideMater,
                                      "NGuide_l"); 
  fPhysiNGuide = new G4PVPlacement(0,                                                                          
                                    G4ThreeVector(0, 0, -fNShieldHeight/2 + fFilter3Height + fNGuideHeight/2), 
                                    fLogicNGuide,                                                              
                                    "NGuide_p",                                                                
                                    fLogicNShield,                                                             
                                    false,                                                                     
                                    0);                                                                        
  
  // Filter 1                        
  G4Tubs* 
  sFilter1 = new G4Tubs("Filer1_s",                    
                        0,                             
                        fFilter1Radius,                
                        fFilter1Height/2,              
                        0.,                              
                        CLHEP::twopi);                 
  fLogicFilter1 = new G4LogicalVolume(sFilter1,        
                                      fFilter1Mater,   
                                      "Filter1_l");    
  fPhysiFilter1 = new G4PVPlacement(0,                                                           
                                    G4ThreeVector(0, 0, -fReflectorHeight/2 + fFilter1Height/2), 
                                    fLogicFilter1,                                               
                                    "Filter1_p",                                                 
                                    fLogicReflector,                                             
                                    false,                                                       
                                    0);  
  fLogicFilter1-> SetVisAttributes(pink_);  
                                                                                               
  // Filter 2
  G4Cons* 
  sFilter2 = new G4Cons("Filter2_s",               
                        0,                         
                        fFilter2Radius_bottom,     
                        0,                         
                        fFilter2Radius_top,        
                        fFilter2Height/2,          
                        0.,                        
                        CLHEP::twopi);             
  fLogicFilter2 = new G4LogicalVolume(sFilter2,     
                                      fFilter2Mater,
                                      "Filter2_l"); 
  fPhysiFilter2 = new G4PVPlacement(0,                                                                          
                                    G4ThreeVector(0, 0, -fNShieldHeight/2 + fFilter3Height + fFilter2Height/2), 
                                    fLogicFilter2,                                                              
                                    "Filter2_p",                                                                
                                    fLogicNShield,                                                              
                                    false,                                                                      
                                    0);                                                                         
  fLogicFilter2-> SetVisAttributes(pink_);
  
  // Filter 3                         
  G4Cons* 
  sFilter3 = new G4Cons("Filter3_s",  
                        0, 
                        fFilter3Radius_bottom, 
                        0, 
                        fFilter3Radius_top,
                        fFilter3Height/2, 
                        0.,
                        CLHEP::twopi);
  fLogicFilter3 = new G4LogicalVolume(sFilter3,      
                                      fFilter3Mater, 
                                      "Filter3_l");  
  fPhysiFilter3 = new G4PVPlacement(0,                                                         
                                    G4ThreeVector(0, 0, -fNShieldHeight/2 + fFilter3Height/2), 
                                    fLogicFilter3,                                             
                                    "Filter3_p",                                               
                                    fLogicNShield,                                             
                                    false,                                                     
                                    0);                                                        
  fLogicFilter3-> SetVisAttributes(pink_);
  
  // Moderator                         
  G4Tubs* 
  sModerator = new G4Tubs("Moderator_s",  
                           0, fModerator_Radius, 0.5*fModerator_Height, 0.,CLHEP::twopi);
  fLogicModerator = new G4LogicalVolume(sModerator,      
                                        fModerator_Mater,
                                        "Moderator_l");  
  fPhysiModerator = new G4PVPlacement(0,                 
                                      G4ThreeVector(0, 0, -fReflectorHeight/2 + fModerator_Height/2 + fFilter1Height),   
                                      fLogicModerator, 
                                      "Moderator_p",   
                                      fLogicReflector, 
                                      false,           
                                      0);              
  fLogicModerator-> SetVisAttributes(orange_);
  
  // neutron DD generator                         
  G4Tubs* 
  sDDGenerator = new G4Tubs("DDGenerator_s",  
                            0, fDDGeneratorRadius, 0.5*fDDGeneratorHeight, 0.,CLHEP::twopi);
  fLogicDDGenerator = new G4LogicalVolume(sDDGenerator,      
                                          fDDGeneratorMater, 
                                          "DDGenerator_l");  
  fPhysiDDGenerator = new G4PVPlacement(0,                   
                                        G4ThreeVector(0, 0, fModerator_Height/2 - fDDGeneratorHeight/2),
                                        fLogicDDGenerator, 
                                        "DDGenerator_p",   
                                        fLogicModerator,   
                                        false,             
                                        0);                
  fLogicDDGenerator-> SetVisAttributes(green_);               
  PrintParameters();
  
  //always return the root volume
  //
  return fPWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The World is " << G4BestUnit(fWorldSize,"Length")
         << " of " << fMaterial->GetName() 
         << "\n \n" << fMaterial << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial;
  if(materialChoice == "Vacuum") pttoMaterial = fVacuum;	
  // search the material by its name
  else pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial) { 
    if(fMaterial != pttoMaterial) {
      fMaterial = pttoMaterial;
      if(fLWorld) { fLWorld->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetWorldMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//SET and GET functions for moderator layers
void DetectorConstruction::GetMaterialTable(){
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetModeratorMaterial(G4String materialChoice)
{
  G4Material* pttoMaterial;
  if(materialChoice == "Vacuum") pttoMaterial = fVacuum;	
  // search the material by its name
  else pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice); 
  if (pttoMaterial) { 
    if(fModerator_Mater != pttoMaterial) {
      fModerator_Mater = pttoMaterial;
      if(fLogicModerator) { fLogicModerator->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetModeratorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::GetModeratorMaterial(){
    G4cout << "Moderator material is " << fModerator_Mater->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetFilter1Material(G4String materialChoice)
{ 
	G4Material* pttoMaterial;
  if(materialChoice == "Vacuum") pttoMaterial = fVacuum;	
  // search the material by its name
  else pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial) { 
    if(fFilter1Mater != pttoMaterial) {
      fFilter1Mater = pttoMaterial;
      if(fLogicFilter1) { fLogicFilter1->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetFilter1Material : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::GetFilter1Material(){
    G4cout << "Filter 1 material is " << fFilter1Mater->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetFilter2Material(G4String materialChoice)
{ 
	G4Material* pttoMaterial;
  if(materialChoice == "Vacuum") pttoMaterial = fVacuum;
  // search the material by its name
  else pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial) { 
    if(fFilter2Mater != pttoMaterial) {
      fFilter2Mater = pttoMaterial;
      if(fLogicFilter2) { fLogicFilter2->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetFilter1Material : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::GetFilter2Material(){
    G4cout << "Filter 2 material is " << fFilter2Mater->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetFilter3Material(G4String materialChoice)
{  
	G4Material* pttoMaterial;
  if(materialChoice == "Vacuum") pttoMaterial = fVacuum;	
  // search the material by its name
  else pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial) { 
    if(fFilter3Mater != pttoMaterial) {
      fFilter3Mater = pttoMaterial;
      if(fFilter3Mater) { fLogicFilter3->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetFilter3Material : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::GetFilter3Material(){
    G4cout << "Filter 3 material is " << fFilter3Mater->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{  
	G4Material* pttoMaterial;
  if(materialChoice == "Vacuum") pttoMaterial = fVacuum;	
  // search the material by its name
  else pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial) { 
    if(fThermalAbsorberMater != pttoMaterial) {
      fThermalAbsorberMater = pttoMaterial;
      if(fThermalAbsorberMater) { fLogicThermalAbsorber->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetModeratorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::GetAbsorberMaterial(){
    G4cout << "Absorber material is " << fThermalAbsorberMater->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetPortRefMaterial(G4String materialChoice)
{  
	G4Material* pttoMaterial;
  if(materialChoice == "Vacuum") pttoMaterial = fVacuum;	
  // search the material by its name
  else pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial) { 
    if(fPortRefMater != pttoMaterial) {
      fPortRefMater = pttoMaterial;
      if(fPortRefMater) { fLogicPortRef->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetModeratorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::GetPortRefMaterial(){
    G4cout << "Port reflector material is " << fPortRefMater->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetReflectorMaterial(G4String materialChoice)
{
	G4Material* pttoMaterial;
  if(materialChoice == "Vacuum") pttoMaterial = fVacuum;	
  // search the material by its name
  else pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial) { 
    if(fReflectorMater != pttoMaterial) {
      fReflectorMater = pttoMaterial;
      if(fReflectorMater) { fLogicReflector->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetModeratorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::GetReflectorMaterial(){
    G4cout << "Neutron reflector material is " << fReflectorMater << G4endl;
}
void DetectorConstruction::SetNGuideMaterial(G4String materialChoice)
{
	G4Material* pttoMaterial;
  if(materialChoice == "Vacuum") pttoMaterial = fVacuum;	
  // search the material by its name
  else pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  if (pttoMaterial) { 
    if(fNGuideMater != pttoMaterial) {
      fNGuideMater = pttoMaterial;
      if(fNGuideMater) { fLogicNGuide->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetModeratorMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::GetNGuideMaterial(){
    G4cout << "Neutron guide material is " << fNGuideMater << G4endl;
}

//END OF SET AND GET FUNCTIONS
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::SetSize(G4double value)
{
  fWorldSize = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Set Dimension Functions
void DetectorConstruction::SetModeratorThickness(G4double value)
{
  fModerator_Thickness = value;
  fModerator_Height = fDDGeneratorHeight + fModerator_Thickness;
  fReflectorHeight = fFilter1Height + fModerator_Height + fReflectorThickness;
  fNShieldHeight = fFilter3Height + fFilter2Height + fReflectorHeight + fNShieldThickness;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetFilter1Height(G4double value)
{
  fFilter1Height = value;
  fReflectorHeight = fFilter1Height + fModerator_Height + fReflectorThickness;
  fNShieldHeight = fFilter3Height + fFilter2Height + fReflectorHeight + fNShieldThickness;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetFilter2Height(G4double value)
{
  fFilter2Height = value;
  fNShieldHeight = fFilter3Height + fFilter2Height + fReflectorHeight + fNShieldThickness;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetFilter3Height(G4double value)
{
  fFilter3Height = value;
  fNShieldHeight = fFilter3Height + fFilter2Height + fReflectorHeight + fNShieldThickness;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetAbsorberHeight(G4double value)
{
  fAbsorberHeight = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetPortRefThickness(G4double value)
{
  fPortRefThickness = value;
  fPortRefInnerRadius =  fPortOuterRadius - fPortRefThickness;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetReflectorThickness(G4double value)
{
  fReflectorThickness = value;
  fReflectorHeight = fFilter1Height + fModerator_Height + fReflectorThickness;
  fReflectorRadius = fModerator_Radius + fReflectorThickness;
  fNShieldHeight = fFilter3Height + fFilter2Height + fReflectorHeight + fNShieldThickness;
  fNShieldRadius = fReflectorRadius + fNShieldThickness;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetNShieldThickness(G4double value)
{
  fNShieldThickness = value;
  fNShieldHeight = fFilter3Height + fFilter2Height + fReflectorHeight + fNShieldThickness;
  fNShieldRadius = fReflectorRadius + fNShieldThickness;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetPortHeight (G4double value) 
{
  fPortHeight = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetSourceCenterX(G4double value)
{
  fCenterX = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
void DetectorConstruction::SetSourceCenterY(G4double value)
{
  fCenterY = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
