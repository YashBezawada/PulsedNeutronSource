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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
    virtual G4VPhysicalVolume* Construct();

    G4Material* 
    MaterialWithSingleIsotope(G4String, G4String, G4double, G4int, G4int);
         
    void SetSize     (G4double);              
    void SetWorldMaterial (G4String);
    
    
    
    //Set and Get Material Functions
    void GetMaterialTable();
    void SetModeratorMaterial (G4String);
    void GetModeratorMaterial ();
    void SetFilter1Material (G4String);
    void GetFilter1Material ();
    void SetFilter2Material (G4String);
    void GetFilter2Material ();
    void SetFilter3Material (G4String);
    void GetFilter3Material ();
    void SetAbsorberMaterial (G4String);
    void GetAbsorberMaterial ();
    void SetPortRefMaterial (G4String);
    void GetPortRefMaterial ();
    void SetReflectorMaterial (G4String);
    void GetReflectorMaterial ();
    void SetNGuideMaterial (G4String);
    void GetNGuideMaterial ();

    //Set Dimension Functions
    void SetModeratorThickness(G4double);
    void SetFilter1Height(G4double);
    void SetFilter2Height(G4double);
    void SetFilter3Height(G4double);
    void SetAbsorberHeight(G4double);
    void SetPortRefThickness(G4double);
    void SetReflectorThickness(G4double);
    void SetNGuideThickness(G4double);
    void SetNShieldThickness(G4double);
    void SetPortHeight (G4double value);
    
    void SetSourceCenterX(G4double value);
    void SetSourceCenterY(G4double value);
    

  public:
  
     const
     G4VPhysicalVolume* GetWorld()              {return fPWorld;};   
     G4VPhysicalVolume* GetPhysiPool()          {return fPhysiPool;};  
     G4VPhysicalVolume* GetPhysiBuffer()        {return fPhysiBuffer;}; 
     G4VPhysicalVolume* GetPhysiCryostat()        {return fPhysiCryostat;}; 
     G4VPhysicalVolume* GetPhysiPort()          {return fPhysiPort;}; 
     G4VPhysicalVolume* GetPhysif()             {return fPhysiNShield;};  
     G4VPhysicalVolume* GetPhysiDDGenerator()   {return fPhysiDDGenerator;};
     G4VPhysicalVolume* GetPhysiReflector()     {return fPhysiReflector;}; 
     G4VPhysicalVolume* GetPhysiNGuide()        {return fPhysiNGuide;}; 
     G4VPhysicalVolume* GetPhysiModerator()     {return fPhysiModerator;};   
     G4VPhysicalVolume* GetPhysiFilter1()       {return fPhysiFilter1;};
     G4VPhysicalVolume* GetPhysiFilter2()       {return fPhysiFilter2;};  
     G4VPhysicalVolume* GetPhysiFilter3()       {return fPhysiFilter3;};  
      
     
     G4LogicalVolume* GetLogicWorld()           {return fLWorld;};     
     G4LogicalVolume* GetLogicPool()            {return fLogicPool;}; 
     G4LogicalVolume* GetLogicBuffer()          {return fLogicBuffer;}; 
     G4LogicalVolume* GetLogicCryostat()          {return fLogicCryostat;};
     G4LogicalVolume* GetLogicPort()            {return fLogicPort;}; 
     G4LogicalVolume* GetLogicShield()          {return fLogicNShield;}; 
     G4LogicalVolume* GetLogicReflector()       {return fLogicReflector;};  
     G4LogicalVolume* GetLogicNGuide()          {return fLogicNGuide;}; 
     G4LogicalVolume* GetLogicDDGenerator()     {return fLogicDDGenerator;}; 
     G4LogicalVolume* GetLogicModerator()       {return fLogicModerator;};   
     G4LogicalVolume* GetLogicFilter1()         {return fLogicFilter1;};
     G4LogicalVolume* GetLogicFilter2()         {return fLogicFilter2;};
     G4LogicalVolume* GetLogicFilter3()         {return fLogicFilter3;}; 
     G4LogicalVolume* GetLogicThermalAbsorber() {return fLogicThermalAbsorber;}; 
     G4LogicalVolume* GetLogicInsulator() {return fLogicInsulator;}; 
                   
                    
     G4double           GetSize()                   {return fWorldSize;};      
     G4Material*        GetMaterial()               {return fMaterial;}; 
     G4double           GetPoolHeight()             {return fPoolHeight;};
     G4double           GetModeratorHeight()        {return fModerator_Height;};
     G4double           GetFilter1Height()          {return fFilter1Height;};
     G4double           GetFilter2Height()          {return fFilter2Height;};
     G4double           GetFilter3Height()          {return fFilter3Height;};
     G4double           GetDDGeneratorHeight()      {return fDDGeneratorHeight;};
     G4double           GetInsulatorHeight()        {return fInsulatorHeight;};
     G4double           GetThermalAbsorberHeight()  {return fAbsorberHeight;};
     G4double           GetPortHeight()  						{return fPortHeight;};
     G4double           GetThermalAbsorborHeight()	{return fAbsorberHeight;};
     G4double						GetClearanceAboveCryostat() {return fClearanceAboveCryostat;};
     G4double						GetSourceCenterX()					{return fCenterX;}
     G4double						GetSourceCenterY()					{return fCenterY;}
     G4double						GetInsulatorThickness()     {return fInsulatorThickness;}
     
     
     void               PrintParameters();
                       
  private:
  
     // world box
     G4VPhysicalVolume* fPWorld;
     G4LogicalVolume*   fLWorld;
     
     G4double           fWorldSize;
     G4Material*        fMaterial;   
     
     G4Material*        fWorldMater;  
     
     // height above cryostat
     G4double						fClearanceAboveCryostat;
     // neutron source centers
     G4double						fCenterX;
     G4double						fCenterY;
     
     
     
     // polypropylene insulator
     G4double           fInsulatorThickness; 
     G4double 					fInsulatorLength;
     G4double           fInsulatorWidth;
     G4double 					fInsulatorHeight;
     G4LogicalVolume*   fLogicInsulator;
     G4VPhysicalVolume* fPhysiInsulator;
     G4Material*        fInsulatorMater;
     
     // Feedthrough port
     G4double						fPortHeight;
  	 G4double						fPortOuterRadius;
  	 G4LogicalVolume*   fLogicPort;
     G4VPhysicalVolume* fPhysiPort;
     G4Material*        fPortMater;
    
     // Feedthrough port reflector
     G4double fPortRefHeight; 
     G4double fPortRefThickness;
     G4double fPortRefOuterRadius; 
     G4double fPortRefInnerRadius;
     G4LogicalVolume* fLogicPortRef;
     G4VPhysicalVolume* fPhysiPortRef;
     G4Material* fPortRefMater;
     
     // stainless steel cryostat
     G4double           fCryostatThickness; 
     G4double fCryostatLength;
  	 G4double fCryostatWidth ;
  	 G4double fCryostatHeight;
     G4LogicalVolume*   fLogicCryostat;
     G4VPhysicalVolume* fPhysiCryostat;
     G4Material*        fCryostatMater;  
     
     
     
     // liquid argon pool
     G4double           fPoolHeight; 
     G4double           fPoolWidth;
     G4double           fPoolLength;
     G4LogicalVolume*   fLogicPool;
     G4VPhysicalVolume* fPhysiPool;
     G4Material*        fPoolMater; 
     
     // gas argon
     G4double           fBufferHeight; 
     G4double           fBufferWidth;
     G4double           fBufferLength;
     G4LogicalVolume*   fLogicBuffer;
     G4VPhysicalVolume* fPhysiBuffer;
     G4Material*        fBufferMater; 
     
     // neutron DD generator
     G4double           fDDGeneratorHeight; //UI
     G4double           fDDGeneratorRadius; //UI
     G4LogicalVolume*   fLogicDDGenerator;
     G4VPhysicalVolume* fPhysiDDGenerator;
     G4Material*        fDDGeneratorMater;  
     
     // Moderator: 
     G4double						fModerator_Thickness;
     G4double           fModerator_Height;
     G4double           fModerator_Radius;
     G4LogicalVolume*   fLogicModerator;
     G4VPhysicalVolume* fPhysiModerator;
     G4Material*        fModerator_Mater; 
     
     // Filter1: 
     G4double           fFilter1Height;
     G4double           fFilter1Radius;
     G4LogicalVolume*   fLogicFilter1;
     G4VPhysicalVolume* fPhysiFilter1;
     G4Material*        fFilter1Mater;

     // Filter2: 
     G4double           fFilter2Height;
     G4double           fFilter2Radius_top;
     G4double           fFilter2Radius_bottom;
     G4LogicalVolume*   fLogicFilter2;
     G4VPhysicalVolume* fPhysiFilter2;
     G4Material*        fFilter2Mater;

     // Filter3:
     G4double           fFilter3Height;
  	 G4double           fFilter3Radius_top;
  	 G4double           fFilter3Radius_bottom; 	 
     G4LogicalVolume*   fLogicFilter3;
     G4VPhysicalVolume* fPhysiFilter3;
     G4Material*        fFilter3Mater;
     
     // Thermal neutron absorber
  	 G4double						fAbsorberHeight;
  	 G4double						fThermalAbsorberRadius;
  	 G4LogicalVolume*   fLogicThermalAbsorber;
     G4VPhysicalVolume* fPhysiThermalAbsorber;
     G4Material*        fThermalAbsorberMater;
     
     // neutron reflector
     G4double           fReflectorThickness;
     G4double						fReflectorHeight;
     G4double						fReflectorRadius;
     G4LogicalVolume*   fLogicReflector;
     G4VPhysicalVolume* fPhysiReflector;
     G4Material*        fReflectorMater; 
     
     // neutron guide
     G4double						fNGuideHeight;
     G4double						fNGuideRadius_top;
     G4double						fNGuideRadius_bottom;
     G4LogicalVolume*   fLogicNGuide;
     G4VPhysicalVolume* fPhysiNGuide;
     G4Material*        fNGuideMater; 
     
     // neutron shield
     G4double						fNShieldThickness;
     G4double           fNShieldHeight;
  	 G4double           fNShieldRadius;
     G4LogicalVolume*   fLogicNShield;
     G4VPhysicalVolume* fPhysiNShield;
     G4Material*        fNShieldMater;    
     
     // define vacuum
     G4Material* fVacuum;
     
     DetectorMessenger* fDetectorMessenger;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes(); 
     // color
     G4VisAttributes* blue_;		
     G4VisAttributes* green_;
     G4VisAttributes* red_;		
     G4VisAttributes* grey_;		
     G4VisAttributes* silver_;
     G4VisAttributes* orange_;
     G4VisAttributes* copper_;
     G4VisAttributes* cyan_; 	
     G4VisAttributes* yellow_;  
     G4VisAttributes* pink_;   
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

