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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id: DetectorMessenger.cc 70755 2013-06-05 12:17:48Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(), 
 fDetector(Det), fTestemDir(0), fDetDir(0), fMaterCmd(0), fSizeCmd(0),
 fIsotopeCmd(0)
{ 
  fTestemDir = new G4UIdirectory("/ddsource/");
  fTestemDir->SetGuidance("commands specific to this example");
  
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/ddsource/det/",broadcast);
  fDetDir->SetGuidance("detector construction commands");
        
  fMaterCmd = new G4UIcmdWithAString("/ddsource/det/setMat",this);
  fMaterCmd->SetGuidance("Select material of the box.");
  fMaterCmd->SetParameterName("choice",false);
  fMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


  //Get and Set Material Functions
  
  fGetMaterialTable = new G4UIcmdWithoutParameter("/ddsource/det/GetMaterialTable", this);
  fGetMaterialTable->SetGuidance("Parameterless, returns a table of materials used.");
  fGetMaterialTable->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fModeratorSetMaterCmd = new G4UIcmdWithAString("/ddsource/det/SetModeratorMaterial",this);
  fModeratorSetMaterCmd->SetGuidance("Select material of the moderator.");
  fModeratorSetMaterCmd->SetParameterName("choice",false);
  fModeratorSetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fModeratorGetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetModeratorMaterial", this);
  fModeratorGetMaterCmd->SetGuidance("Parameterless, returns moderator material.");
  fModeratorGetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFilter1SetMaterCmd = new G4UIcmdWithAString("/ddsource/det/SetFilter1Material",this);
  fFilter1SetMaterCmd->SetGuidance("Select material of the filter 1.");
  fFilter1SetMaterCmd->SetParameterName("choice",false);
  fFilter1SetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFilter1GetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetFilter1Material", this);
  fFilter1GetMaterCmd->SetGuidance("Parameterless, returns filter 1 material.");
  fFilter1GetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFilter2SetMaterCmd = new G4UIcmdWithAString("/ddsource/det/SetFilter2Material",this);
  fFilter2SetMaterCmd->SetGuidance("Select material of the filter 2.");
  fFilter2SetMaterCmd->SetParameterName("choice",false);
  fFilter2SetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFilter2GetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetFilter2Material", this);
  fFilter2GetMaterCmd->SetGuidance("Parameterless, returns filter 2 material.");
  fFilter2GetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFilter3SetMaterCmd = new G4UIcmdWithAString("/ddsource/det/SetFilter3Material",this);
  fFilter3SetMaterCmd->SetGuidance("Select material of the filter 3.");
  fFilter3SetMaterCmd->SetParameterName("choice",false);
  fFilter3SetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFilter3GetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetFilter3Material", this);
  fFilter3GetMaterCmd->SetGuidance("Parameterless, returns Filter material 3.");
  fFilter3GetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAbsorberSetMaterCmd = new G4UIcmdWithAString("/ddsource/det/SetAbsorberMaterial",this);
  fAbsorberSetMaterCmd->SetGuidance("Select material of the Absorber.");
  fAbsorberSetMaterCmd->SetParameterName("choice",false);
  fAbsorberSetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAbsorberGetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetAbsorberMaterial", this);
  fAbsorberGetMaterCmd->SetGuidance("Parameterless, returns Absorber material.");
  fAbsorberGetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fPortRefSetMaterCmd = new G4UIcmdWithAString("/ddsource/det/SetPortRefMaterial",this);
  fPortRefSetMaterCmd->SetGuidance("Select material of the port reflector.");
  fPortRefSetMaterCmd->SetParameterName("choice",false);
  fPortRefSetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fPortRefGetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetPortRefMaterial", this);
  fPortRefGetMaterCmd->SetGuidance("Parameterless, returns material of the port reflector.");
  fPortRefGetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fReflectorSetMaterCmd = new G4UIcmdWithAString("/ddsource/det/SetReflectorMaterial",this);
  fReflectorSetMaterCmd->SetGuidance("Select material of the reflector.");
  fReflectorSetMaterCmd->SetParameterName("choice",false);
  fReflectorSetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fReflectorGetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetReflectorMaterial", this);
  fReflectorGetMaterCmd->SetGuidance("Parameterless, returns material of the reflector.");
  fReflectorGetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);



  //Set Dimension Commands
  fModeratorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetModeratorThickness",this);
  fModeratorThicknessCmd->SetGuidance("Set Height of the moderator");
  fModeratorThicknessCmd->SetParameterName("ModeratorThickness",false);
  fModeratorThicknessCmd->SetRange("ModeratorThickness>0.");
  fModeratorThicknessCmd->SetUnitCategory("Length");
  fModeratorThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFilter1HeightCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetFilter1Height",this);
  fFilter1HeightCmd->SetGuidance("Set Height of the 1st filter");
  fFilter1HeightCmd->SetParameterName("Filter1Height",false);
  fFilter1HeightCmd->SetRange("Filter1Height>0.");
  fFilter1HeightCmd->SetUnitCategory("Length");
  fFilter1HeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFilter2HeightCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetFilter2Height",this);
  fFilter2HeightCmd->SetGuidance("Set Height of the 2nd filter");
  fFilter2HeightCmd->SetParameterName("Filter2Height",false);
  fFilter2HeightCmd->SetRange("Filter2Height>0.");
  fFilter2HeightCmd->SetUnitCategory("Length");
  fFilter2HeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFilter3HeightCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetFilter3Height",this);
  fFilter3HeightCmd->SetGuidance("Set Height of the 3rd filter");
  fFilter3HeightCmd->SetParameterName("Filter3Height",false);
  fFilter3HeightCmd->SetRange("Filter3Height>0.");
  fFilter3HeightCmd->SetUnitCategory("Length");
  fFilter3HeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fAbsorberHeightCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetAbsorberHeight",this);
  fAbsorberHeightCmd->SetGuidance("Set Height of the thermal neutron absorber");
  fAbsorberHeightCmd->SetParameterName("AbsorberHeight",false);
  fAbsorberHeightCmd->SetRange("AbsorberHeight>0.");
  fAbsorberHeightCmd->SetUnitCategory("Length");
  fAbsorberHeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fPortRefThickCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetPortRefThickness",this);
  fPortRefThickCmd->SetGuidance("Set thickness of the port reflector");
  fPortRefThickCmd->SetParameterName("RefThickness",false);
  fPortRefThickCmd->SetRange("RefThickness>0.");
  fPortRefThickCmd->SetUnitCategory("Length");
  fPortRefThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fReflectorThickCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetReflectorThickness",this);
  fReflectorThickCmd->SetGuidance("Set thickness of the reflector");
  fReflectorThickCmd->SetParameterName("RefThickness",false);
  fReflectorThickCmd->SetRange("RefThickness>0.");
  fReflectorThickCmd->SetUnitCategory("Length");
  fReflectorThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fNShieldThickCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetShieldThickness",this);
  fNShieldThickCmd->SetGuidance("Set thickness of the shield");
  fNShieldThickCmd->SetParameterName("ShieldThickness",false);
  fNShieldThickCmd->SetRange("ShieldThickness>0.");
  fNShieldThickCmd->SetUnitCategory("Length");
  fNShieldThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fReflectorThickCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetPortHeight",this);
  fReflectorThickCmd->SetGuidance("Set Clearance Above Cryostat");
  fReflectorThickCmd->SetParameterName("ClearanceAboveCryostat",false);
  fReflectorThickCmd->SetRange("ClearanceAboveCryostat>0.");
  fReflectorThickCmd->SetUnitCategory("Length");
  fReflectorThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fReflectorThickCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetSourceCenterX",this);
  fReflectorThickCmd->SetGuidance("Set Neutron Source X position");
  fReflectorThickCmd->SetParameterName("SourceCenterX",false);
  fReflectorThickCmd->SetRange("SourceCenterX>-100");
  fReflectorThickCmd->SetUnitCategory("Length");
  fReflectorThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fReflectorThickCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetSourceCenterY",this);
  fReflectorThickCmd->SetGuidance("Set Neutron Source Y position");
  fReflectorThickCmd->SetParameterName("SourceCenterY",false);
  fReflectorThickCmd->SetRange("SourceCenterY>-100.");
  fReflectorThickCmd->SetUnitCategory("Length");
  fReflectorThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
 
  
  fSizeCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/SetSize",this);
  fSizeCmd->SetGuidance("Set size of the box");
  fSizeCmd->SetParameterName("Size",false);
  fSizeCmd->SetRange("Size>0.");
  fSizeCmd->SetUnitCategory("Length");
  fSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fIsotopeCmd = new G4UIcommand("/ddsource/det/SetIsotopeMat",this);
  fIsotopeCmd->SetGuidance("Build and select a material with single isotope");
  fIsotopeCmd->SetGuidance("  symbol of isotope, Z, A, density of material");
  //
  G4UIparameter* symbPrm = new G4UIparameter("isotope",'s',false);
  symbPrm->SetGuidance("isotope symbol");
  fIsotopeCmd->SetParameter(symbPrm);
  //      
  G4UIparameter* ZPrm = new G4UIparameter("Z",'i',false);
  ZPrm->SetGuidance("Z");
  ZPrm->SetParameterRange("Z>0");
  fIsotopeCmd->SetParameter(ZPrm);
  //      
  G4UIparameter* APrm = new G4UIparameter("A",'i',false);
  APrm->SetGuidance("A");
  APrm->SetParameterRange("A>0");
  fIsotopeCmd->SetParameter(APrm);  
  //    
  G4UIparameter* densityPrm = new G4UIparameter("density",'d',false);
  densityPrm->SetGuidance("density of material");
  densityPrm->SetParameterRange("density>0.");
  fIsotopeCmd->SetParameter(densityPrm);
  //
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of density");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("g/cm3"));
  unitPrm->SetParameterCandidates(unitList);
  fIsotopeCmd->SetParameter(unitPrm);
  //
  fIsotopeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterCmd;

  //Material Commands
  delete fModeratorSetMaterCmd;
  delete fModeratorGetMaterCmd;
  delete fFilter1SetMaterCmd;
  delete fFilter1GetMaterCmd;
  delete fFilter2SetMaterCmd;
  delete fFilter2GetMaterCmd;
  delete fFilter3SetMaterCmd;
  delete fFilter3GetMaterCmd;
  delete fAbsorberSetMaterCmd;
  delete fAbsorberGetMaterCmd;
  delete fPortRefSetMaterCmd;
  delete fPortRefGetMaterCmd;
  delete fReflectorSetMaterCmd;
  delete fReflectorGetMaterCmd;

  delete fModeratorThicknessCmd;
  delete fFilter1HeightCmd;
  delete fFilter2HeightCmd;
  delete fFilter3HeightCmd;
  delete fPortRefThickCmd;
  delete fReflectorThickCmd;
  delete fNShieldThickCmd;

  delete fSizeCmd;
  delete fIsotopeCmd;
  delete fDetDir;
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fMaterCmd )
   { fDetector->SetWorldMaterial(newValue);}

  
  //Set and Get Material Commands
  if( command == fModeratorSetMaterCmd )
   { fDetector->SetModeratorMaterial(newValue);}

  if( command == fModeratorGetMaterCmd )
   { fDetector->GetModeratorMaterial();}

  if( command == fFilter1SetMaterCmd )
   { fDetector->SetFilter1Material(newValue);}

  if( command == fFilter1GetMaterCmd )
   { fDetector->GetFilter1Material();}

  if( command == fFilter2SetMaterCmd )
   { fDetector->SetFilter2Material(newValue);}

  if( command == fFilter2GetMaterCmd )
   { fDetector->GetFilter2Material();}

  if( command == fFilter3SetMaterCmd )
   { fDetector->SetFilter3Material(newValue);}

  if( command == fFilter3GetMaterCmd )
   { fDetector->GetFilter3Material();}

  if( command == fAbsorberSetMaterCmd )
   { fDetector->SetAbsorberMaterial(newValue);}

  if( command == fAbsorberGetMaterCmd )
   { fDetector->GetAbsorberMaterial();}

  if( command == fPortRefSetMaterCmd )
   { fDetector->SetPortRefMaterial(newValue);}

  if( command == fPortRefGetMaterCmd )
   { fDetector->GetPortRefMaterial();}

  if( command == fReflectorSetMaterCmd )
   { fDetector->SetReflectorMaterial(newValue);}

  if( command == fReflectorGetMaterCmd )
   { fDetector->GetReflectorMaterial();}
   



   
  //Set Dimension Commands
  if( command == fModeratorThicknessCmd )
   { fDetector->SetModeratorThickness(fSizeCmd->GetNewDoubleValue(newValue));}
  
  if( command == fFilter1HeightCmd )
   { fDetector->SetFilter1Height(fSizeCmd->GetNewDoubleValue(newValue));}

  if( command == fFilter2HeightCmd )
   { fDetector->SetFilter2Height(fSizeCmd->GetNewDoubleValue(newValue));}

  if( command == fFilter3HeightCmd )
   { fDetector->SetFilter3Height(fSizeCmd->GetNewDoubleValue(newValue));}
   
  if( command == fAbsorberHeightCmd )
   { fDetector->SetAbsorberHeight(fSizeCmd->GetNewDoubleValue(newValue));}

  if( command == fPortRefThickCmd )
   { fDetector->SetPortRefThickness(fSizeCmd->GetNewDoubleValue(newValue));}

  if( command == fReflectorThickCmd )
   { fDetector->SetReflectorThickness(fSizeCmd->GetNewDoubleValue(newValue));}

  if( command == fNShieldThickCmd )
   { fDetector->SetNShieldThickness(fSizeCmd->GetNewDoubleValue(newValue));}
   if( command == fClearanceAboveCryostatCmd )
   { fDetector->SetPortHeight(fSizeCmd->GetNewDoubleValue(newValue));}     
   
  if( command == fSourceCenterXCmd )
   { fDetector->SetSourceCenterX(fSizeCmd->GetNewDoubleValue(newValue));}     
  if( command == fSourceCenterYCmd )
   { fDetector->SetSourceCenterY(fSizeCmd->GetNewDoubleValue(newValue));}     
  

  if (command == fIsotopeCmd)
   {
     G4int Z; G4int A; G4double dens;
     G4String name, unt;
     std::istringstream is(newValue);
     is >> name >> Z >> A >> dens >> unt;
     dens *= G4UIcommand::ValueOf(unt);
     fDetector->MaterialWithSingleIsotope (name,name,dens,Z,A);
     fDetector->SetWorldMaterial(name);    
   }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......  fFilterGetMaterCmd->SetGuidance("Parameterless, returns Filter material.");
