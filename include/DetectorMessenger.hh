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
/// \file DetectorMessenger.hh
/// \brief Definition of the DetectorMessenger class
//
// $Id: DetectorMessenger.hh 67103 2013-01-31 18:18:03Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger: public G4UImessenger
{
  public:
  
    DetectorMessenger(DetectorConstruction* );
   ~DetectorMessenger();
    
    virtual void SetNewValue(G4UIcommand*, G4String);
    
  private:
  
    DetectorConstruction*      fDetector;
    
    G4UIdirectory*             fTestemDir;
    G4UIdirectory*             fDetDir;
    G4UIcmdWithAString*        fMaterCmd;

    //Set and Get Commands
    G4UIcmdWithoutParameter*   fGetMaterialTable;
    G4UIcmdWithAString*        fModeratorSetMaterCmd;
    G4UIcmdWithoutParameter*   fModeratorGetMaterCmd;
    G4UIcmdWithAString*        fFilter1SetMaterCmd;
    G4UIcmdWithoutParameter*   fFilter1GetMaterCmd;
    G4UIcmdWithAString*        fFilter2SetMaterCmd;
    G4UIcmdWithoutParameter*   fFilter2GetMaterCmd;
    G4UIcmdWithAString*        fFilter3SetMaterCmd;
    G4UIcmdWithoutParameter*   fFilter3GetMaterCmd;
    G4UIcmdWithAString*        fAbsorberSetMaterCmd;
    G4UIcmdWithoutParameter*   fAbsorberGetMaterCmd;
    G4UIcmdWithAString*        fPortRefSetMaterCmd;
    G4UIcmdWithoutParameter*   fPortRefGetMaterCmd;
    G4UIcmdWithAString*        fReflectorSetMaterCmd;
    G4UIcmdWithoutParameter*   fReflectorGetMaterCmd;

    //Set Dimesion Commands
    G4UIcmdWithADoubleAndUnit* fModeratorThicknessCmd;
    G4UIcmdWithADoubleAndUnit* fFilter1HeightCmd;
    G4UIcmdWithADoubleAndUnit* fFilter2HeightCmd;
    G4UIcmdWithADoubleAndUnit* fFilter3HeightCmd;
    G4UIcmdWithADoubleAndUnit* fAbsorberHeightCmd;
    G4UIcmdWithADoubleAndUnit* fPortRefThickCmd;
    G4UIcmdWithADoubleAndUnit* fReflectorThickCmd;
    G4UIcmdWithADoubleAndUnit* fNShieldThickCmd;
    G4UIcmdWithADoubleAndUnit* fClearanceAboveCryostatCmd;
    G4UIcmdWithADoubleAndUnit* fSourceCenterXCmd;
    G4UIcmdWithADoubleAndUnit* fSourceCenterYCmd;

    G4UIcmdWithADoubleAndUnit* fSizeCmd;
    G4UIcmdWithADoubleAndUnit* fFilterHeightCmd;
    
    
    G4UIcommand*               fIsotopeCmd;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

