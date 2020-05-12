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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 76293 2013-11-08 13:11:23Z gcosmo $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "Run.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction()
{  
                            
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  fNNeutronExit_Generator = 0;
  
  fNNeutronEnter_Moderator = 0;  
  fNNeutronExit_Moderator = 0;   
                              
  fNNeutronEnter_Filter1 = 0;  
  fNNeutronExit_Filter1 = 0;   
                              
  fNNeutronEnter_Filter2 = 0;      
  fNNeutronExit_Filter2 = 0; 

  fNNeutronEnter_Filter3 = 0;      
  fNNeutronExit_Filter3 = 0;      
                              
  fNNeutronEnter_ThermalAbsorber = 0;  
  fNNeutronExit_ThermalAbsorber = 0; 
  
  fNNeutronEnter_Port = 0;  
  fNNeutronExit_Port = 0; 
  
  fNNeutronEnter_Cryostat = 0;
  fNNeutronExit_Cryostat = 0;
  
  fNNeutronEnter_ArBuffer = 0;  
  fNNeutronExit_ArBuffer = 0; 
  
  fNNeutronEnter_LArPool = 0;  
  fNNeutronExit_LArPool = 0;

  fGammaEnter_World = 0;
  fGammaEnter_Shield = 0;
  fGammaEnter_Reflector = 0;
  fGammaEnter_LArPool = 0;
  fGammaExit_LArPool = 0;
 
  fNNeutronExit_Shield = 0;
  fNNeutronEnter_World = 0;
  
  fIfFallToAntiResonance = 0;
  fIfFallToNarrowerAntiResonance = 0;
  fNStepsInAntiResonance = 0;
  fTrackLengthInAntiResonance = 0;
  fEnergyLossInAntiResnance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
//  G4AnalysisManager::Instance()->FillNtupleIColumn(2, 0, fIfFallToAntiResonance); 
//  G4AnalysisManager::Instance()->FillNtupleIColumn(2, 1, fNStepsInAntiResonance); 
//  G4AnalysisManager::Instance()->FillNtupleDColumn(2, 2, fTrackLengthInAntiResonance); 
//  G4AnalysisManager::Instance()->FillNtupleDColumn(2, 3, fEnergyLossInAntiResnance); 
//  G4AnalysisManager::Instance()->FillNtupleIColumn(2, 4, fIfFallToNarrowerAntiResonance); 
//  G4AnalysisManager::Instance()->AddNtupleRow(2);
//----------------------------------------------------------------- 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


