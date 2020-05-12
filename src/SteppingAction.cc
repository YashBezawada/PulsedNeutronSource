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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "Run.hh"
#include "TrackingAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
                           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, TrackingAction* TrAct, EventAction* event)
: G4UserSteppingAction(),fDetector(det),fTrackingAction(TrAct), fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{

  // Get process information 
  const G4StepPoint* endPoint = step->GetPostStepPoint();
  const G4VProcess* process   = endPoint->GetProcessDefinedStep();
  G4String processName = process->GetProcessName();  
    
  // Count processes  
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->CountProcesses(process);
  
  // Get step information
  const G4StepPoint* prePoint = step->GetPreStepPoint();   
  const G4StepPoint* postPoint = step->GetPostStepPoint();
  const G4VPhysicalVolume* prePhysical = prePoint->GetPhysicalVolume();
  const G4VPhysicalVolume* postPhysical = postPoint->GetPhysicalVolume();
  G4double x = postPoint->GetPosition().x(), y = postPoint->GetPosition().y(), z = postPoint->GetPosition().z(); 		
  G4ThreeVector p = step->GetTrack()->GetMomentumDirection();
  G4double theta = GetTheta(p.z());
  G4double phi = GetPhi(p.x(), p.y());
  G4double preEkin = prePoint->GetKineticEnergy()/CLHEP::MeV;
  G4double postEkin = postPoint->GetKineticEnergy()/CLHEP::MeV;	
  G4double preTime = prePoint->GetGlobalTime()/CLHEP::us;
  G4double postTime = postPoint->GetGlobalTime()/CLHEP::us;
  
  // Get track information
  G4Track* track = step->GetTrack();
  G4double ekin = track->GetKineticEnergy();
  G4double trackl = track->GetTrackLength();
  G4double time = track->GetGlobalTime(); 
  
  // Sanity checks
  if(prePhysical == 0 || postPhysical == 0) return;  // The track does not exist  
  if(prePhysical->GetCopyNo() == -1 && postPhysical->GetCopyNo() == -1) return; // Both steps are in the World
  
  // Get logical volume
  const G4LogicalVolume* preLogical = prePoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  const G4LogicalVolume* postLogical = postPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  	
  // Get particle name
  G4String particleName = step->GetTrack()->GetDefinition()->GetParticleName();

  // incident neutron
  //
  if (step->GetTrack()->GetTrackID() == 1) { 
    G4double ekin  = postPoint->GetKineticEnergy();
    G4double trackl = step->GetTrack()->GetTrackLength();
    G4double time   = step->GetTrack()->GetGlobalTime();           
    fTrackingAction->UpdateTrackInfo(ekin,trackl,time);
  }
  
  // neutron step energy in liquid argon
  if(particleName == "neutron" && preLogical == fDetector->GetLogicPool()) {
    G4AnalysisManager::Instance()->FillH1(6,ekin);	
  }
  
//  // neutron falling into antiresonance
//  if(particleName == "neutron" && postLogical == fDetector->GetLogicPool() ) {
//  	if(ekin>47*CLHEP::keV && ekin<62*CLHEP::keV) {
//	    fEventAction->fIfFallToAntiResonance = 1;
//	    fEventAction->fNStepsInAntiResonance++;
//	    fEventAction->fTrackLengthInAntiResonance += step->GetStepLength();
//	    fEventAction->fEnergyLossInAntiResnance += step->GetTotalEnergyDeposit();
//	  }
//	  if(ekin>54*CLHEP::keV && ekin<60*CLHEP::keV) {
//	  	fEventAction->fIfFallToNarrowerAntiResonance = 1;
//	  }
//  }
  
//  // neutron energy loss per step
//  if(particleName == "neutron" && preLogical == fDetector->GetLogicPool() && postLogical == fDetector->GetLogicPool() ) {
//    G4AnalysisManager::Instance()->FillNtupleDColumn(3, 0, ekin);	
//    G4AnalysisManager::Instance()->FillNtupleDColumn(3, 1, step->GetTotalEnergyDeposit());
//  }
  
  // neutron across boundary
  G4int voltag = 0;
  if(particleName == "neutron" && postPoint->GetStepStatus() == fGeomBoundary) {
  	
    // neutrons from DD generator to moderator
    voltag = 0;
    if(preLogical == fDetector->GetLogicDDGenerator() && postLogical ==  fDetector->GetLogicModerator())	{
    	fEventAction->AddNeutronExit_Generator();
      if(fEventAction->GetNNeutronExit_Generator() == 1) {
      	G4AnalysisManager::Instance()->FillH1(7,ekin);
      	G4AnalysisManager::Instance()->FillNtupleDColumn(0, 0, ekin);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 1, theta);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 2, phi);
      	G4AnalysisManager::Instance()->FillNtupleIColumn(0, 3, voltag); 
  	    G4AnalysisManager::Instance()->AddNtupleRow(0);
      }
    }
    
    // neutrons from moderator to filter 1
    voltag = 1;
    if (preLogical == fDetector->GetLogicModerator() && postLogical == fDetector->GetLogicFilter1()) {
      fEventAction->AddNeutronEnter_Filter1();
      if(fEventAction->GetNNeutronEnter_Filter1() == 1) {
      	G4AnalysisManager::Instance()->FillH1(8,ekin);
      	G4AnalysisManager::Instance()->FillNtupleDColumn(0, 0, ekin);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 1, theta);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 2, phi);
  	    G4AnalysisManager::Instance()->FillNtupleIColumn(0, 3, voltag); 
  	    G4AnalysisManager::Instance()->AddNtupleRow(0); 
      }
      if(ekin >= 50*CLHEP::keV && ekin < 1*CLHEP::MeV && fEventAction->GetNNeutronEnter_Filter1() == 1) run->AddARCount();
    }   
    
    // neutrons from filter 1 to filter 2
    voltag = 2;
    if (preLogical == fDetector->GetLogicFilter1() && postLogical == fDetector->GetLogicFilter2()) {
      fEventAction->AddNeutronEnter_Filter2();
      if(fEventAction->GetNNeutronEnter_Filter2() == 1) {
      	G4AnalysisManager::Instance()->FillH1(9,ekin);
      	G4AnalysisManager::Instance()->FillNtupleDColumn(0, 0, ekin);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 1, theta);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 2, phi);
  	    G4AnalysisManager::Instance()->FillNtupleIColumn(0, 3, voltag); 
  	    G4AnalysisManager::Instance()->AddNtupleRow(0); 
      }
    }
    
    // neutrons from filter 2 to filter 3
    voltag = 3;
    if (preLogical == fDetector->GetLogicFilter2() && postLogical == fDetector->GetLogicFilter3()) {
      fEventAction->AddNeutronEnter_Filter3();
      if(fEventAction->GetNNeutronEnter_Filter3() == 1) {
      	G4AnalysisManager::Instance()->FillH1(10,ekin);
      	G4AnalysisManager::Instance()->FillNtupleDColumn(0, 0, ekin);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 1, theta);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 2, phi);
  	    G4AnalysisManager::Instance()->FillNtupleIColumn(0, 3, voltag); 
  	    G4AnalysisManager::Instance()->AddNtupleRow(0); 
      }
    }
    
    // neutrons from filter 3 to port
    voltag = 4;
    if (preLogical == fDetector->GetLogicFilter3() && postLogical == fDetector->GetLogicPort()) {
      fEventAction->AddNeutronEnter_Port();
      if(fEventAction->GetNNeutronEnter_Port() == 1) {
      	G4AnalysisManager::Instance()->FillH1(11,ekin);
      	G4AnalysisManager::Instance()->FillNtupleDColumn(0, 0, ekin);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 1, theta);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 2, phi);
  	    G4AnalysisManager::Instance()->FillNtupleIColumn(0, 3, voltag); 
  	    G4AnalysisManager::Instance()->AddNtupleRow(0);
      }
    }
    
    // neutrons from  absorber to cryostat
    voltag = 5;
    if (preLogical == fDetector->GetLogicThermalAbsorber() && postLogical == fDetector->GetLogicCryostat()) {
      fEventAction->AddNeutronEnter_Cryostat();
      if(fEventAction->GetNNeutronEnter_Cryostat() == 1) {
      	G4AnalysisManager::Instance()->FillH1(12,ekin);
      	G4AnalysisManager::Instance()->FillNtupleDColumn(0, 0, ekin);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 1, theta);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 2, phi);
  	    G4AnalysisManager::Instance()->FillNtupleIColumn(0, 3, voltag); 
  	    G4AnalysisManager::Instance()->AddNtupleRow(0);
      }
    }
    
    // neutrons from cryostat memberane to gas argon buffer
    voltag = 6;
    if (preLogical == fDetector->GetLogicCryostat() && postLogical == fDetector->GetLogicBuffer()) {  
      fEventAction->AddNeutronEnter_ArBuffer();        
      if(fEventAction->GetNNeutronEnter_ArBuffer() == 1) {
      	G4AnalysisManager::Instance()->FillH1(13,ekin);
      	G4AnalysisManager::Instance()->FillNtupleDColumn(0, 0, ekin);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 1, theta);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 2, phi);
  	    G4AnalysisManager::Instance()->FillNtupleIColumn(0, 3, voltag); 
  	    G4AnalysisManager::Instance()->AddNtupleRow(0);
      }
    }
    
    // neutrons entering liquid Argon pool
    voltag = 7;
    if (preLogical == fDetector->GetLogicBuffer() && postLogical == fDetector->GetLogicPool()) {  
      fEventAction->AddNeutronEnter_LArPool();      
      if(fEventAction->GetNNeutronEnter_LArPool() == 1) {
      	G4AnalysisManager::Instance()->FillH1(14,ekin);
      	G4AnalysisManager::Instance()->FillNtupleDColumn(0, 0, ekin);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 1, theta);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 2, phi);
  	    G4AnalysisManager::Instance()->FillNtupleIColumn(0, 3, voltag); 
  	    G4AnalysisManager::Instance()->AddNtupleRow(0);
      }
    }
    
    // neutrons exiting liquid argon
    voltag = 8;
    if (preLogical == fDetector->GetLogicPool() && postLogical != fDetector->GetLogicPool()) { 
    	fEventAction->fNNeutronExit_LArPool++; 
    	if(fEventAction->fNNeutronExit_LArPool == 1) {
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 0, ekin);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 1, theta);
        G4AnalysisManager::Instance()->FillNtupleDColumn(0, 2, phi);
  	    G4AnalysisManager::Instance()->FillNtupleIColumn(0, 3, voltag); 
  	    G4AnalysisManager::Instance()->AddNtupleRow(0);
  	  }
    }
  }  // end of neutron crossing
  
  // captured neutron
  if(particleName == "neutron" && processName == "nCapture") {
  	
  	// fill ntuple
    G4AnalysisManager::Instance()->FillNtupleDColumn(1, 0, x/1000); // ID, column, value
    G4AnalysisManager::Instance()->FillNtupleDColumn(1, 1, y/1000); // ID, column, value
    G4AnalysisManager::Instance()->FillNtupleDColumn(1, 2, z/1000); // ID, column, value
    G4AnalysisManager::Instance()->FillNtupleDColumn(1, 3, postTime); // ID, column, value
    G4AnalysisManager::Instance()->FillNtupleDColumn(1, 4, preEkin); // ID, column, value
    if(postLogical == fDetector->GetLogicPool())	{
    	G4AnalysisManager::Instance()->FillNtupleIColumn(1, 5, 1); //in liquid argon
    	// fill histograms
      G4AnalysisManager::Instance()->FillH1(20,time);
      G4AnalysisManager::Instance()->FillH2(0,x, y);
  	  G4AnalysisManager::Instance()->FillH2(1,x, z);	
    }
    else G4AnalysisManager::Instance()->FillNtupleIColumn(1, 5, 0); 
    
    G4AnalysisManager::Instance()->AddNtupleRow(1);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double SteppingAction::GetTheta(G4double dz){
    
    G4double theta = CLHEP::pi - std::acos(dz);
    theta = theta * (180./CLHEP::pi);	

    return theta*CLHEP::degree;
}

G4double SteppingAction::GetPhi(G4double dx, G4double dy) {
    
    G4double phi = 0.;
    if(dx ==0.0) {
      if(dy>0) phi = CLHEP::pi/2;
      else if(dy<0) phi = -CLHEP::pi/2;
      else phi = 0;
      //else std::cout<<"Error: SteppingAction::GetPhi()"<<std::endl;
    }
    else if( dx!=0.0 ){
      phi = atan(dy/dx);
      if( dx<0.0 ){
        if( dy>0.0 ) phi += CLHEP::pi;
        if( dy<0.0 ) phi -= CLHEP::pi;
      }  
      phi = phi * (180./CLHEP::pi);
    }
    return phi*CLHEP::degree;
}