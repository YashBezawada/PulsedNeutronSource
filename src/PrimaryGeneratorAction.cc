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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "HistoManager.hh"
#include "Run.hh"
#include <math.h>
#include "g4root.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det, size_t Size)
: G4VUserPrimaryGeneratorAction(),fParticleGun(0), fDetector(det)
{
  size = Size;
  probDist = new G4double[size];
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle); 
  
  // default particle kinematic

  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleEnergy(2.5*MeV);
  G4double x = fDetector->GetSourceCenterX();
  G4double y = fDetector->GetSourceCenterY();
  G4double z = fDetector->GetInsulatorHeight()/2 - fDetector->GetInsulatorThickness() + fDetector->GetPortHeight() + fDetector->GetFilter3Height() + fDetector->GetFilter2Height() 
               + fDetector->GetFilter1Height() + fDetector->GetModeratorHeight() - fDetector->GetDDGeneratorHeight() + 2.0*cm;
  fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));

  fMessenger = new G4GenericMessenger(this,"/primary/", "...doc...");
  fMessenger->DeclareMethod("updateGunPosition", &PrimaryGeneratorAction::UpdateGunPosition, "...doc...");

  populateProbDist();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete [] probDist;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
   G4double theta = std::abs(DDrandom()*degree);
   G4AnalysisManager::Instance()->FillH1(15,theta);
   
// direction generator (J. Wang Feb 5, 2019)
   G4double cosTheta = -std::cos(theta);
   G4double phi = twopi*G4UniformRand();
   G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
   G4double dx = sinTheta*std::cos(phi),
            dy = sinTheta*std::sin(phi),
            dz = cosTheta;
   	
   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1)); // along z direction
   //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dx,dy,dz)); // realistic angular distribution
   fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::UpdateGunPosition() {
  G4double x = fDetector->GetSourceCenterX();
  G4double y = fDetector->GetSourceCenterY();
  G4double z = fDetector->GetInsulatorHeight()/2 - fDetector->GetInsulatorThickness() + fDetector->GetPortHeight() + fDetector->GetFilter3Height() + fDetector->GetFilter2Height() + fDetector->GetFilter1Height() + fDetector->GetModeratorHeight() - fDetector->GetDDGeneratorHeight() + 2.0*cm;
  fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));	
}

void PrimaryGeneratorAction::populateProbDist(){

    G4double thetaMin = -180;
    G4double thetaMax = 180;

    G4double A = 0.237;
    G4double B = -0.051;
    G4double C = -0.130;

    for(G4int i=0; i<size; i++){
        G4double theta = thetaMin + ((thetaMax - thetaMin)/size)*i;
        if(theta < -90 || theta > 90){
            probDist[i] = A + B + C;
        }
        else{
            G4double cosine = std::cos(theta *(pi/180) + (pi/2));
            probDist[i] = A + B * std::pow(cosine,2) + C * std::pow(cosine,4);
        }
    }
}

G4double PrimaryGeneratorAction::DDrandom(){
    G4RandGeneral *randDistribution = new G4RandGeneral(probDist, size);
    G4double random = randDistribution->shoot();
    G4double angle = -180 + (360 * random);
    return angle;
}
