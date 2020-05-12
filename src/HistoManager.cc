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
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 67909 2013-03-12 18:51:09Z vnivanch $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "g4root.hh"
#include <math.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("DDsource.root")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     //enable inactivation of histograms
  
  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;
  G4double xmin = -180.;
  G4double xmax = 180.; 
  
  //============== 1D Histos =================
  
  //Histogram 0 - nb of collisions above 1 eV
  G4int ih = analysisManager->CreateH1("nCollision1", "incident neutron: nb of collisions above 1 eV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 1 - total track length above 1 eV
  ih = analysisManager->CreateH1("tracklength1", "incident neutron: total track length above 1 eV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 2 - time of light above 1 eV
  ih = analysisManager->CreateH1("tof1", "incident neutron: time of flight above 1 eV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 3 - nb of collisions bellow 1 eV
  ih = analysisManager->CreateH1("nCollision2", "incident neutron: nb of collisions below 1 eV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 4 - total track length below 1 eV
  ih = analysisManager->CreateH1("tracklength2", "incident neutron: total track length below 1 eV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 5 - time of flight below 1 eV
  ih = analysisManager->CreateH1("tof2", "incident neutron: time of flight below 1 eV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 6 - energy distribution for all steps
  ih = analysisManager->CreateH1("eStep", "Neutron energy distribution for all steps", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 7 - energy spectrum of neutrons from generator
  ih = analysisManager->CreateH1("eyGen", "Energy Spectrum from Generator", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);
  
  //Histogram 8 - energy spectrum of neutrons from moderator to filter1
  ih = analysisManager->CreateH1("eModerator", "Energy Spectrum from Moderator", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 9 - energy spectrum of neutrons from filter 1 to filter2
  ih = analysisManager->CreateH1("eFilter1",  "Energy Spectrum from Filter 1", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 10 - energy spectrum of neutrons from filter 2 to filter3
  ih = analysisManager->CreateH1("eFilter2", "Energy Spectrum from Filter 2", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 11 - energy spectrum of neutrons from filter 3 to port
  ih = analysisManager->CreateH1("eFilter3", "Energy Spetrum from Fileter 3 to port", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);
  
  //Histogram 12 - energy spectrum of neutrons from thermal absorber to cryostat
  ih = analysisManager->CreateH1("eAbsorber", "Energy Spectrum from Thermal Absorber to cryostat", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);
  
  //Histogram 13 - energy spectrum of neutrons from cryostat membrane to gas argon buffer
  ih = analysisManager->CreateH1("eMembrane", "Energy Spectrum from cryostat membrane to argon gas buffer", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 14 - energy spectrum of neutrons entering liquid argon pool (from gas buffer)
  ih = analysisManager->CreateH1("eArbuffer", "Neutrons entering liquid Argon pool", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);
  
  //Histogram 15 - DD gun neutron initial momentum angle
  ih = analysisManager->CreateH1("AngleGen", "DD gun neutrons initial momentum angle", 360, xmin, xmax, "degree");
  analysisManager->SetH1Activation(ih, true);
  
  //Histogram 16 - Theta angle for neutrons leaving the absorber
  ih = analysisManager->CreateH1("ThetaAbsorber", "Theta angle for neutrons leaving the Li absorber", 360, xmin, xmax, "degree");
  analysisManager->SetH1Activation(ih, true);
  
  //Histogram 17 - Phi angle for neutrons leaving the absorber
  ih = analysisManager->CreateH1("PhiAbsorber", "Phi angle for neutrons leaving the thermal absorber", 360, xmin, xmax, "degree");
  analysisManager->SetH1Activation(ih, true);
  
  //Histogram 18 - Theta angle for neutrons entering LAr Pool
  ih = analysisManager->CreateH1("ThetaToLAr", "Theta angle for neutrons entering LAr Pool", 360, -180, 180, "degree");
  analysisManager->SetH1Activation(ih, true);

  //Histogram 19 - Phi angle for neutrons entering LAr Pool
  ih = analysisManager->CreateH1("PhiToLAr", "Phi angle for neutrons entering LAr Pool", 360, -180, 180, "degree");
  analysisManager->SetH1Activation(ih, true);
  
  //Histogram 20 - neutron capturer time
  ih = analysisManager->CreateH1("nCapTime", "Neutron capture time", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);
  
  
  //============== 2D Histos =================
  
  ih = analysisManager->CreateH2("nCap_y_x","neutron capture position in in LAr TPC (top view, y:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);
  ih = analysisManager->CreateH2("nCap_z_x","neutron capture position in LAr TPC (side view, z:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);
  
  //============== nTuples =================
  
  // ID=0, neutron transport
  analysisManager->CreateNtuple("ntran", "Neutron transport"); //id = 0
  analysisManager->CreateNtupleDColumn("ekin");//column 0
  analysisManager->CreateNtupleDColumn("theta");//column 1
  analysisManager->CreateNtupleDColumn("phi");//column 2
  analysisManager->CreateNtupleIColumn("voltag"); //column 3
  
  // ID=1, neutron capture
  analysisManager->CreateNtuple("ncap", "Neutron captures"); //id = 1
  analysisManager->CreateNtupleDColumn("x");       //column 0
  analysisManager->CreateNtupleDColumn("y");       //column 1
  analysisManager->CreateNtupleDColumn("z");       //column 2
  analysisManager->CreateNtupleDColumn("time");    //column 3
  analysisManager->CreateNtupleDColumn("ekin");    //column 4
  analysisManager->CreateNtupleIColumn("voltag");  //column 5 
  analysisManager->FinishNtuple();  
  
//  // Id=2, neutron in antiresonance
//  analysisManager->CreateNtuple("antiresonance", "Neutron in antiresonance"); //id = 2
//  analysisManager->CreateNtupleIColumn("ifInAR");       //column 0
//  analysisManager->CreateNtupleIColumn("stepsInAR");       //column 1
//  analysisManager->CreateNtupleDColumn("tracklenInAR");       //column 2
//  analysisManager->CreateNtupleDColumn("elossInAR");       //column 3
//  analysisManager->CreateNtupleIColumn("ifInNarrowAR");       //column 4
//  analysisManager->FinishNtuple();  
//  
//  // ID=3, neutron energy loss per step
//  analysisManager->CreateNtuple("eloss", "Neutron eloss"); //id = 3
//  analysisManager->CreateNtupleDColumn("ekin");
//  analysisManager->CreateNtupleDColumn("eloss");  
//  analysisManager->FinishNtuple();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
