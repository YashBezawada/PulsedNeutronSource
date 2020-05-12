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
/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void AddNeutronExit_Generator() {fNNeutronExit_Generator++;};
    
    void AddNeutronEnter_Moderator() {fNNeutronEnter_Moderator++;};
    void AddNeutronExit_Moderator() {fNNeutronExit_Moderator++;};
    
    void AddNeutronEnter_Filter1() {fNNeutronEnter_Filter1++;};
    void AddNeutronExit_Filter1() {fNNeutronExit_Filter1++;};
    
    void AddNeutronEnter_Filter2() {fNNeutronEnter_Filter2++;};
    void AddNeutronExit_Filter2() {fNNeutronExit_Filter2++;};

    void AddNeutronEnter_Filter3() {fNNeutronEnter_Filter3++;};
    void AddNeutronExit_Filter3() {fNNeutronExit_Filter3++;};
    
    void AddNeutronEnter_ThermalAbsorber() {fNNeutronEnter_ThermalAbsorber++;};
    void AddNeutronExit_ThermalAbsorber() {fNNeutronExit_ThermalAbsorber++;};
    
    void AddNeutronEnter_Port() {fNNeutronEnter_Port++;};
    void AddNeutronExit_Port() {fNNeutronExit_Port++;};
    
    void AddNeutronEnter_Cryostat() {fNNeutronEnter_Cryostat++;};
    void AddNeutronExit_Cryostat() {fNNeutronExit_Cryostat++;};
    
    void AddNeutronEnter_ArBuffer() {fNNeutronEnter_ArBuffer++;};
    void AddNeutronExit_ArBuffer() {fNNeutronExit_ArBuffer++;};
    
    void AddNeutronEnter_LArPool() {fNNeutronEnter_LArPool++;};
    void AddNeutronExit_LArPool() {fNNeutronExit_LArPool++;};

    void AddGammaEnter_World() {fGammaEnter_World++;};
    void AddGammaEnter_Shield() {fGammaEnter_Shield++;};
    void AddGammaExit_Shield() {fGammaExit_Shield++;};
    void AddGammaEnter_Reflector() {fGammaEnter_Reflector++;};
    void AddGammaEnter_LArPool() {fGammaEnter_LArPool++;};
    void AddGammaExit_LArPool() {fGammaExit_LArPool++;};

    void AddNeutronExit_Shield() {fNNeutronExit_Shield++;};
    void AddNeutronEnter_World() {fNNeutronEnter_World++;};

    G4int GetNNeutronExit_Generator() {return fNNeutronExit_Generator;};
    
    G4int GetNNeutronEnter_Moderator() {return fNNeutronEnter_Moderator;};
    G4int GetNNeutronExit_Moderator() {return fNNeutronExit_Moderator;};

    G4int GetNNeutronEnter_Filter1() {return fNNeutronEnter_Filter1;};
    G4int GetNNeutronExit_Filter1() {return fNNeutronExit_Filter1;};

    G4int GetNNeutronEnter_Filter2() {return fNNeutronEnter_Filter2;};
    G4int GetNNeutronExit_Filter2() {return fNNeutronExit_Filter2;};

    G4int GetNNeutronEnter_Filter3() {return fNNeutronEnter_Filter3;};
    G4int GetNNeutronExit_Filter3() {return fNNeutronExit_Filter3;};

    G4int GetNNeutronEnter_ThermalAbsorber() {return fNNeutronEnter_ThermalAbsorber;};
    G4int GetNNeutronExit_ThermalAbsorber() {return fNNeutronExit_ThermalAbsorber;};
    
    G4int GetNNeutronEnter_Port() {return fNNeutronEnter_Port;};
    G4int GetNNeutronExit_Port() {return fNNeutronExit_Port;};
    
    G4int GetNNeutronEnter_Cryostat() {return fNNeutronEnter_Cryostat;};
    G4int GetNNeutronExit_Cryostat() {return fNNeutronExit_Cryostat;};
    
    G4int GetNNeutronEnter_ArBuffer() {return fNNeutronEnter_ArBuffer;};
    G4int GetNNeutronExit_ArBuffer() {return fNNeutronExit_ArBuffer;};
    
    G4int GetNNeutronEnter_LArPool() {return fNNeutronEnter_LArPool;};
    G4int GetNNeutronExit_LArPool() {return fNNeutronExit_LArPool;};

    G4int GetGammaEnter_World() {return fGammaEnter_World;};
    G4int GetGammaEnter_Shield() {return fGammaEnter_Shield;};
    G4int GetGammaExit_Shield() {return fGammaExit_Shield;};
    G4int GetGammaEnter_Reflector() {return fGammaEnter_Reflector;};
    G4int GetGammaEnter_LArPool() {return fGammaEnter_LArPool;};
    G4int GetGammaExit_LArPool() {return fGammaExit_LArPool;};

    G4int GetNNeutronExit_Shield() {return fNNeutronExit_Shield;};
    G4int GetNNeutronEnter_World() {return fNNeutronEnter_World;};   
    
    G4int			  fIfFallToAntiResonance;
    G4int	      fIfFallToNarrowerAntiResonance;
    G4int				fNStepsInAntiResonance;
    G4double		fTrackLengthInAntiResonance;
    G4double		fEnergyLossInAntiResnance; 
                
  public:                  
    G4int       fNNeutronExit_Generator;

    G4int				fNNeutronEnter_Moderator;
    G4int				fNNeutronExit_Moderator;
    
    G4int				fNNeutronEnter_Filter1;
    G4int				fNNeutronExit_Filter1;
    
    G4int				fNNeutronEnter_Filter2;
    G4int				fNNeutronExit_Filter2;

    G4int				fNNeutronEnter_Filter3;
    G4int				fNNeutronExit_Filter3;
    
    G4int				fNNeutronEnter_ThermalAbsorber;
    G4int				fNNeutronExit_ThermalAbsorber;  
    
    G4int				fNNeutronEnter_Port;
    G4int				fNNeutronExit_Port;  
    
    G4int				fNNeutronEnter_Cryostat;
    G4int				fNNeutronExit_Cryostat;
    
    G4int				fNNeutronEnter_ArBuffer;
    G4int				fNNeutronExit_ArBuffer;  
    
    G4int				fNNeutronEnter_LArPool;
    G4int				fNNeutronExit_LArPool;

    G4int               fGammaEnter_World;
    G4int               fGammaEnter_Shield;
    G4int               fGammaExit_Shield;
    G4int               fGammaEnter_Reflector;
    G4int               fGammaEnter_LArPool;
    G4int               fGammaExit_LArPool;
    
    G4int               fNNeutronExit_Shield;
    G4int               fNNeutronEnter_World;
    
    

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
