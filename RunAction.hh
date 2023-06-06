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
//
/// \file RunAction.hh
/// \brief Definition of the B2::RunAction class
#pragma once
#ifndef RUNACTION_HH
#define RUNACTION_HH

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Run.hh"

#include "G4AnalysisManager.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include "DG4LTSettingsManager.hh"
//#include "DG4LTPulseGenerator.hh"

#include "G4DLTPGenDefines.hh"
#include "MyEventAction.hh"
#include "EventStreamManager.hh"

class G4Run;
class EventInfo;

/// Run action class

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction() override;
	

	
    void BeginOfRunAction(const G4Run* run) override;
    void   EndOfRunAction(const G4Run* run) override;
    
    void addEventInfo( EventInfo info) {
		m_infoVec.push_back(info);
    EventStreamManager::sharedInstance()->writeEvent(&info, sizeof(EventInfo));
	} 
    std::vector<EventInfo> m_infoVec;
    
    
 
};



#endif
