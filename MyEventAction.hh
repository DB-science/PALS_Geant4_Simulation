#pragma once
/// \file EventAction.hh
/// \brief Definition of the B2::EventAction class

#ifndef MYEVENTACTION_HH
#define MYEVENTACTION_HH



#include "G4UserEventAction.hh"
#include "G4THitsMap.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4AnalysisManager.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorInfo.hh"
#include "EventInfo.hh"


class RunAction;
class MyEventAction : public G4UserEventAction
{
  public:
    MyEventAction(RunAction *runAction, PrimaryGeneratorAction *primaryGenAction);
    ~MyEventAction();

    void  BeginOfEventAction(const G4Event* );
    void    EndOfEventAction(const G4Event* );

   
     void AddEdep (G4int iVol, G4double Edep, G4double time, G4double interActionX, G4double interActionY, G4double interActionZ);
     void AddPos(G4double pos) { fPos += pos;}
    


 private:

 G4THitsMap<G4double>* GetHitsCollection(G4int hcID,
                                          const G4Event* event) const;
  G4double GetSum(G4THitsMap<G4double>* hitsMap) const;

  EventInfo m_eventInfo;
  DetectorInfo m_detector1;
  DetectorInfo m_detector2;
  
  G4double fEdep1,   fEdep2;
  G4double fPos;
  G4double fTime0;
  G4double ArrCF1;
  G4double ArrCF2;
  
  G4double MeanArr1;
  G4double MeanArr2;
  G4int NumberOfCounts1;
  G4int NumberOfCounts2;

  G4double MedianArr1;
  G4double MedianArr2;

  G4double MaxPulse1;
  G4double MaxPulse2;

  G4double finterAction1X, finterAction2X;
  G4double finterAction1Y, finterAction2Y;
  G4double finterAction1Z, finterAction2Z;
 

  G4int fHitID;

   RunAction* m_runAction;
   PrimaryGeneratorAction *m_primaryGenAction;
};

#endif
