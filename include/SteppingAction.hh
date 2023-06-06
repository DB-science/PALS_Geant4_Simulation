#pragma once
#ifndef STEPPINGACTION_HH
#define STEPPINGACTION_HH

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "MyEventAction.hh"

class SteppingAction : public G4UserSteppingAction
{
public:
	SteppingAction(MyEventAction* eventAction);
	~SteppingAction();

	void UserSteppingAction(const G4Step*);


private:
	MyEventAction *fEventAction;
	G4double InterActionX;
	G4double InterActionY;
	G4double InterActionZ;

};

#endif
