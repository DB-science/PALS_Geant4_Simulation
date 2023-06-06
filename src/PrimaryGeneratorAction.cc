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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the B2::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGenerator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fPrimaryGenerator(0)
{ 
  fPrimaryGenerator = new PrimaryGenerator();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{ 
   delete fPrimaryGenerator;
}

double PrimaryGeneratorAction::getTheoLifetime() const {
	return fPrimaryGenerator->getTheoLifetime();
	
}
    double PrimaryGeneratorAction::getStartTime() const {
		return fPrimaryGenerator->getStartTime();
	}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  fPrimaryGenerator->GeneratePrimaryVertex(event);
}

float PrimaryGeneratorAction::getStartAngleX() const {
	return fPrimaryGenerator->getStartAngleX();;
}
float PrimaryGeneratorAction::getStoppAngle1X() const {
	return fPrimaryGenerator->getStoppAngle1X();
}
float PrimaryGeneratorAction::getStoppAngle2X() const {
	return fPrimaryGenerator->getStoppAngle2X();
}
float PrimaryGeneratorAction::getStartAngleY() const {
	return fPrimaryGenerator->getStartAngleY();
}
float PrimaryGeneratorAction::getStoppAngle1Y() const {
	return fPrimaryGenerator->getStoppAngle1Y();
}
float PrimaryGeneratorAction::getStoppAngle2Y() const {
	return fPrimaryGenerator->getStoppAngle2Y();
}
float PrimaryGeneratorAction::getStartAngleZ() const {
	return fPrimaryGenerator->getStartAngleZ();
}
float PrimaryGeneratorAction::getStoppAngle1Z() const {
	return fPrimaryGenerator->getStoppAngle1Z();
}
float PrimaryGeneratorAction::getStoppAngle2Z() const {
	return fPrimaryGenerator->getStoppAngle2Z();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



