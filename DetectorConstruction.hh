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
/// \file DetectorConstruction.hh
/// \brief Definition of the B2a::DetectorConstruction class
#pragma once
#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "tls.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4VSolid.hh"


class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;



class DetectorMessenger;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;

    G4LogicalVolume *GetScoringVolume() const { return fScoringVolume;}

  public:
    G4VPhysicalVolume* Construct() override;
    G4LogicalVolume *logicDetector, *logicFoilO,*logicFoilU,*logicSubstract,*logicDetectorSubstract, *logicFoilL, *logicFoilR, *logicFoilD, *logicDetector1, *logicFoilO1,*logicFoilU1, *logicFoilL1, *logicFoilR1, *logicFoilD1;
    void ConstructSDandField() override;
    G4LogicalVolume *fScoringVolume;
    G4VPhysicalVolume *physFoil;
    // Set methods
   /* void SetTargetMaterial (G4String );
    void SetChamberMaterial(G4String );*/
    void SetMaxStep (G4double );
    void SetCheckOverlaps(G4bool );
  
  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    const G4VPhysicalVolume* GetCAbsorberPV() const;
  const G4VPhysicalVolume* GetNAbsorberPV() const;
  const G4VPhysicalVolume* GetNEAbsorberPV() const;
  const G4VPhysicalVolume* GetEAbsorberPV() const;
  const G4VPhysicalVolume* GetSEAbsorberPV() const;
  const G4VPhysicalVolume* GetSAbsorberPV() const;
  const G4VPhysicalVolume* GetSWAbsorberPV() const;
  const G4VPhysicalVolume* GetWAbsorberPV() const;
  const G4VPhysicalVolume* GetNWAbsorberPV() const;

    G4VPhysicalVolume*   cAbsorberPV; // the absorber physical volume
G4VPhysicalVolume *physFoilO, *physFoilU,*physFoilR, *physFoilL,*physFoilD;

  G4double fWorld_x;
  G4double fWorld_y;
  G4double fWorld_z;

  G4bool fVerbose;
  G4bool fDumpGdml;
 
 // static data members
    static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;
                                         // magnetic field messenger
   
    G4OpticalSurface *mirrorSurface;
    G4UserLimits* fStepLimit = nullptr; // pointer to user step limits

    DetectorMessenger* fMessenger = nullptr; // messenger

    G4bool fCheckOverlaps = true; // option to activate checking of volumes overlaps
};



#endif
