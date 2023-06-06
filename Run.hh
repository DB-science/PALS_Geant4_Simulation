#pragma once
#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
 public:
  Run();
  ~Run();

  void SetPrimary(G4ParticleDefinition* particle, G4double energy);

  

  void Merge(const G4Run*) override;
  void EndOfRun();

 private:
  G4ParticleDefinition* fParticle;
  
  
  G4double fEnergy;




};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
