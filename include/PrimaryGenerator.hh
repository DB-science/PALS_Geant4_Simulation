/// \file PrimaryGenerator.hh
/// \brief Definition of the PrimaryGenerator class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#pragma once
#ifndef PrimaryGenerator_h
#define PrimaryGenerator_h 1

#include "G4VPrimaryGenerator.hh"
#include <random>



typedef enum class DG4LTPGenParticleType : int {
    UNKNOWN = -1,
    START   = 0,
    STOP_1  = 1,
    STOP_2  = 2
};

typedef enum class D4LTGenTrackType : int {
    UNKNWON_VOLUME  = -1,
    DET_A_INELASTIC = 0,
    DET_B_INELASTIC = 1,
    DET_A_ELASTIC   = 2,
    DET_B_ELASTIC   = 3,
};

typedef enum class D4LTGenEventType : int {
    
};

class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGenerator : public G4VPrimaryGenerator
{
  public:
    PrimaryGenerator();    
   ~PrimaryGenerator();

  public:
    virtual void GeneratePrimaryVertex(G4Event* event);
    
    double getTheoLifetime() const;
    double getStartTime() const;
    float getStartAngleX() const;
    float getStoppAngle1X() const;
    float getStoppAngle2X() const;
    float getStartAngleY() const;
    float getStoppAngle1Y() const;
    float getStoppAngle2Y() const;
    float getStartAngleZ() const;
    float getStoppAngle1Z() const;
    float getStoppAngle2Z() const;
    
    inline DG4LTPGenParticleType particleType() const {
        return primaryGenType;
        
    }
    double irand(int min, int max) 
{
    return ((double)rand() / ((double)RAND_MAX + 1.0)) * (max - min) + min;
}
double irand1(int min, int max) 
{
    return ((double)rand() / ((double)RAND_MAX + 1.0)) * (max - min) + min;
}
    G4double fLifetime;
    G4double timeA;
    float StartWinkelX;
    float StartWinkelY;
    float StartWinkelZ;
    float StoppWinkel1X;
    float StoppWinkel1Y;
    float StoppWinkel1Z;
    float StoppWinkel2X;
    float StoppWinkel2Y;
    float StoppWinkel2Z;
  
    
  
   std::default_random_engine distribution_generator;
   std::piecewise_constant_distribution<double> distribution;
   private:
   DG4LTPGenParticleType primaryGenType;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
