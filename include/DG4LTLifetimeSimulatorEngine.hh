#pragma once

#ifndef DG4LTLIFETIMESIMULATORENGINE_HH
#define DG4LTLIFETIMESIMULATORENGINE_HH
#include <random>
#include <array>
#include <chrono>

#include "DG4LTSettingsManager.hh"
#include "G4DLTPGenDefines.hh"

class DG4LTLifetimeSimulatorEngine 
{
    std::default_random_engine m_decayWeightingGenerator;
    std::piecewise_constant_distribution<double> m_decayWeightingDistribution;

    std::default_random_engine m_individualLTGenerator;
    std::exponential_distribution<double> m_distributionLT1;
    std::exponential_distribution<double> m_distributionLT2;
    std::exponential_distribution<double> m_distributionLT3;
    std::exponential_distribution<double> m_distributionLT4;
    std::exponential_distribution<double> m_distributionLT5;

public:
    static DG4LTLifetimeSimulatorEngine* sharedInstance();

    void init(const std::string& isotopeName = "22Na");
    double getRandomLifetime();

private:
    DG4LTLifetimeSimulatorEngine();
    ~DG4LTLifetimeSimulatorEngine();
};

#endif
