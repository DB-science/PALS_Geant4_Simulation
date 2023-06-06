#include "DG4LTLifetimeSimulatorEngine.hh"

static DG4LTLifetimeSimulatorEngine* __sharedInstanceDG4LTLifetimeSimulatorEngine = nullptr;

DG4LTLifetimeSimulatorEngine* DG4LTLifetimeSimulatorEngine::sharedInstance() {
    if (!__sharedInstanceDG4LTLifetimeSimulatorEngine)
        __sharedInstanceDG4LTLifetimeSimulatorEngine = new DG4LTLifetimeSimulatorEngine;

    return __sharedInstanceDG4LTLifetimeSimulatorEngine;
}

DG4LTLifetimeSimulatorEngine::DG4LTLifetimeSimulatorEngine() {}

DG4LTLifetimeSimulatorEngine::~DG4LTLifetimeSimulatorEngine() {
    if (__sharedInstanceDG4LTLifetimeSimulatorEngine)
        delete __sharedInstanceDG4LTLifetimeSimulatorEngine;
}

void DG4LTLifetimeSimulatorEngine::init(const std::string& isotopeName) {
    if (isotopeName == "60Co") {
        DG4LTSimulationParameterManager::sharedInstance()->m_lt1Enabled = true;
        DG4LTSimulationParameterManager::sharedInstance()->m_lt2Enabled = false;
        DG4LTSimulationParameterManager::sharedInstance()->m_lt3Enabled = false;
        DG4LTSimulationParameterManager::sharedInstance()->m_lt4Enabled = false;
        DG4LTSimulationParameterManager::sharedInstance()->m_lt5Enabled = false;

        DG4LTSimulationParameterManager::sharedInstance()->m_I1 = 1.;
        DG4LTSimulationParameterManager::sharedInstance()->m_I2 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_I3 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_I4 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_I5 = 0.;

        DG4LTSimulationParameterManager::sharedInstance()->m_tau1 = 0.003;
        DG4LTSimulationParameterManager::sharedInstance()->m_tau2 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_tau3 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_tau4 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_tau5 = 0.;
    }
    else if (isotopeName == "207Bi") {
        DG4LTSimulationParameterManager::sharedInstance()->m_lt1Enabled = true;
        DG4LTSimulationParameterManager::sharedInstance()->m_lt2Enabled = false;
        DG4LTSimulationParameterManager::sharedInstance()->m_lt3Enabled = false;
        DG4LTSimulationParameterManager::sharedInstance()->m_lt4Enabled = false;
        DG4LTSimulationParameterManager::sharedInstance()->m_lt5Enabled = false;

        DG4LTSimulationParameterManager::sharedInstance()->m_I1 = 1.;
        DG4LTSimulationParameterManager::sharedInstance()->m_I2 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_I3 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_I4 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_I5 = 0.;

        DG4LTSimulationParameterManager::sharedInstance()->m_tau1 = 0.182;
        DG4LTSimulationParameterManager::sharedInstance()->m_tau2 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_tau3 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_tau4 = 0.;
        DG4LTSimulationParameterManager::sharedInstance()->m_tau5 = 0.;
    }

    unsigned int time_ui = static_cast<unsigned int>(time(NULL));

    m_decayWeightingGenerator.seed(time_ui);
    m_individualLTGenerator.seed(time_ui + 1);

    std::array<double, 5> weights {
        (DG4LTSimulationParameterManager::sharedInstance()->m_lt1Enabled ? DG4LTSimulationParameterManager::sharedInstance()->m_I1 : 0.0),
        (DG4LTSimulationParameterManager::sharedInstance()->m_lt2Enabled ? DG4LTSimulationParameterManager::sharedInstance()->m_I2 : 0.0),
        (DG4LTSimulationParameterManager::sharedInstance()->m_lt3Enabled ? DG4LTSimulationParameterManager::sharedInstance()->m_I3 : 0.0),
        (DG4LTSimulationParameterManager::sharedInstance()->m_lt4Enabled ? DG4LTSimulationParameterManager::sharedInstance()->m_I4 : 0.0),
        (DG4LTSimulationParameterManager::sharedInstance()->m_lt5Enabled ? DG4LTSimulationParameterManager::sharedInstance()->m_I5 : 0.0)
    };

    std::array<double, 6> intervals {
        0.0,
        1.0,
        2.0,
        3.0,
        4.0,
        5.0
    };

    m_decayWeightingDistribution = std::piecewise_constant_distribution<double>(intervals.begin(), intervals.end(), weights.begin());

    m_distributionLT1 = std::exponential_distribution<double>(1.0f / DG4LTSimulationParameterManager::sharedInstance()->m_tau1);
    m_distributionLT2 = std::exponential_distribution<double>(1.0f / DG4LTSimulationParameterManager::sharedInstance()->m_tau2);
    m_distributionLT3 = std::exponential_distribution<double>(1.0f / DG4LTSimulationParameterManager::sharedInstance()->m_tau3);
    m_distributionLT4 = std::exponential_distribution<double>(1.0f / DG4LTSimulationParameterManager::sharedInstance()->m_tau4);
    m_distributionLT5 = std::exponential_distribution<double>(1.0f / DG4LTSimulationParameterManager::sharedInstance()->m_tau5);
}

double DG4LTLifetimeSimulatorEngine::getRandomLifetime() {
    const int index = static_cast<int>(m_decayWeightingDistribution(m_decayWeightingGenerator));

    if (index == 0)
        return m_distributionLT1(m_individualLTGenerator);
    else if (index == 1)
        return m_distributionLT2(m_individualLTGenerator);
    else if (index == 2)
        return m_distributionLT3(m_individualLTGenerator);
    else if (index == 3)
        return m_distributionLT4(m_individualLTGenerator);
    else if (index == 4)
        return m_distributionLT5(m_individualLTGenerator);
}
