#pragma once
#ifndef DG4LTSETTINGSMANAGER_HH
#define DG4LTSETTINGSMANAGER_HH
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <iostream>
#include <ostream>
#include <iterator>
#include <vector>

#include "G4DLTPGenDefines.hh"

class DG4LTSimulationParameterManager {
	std::string m_fileName;
public:
	/* lifetime spectrum exponential decays */
	bool m_lt1Enabled;
	bool m_lt2Enabled;
	bool m_lt3Enabled;
	bool m_lt4Enabled;
	bool m_lt5Enabled;

	double m_tau1, m_I1;
	double m_tau2, m_I2;
	double m_tau3, m_I3;
	double m_tau4, m_I4;
	double m_tau5, m_I5;

	/* isotope */
	std::string m_isotopeName;

public:
	static DG4LTSimulationParameterManager* sharedInstance();

	bool readFile();
	bool registerFile(const std::string& fileName);

private:
	DG4LTSimulationParameterManager();
	~DG4LTSimulationParameterManager();

	double parseForValue(const std::string& line);
	std::string parseForString(const std::string& line);
};

#endif
