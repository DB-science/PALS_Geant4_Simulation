#include "DG4LTSettingsManager.hh"

static DG4LTSimulationParameterManager* __sharedInstanceDG4LTSettingsManager = nullptr;

DG4LTSimulationParameterManager* DG4LTSimulationParameterManager::sharedInstance() {
	if (!__sharedInstanceDG4LTSettingsManager)
		__sharedInstanceDG4LTSettingsManager = new DG4LTSimulationParameterManager;

	return __sharedInstanceDG4LTSettingsManager;
}

DG4LTSimulationParameterManager::DG4LTSimulationParameterManager() {
	/* lifetime spectrum exponential decays */
	m_lt1Enabled = true;
	m_lt2Enabled = false;
	m_lt3Enabled = false;
	m_lt4Enabled = false;
	m_lt5Enabled = false;

	m_tau1 = 0.106, m_I1 = 0.85; 
	m_tau2 = 0.385, m_I2 = 0.172;
	m_tau3 = 2.750,   m_I3 = 0.003;
	m_tau4 = 0.0,   m_I4 = 0.0;
	m_tau5 = 0.0,   m_I5 = 0.0;

	/* isotope */
	m_isotopeName      = "22Na";
}

DG4LTSimulationParameterManager::~DG4LTSimulationParameterManager() {}

double DG4LTSimulationParameterManager::parseForValue(const std::string& line) {
	std::istringstream iss(line);
	std::vector<std::string> results((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());

	if (results.empty())
		return 0.0;

	const std::string strVal = results.at(results.size() - 1);
	const double val = ::atof(strVal.c_str());

	//G4cout << "line: \"" << line << "\" >> parsed value: \"" << val << "\""<< G4endl;

	return val;
}

std::string DG4LTSimulationParameterManager::parseForString(const std::string& line) {
	std::istringstream iss(line);
	std::vector<std::string> results((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());

	if (results.empty())
		return "";

	const std::string strVal = results.at(results.size() - 1);
	
	//G4cout << "line: \"" << line << "\" >> parsed value: \"" << strVal << "\"" << G4endl;

	return strVal;
}

bool DG4LTSimulationParameterManager::readFile() {
	return registerFile(m_fileName);

}

bool DG4LTSimulationParameterManager::registerFile(const std::string& fileName) {
	std::ifstream infile(fileName);

	G4cout << "\nstart reading input file \"" << fileName << "\"" << G4endl;
	//G4cout << "\n---------------------------------------------------------------- \"" << G4endl;

	if (!infile.is_open()) {
		G4cout << "\nfatal error: file does not exist or cannot be opened!" << G4endl;
		//G4cout << "\n---------------------------------------------------------------- \"" << G4endl;

		return false;
	}

	m_fileName = fileName;

	for (std::string line; std::getline(infile, line); ) {
		if (line.find('!') != std::string::npos)
			continue;

		if (line.find('#') != std::string::npos) {
			// exp. decays
			// 1
			if (line.find("ltDecay1_enabled") != std::string::npos) {
				m_lt1Enabled = parseForValue(line);
			}
			else if (line.find("ltDecay1_tau") != std::string::npos) {
				m_tau1 = parseForValue(line);
				G4cout<<"tau"<<m_tau1<<G4endl;
			}
			else if (line.find("ltDecay1_intensity") != std::string::npos) {
				m_I1 = parseForValue(line);
			}
			// 2
			else if (line.find("ltDecay2_enabled") != std::string::npos) {
				m_lt2Enabled = parseForValue(line);
			}
			else if (line.find("ltDecay2_tau") != std::string::npos) {
				m_tau2 = parseForValue(line);
			}
			else if (line.find("ltDecay2_intensity") != std::string::npos) {
				m_I2 = parseForValue(line);
			}
			// 3
			else if (line.find("ltDecay3_enabled") != std::string::npos) {
				m_lt3Enabled = parseForValue(line);
			}
			else if (line.find("ltDecay3_tau") != std::string::npos) {
				m_tau3 = parseForValue(line);
			}
			else if (line.find("ltDecay3_intensity") != std::string::npos) {
				m_I3 = parseForValue(line);
			}
			// 4
			else if (line.find("ltDecay4_enabled") != std::string::npos) {
				m_lt4Enabled = parseForValue(line);
			}
			else if (line.find("ltDecay4_tau") != std::string::npos) {
				m_tau4 = parseForValue(line);
			}
			else if (line.find("ltDecay4_intensity") != std::string::npos) {
				m_I4 = parseForValue(line);
			}
			// 5
			else if (line.find("ltDecay5_enabled") != std::string::npos) {
				m_lt5Enabled = parseForValue(line);
			}
			else if (line.find("ltDecay5_tau") != std::string::npos) {
				m_tau5 = parseForValue(line);
			}
			else if (line.find("ltDecay5_intensity") != std::string::npos) {
				m_I5 = parseForValue(line);
			}
			// isotope
			else if (line.find("name_of_isotope") != std::string::npos) {
				m_isotopeName = parseForString(line);
			}
		}
	}

	//G4cout << "\n---------------------------------------------------------- \"" << G4endl;

	infile.close();

	return true;
}

