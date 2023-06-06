#pragma once

typedef struct {
	enum Type : int {
		Na22  = 0,
		Co60  = 1,
		Bi207 = 2
	};
} G4DIsotope;

typedef struct {
	G4DIsotope::Type m_isotope;
	double m_activityInBq;
} G4IsotopeInfo;

