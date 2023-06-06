#pragma once

//#define __ENABLE_OPTICAL_PHYSICS

#ifdef __ENABLE_OPTICAL_PHYSICS
#define __USE_OPENMP
#endif

#include <limits>

class G4DLTPGenBasicFunctions {
public:
	static bool isEqual(double d1, double d2) {
		return abs(d1 - d2) < std::numeric_limits<double>::epsilon();
	}
};

#define DDELETE_SAFETY(__var__) if (__var__) { delete __var__; __var__ = nullptr; }
