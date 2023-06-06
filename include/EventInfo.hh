#pragma once

#ifndef EVENTINFO_HH
#define EVENTINFO_HH

#include "DetectorInfo.hh"
#include "globals.hh"
using namespace std;



class EventInfo {
public:
    EventInfo(){clear();}
    ~EventInfo(){}
    
    void clear() {
		m_lifetime = 0.;
		m_startTime = 0.;
        m_StartWinkelX = 0.;
        m_StartWinkelY= 0.;
        m_StartWinkelZ= 0.;
        m_StoppWinkel1X= 0.;
        m_StoppWinkel1Y= 0.;
        m_StoppWinkel1Z= 0.;
        m_StoppWinkel2X= 0.;
        m_StoppWinkel2Y= 0.;
        m_StoppWinkel2Z= 0.;
		m_detector1.clear();
		m_detector2.clear();
	}

	void setInfo(double lifetime, double startTime, float StartWinkelX, float StartWinkelY, float StartWinkelZ, float StoppWinkel1X,
    float StoppWinkel1Y, float StoppWinkel1Z, float StoppWinkel2X, float StoppWinkel2Y, float StoppWinkel2Z) {
		m_lifetime = lifetime;
		m_startTime = startTime;
        m_StartWinkelX = StartWinkelX;
        m_StartWinkelY = StartWinkelY;
        m_StartWinkelZ = StartWinkelZ;
        m_StoppWinkel1X = StoppWinkel1X;
        m_StoppWinkel1Y = StoppWinkel1Y;
        m_StoppWinkel1Z = StoppWinkel1Z;
        m_StoppWinkel2X = StoppWinkel2X;
        m_StoppWinkel2Y = StoppWinkel2Y;
        m_StoppWinkel2Z = StoppWinkel2Z;
	}
	
    void attach(const DetectorInfo& detector1, const DetectorInfo& detector2) {
        m_detector1 = detector1;
        m_detector2 = detector2;
    }

    DetectorInfo info1() const {
        return m_detector1;
    }

    DetectorInfo info2() const {
        return m_detector2;
    }
    
    
     double lifetime() const {
        return m_lifetime;
    }

    double startTime() const {
        return m_startTime;
    }

    float StartWinkelX() const {
        return m_StartWinkelX;
    }

    float StartWinkelY() const {
        return m_StartWinkelY;
    }

    float StartWinkelZ() const {
        return m_StartWinkelZ;
    }

    float StoppWinkel1X() const {
        return m_StoppWinkel1X;
    }
    float StoppWinkel1Y() const {
        return m_StoppWinkel1Y;
    }
    float StoppWinkel1Z() const {
        return m_StoppWinkel1Z;
    }

     float StoppWinkel2X() const {
        return m_StoppWinkel2X;
    }
    float StoppWinkel2Y() const {
        return m_StoppWinkel2Y;
    }
    float StoppWinkel2Z() const {
        return m_StoppWinkel2Z;
    }

private:
    DetectorInfo m_detector1;
    DetectorInfo m_detector2;
	
    double m_lifetime;
    double m_startTime;
    float m_StartWinkelX;
    float m_StartWinkelY;
    float m_StartWinkelZ;
    float m_StoppWinkel1X;
    float m_StoppWinkel1Y;
    float m_StoppWinkel1Z;
    float m_StoppWinkel2X;
    float m_StoppWinkel2Y;
    float m_StoppWinkel2Z;
};

#endif
