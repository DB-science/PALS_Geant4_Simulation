#pragma once
#ifndef EVENTSTREAMMANAGER_H
#define EVENTSTREAMMANAGER_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <fstream>
#include "EventInfo.hh"

using namespace std;
//class EventInfo;
#define DDELETE_SAFETY(__var__) if (__var__) { delete __var__; __var__ = nullptr; }

class EventStreamManager {
	EventStreamManager ();
	virtual ~EventStreamManager ();

	ofstream *m_file;
	string m_fileName;
	bool m_isArmed;

	int64_t m_contentInBytes;

public:
	static EventStreamManager *sharedInstance();

	/* this function creates the stream binary file */

	bool start(const string& fileName);

	/* this function closes the stream binary file */

	void stopAndSave();

	/* int 'nBytes': size of event in bytes */

	bool writeEvent(EventInfo  *event, int nBytes);

	inline bool isArmed()    const { return m_isArmed; }
	inline string fileName() const { return m_fileName; }

	int64_t streamedContentInBytes() const { return m_contentInBytes; }
};

#endif // EVENTSTREAMMANAGER_H
