#include "EventStreamManager.hh"

EventStreamManager *__sharedInstanceEventStreamManager = nullptr;

EventStreamManager::EventStreamManager() {
    m_file           = nullptr;
    m_isArmed        = false;
    m_fileName       = "";
    m_contentInBytes = 0;
}

EventStreamManager::~EventStreamManager() {
    if (m_file && isArmed())
        stopAndSave();

    DDELETE_SAFETY(m_file);
}

EventStreamManager *EventStreamManager::sharedInstance() {
    if (!__sharedInstanceEventStreamManager)
        __sharedInstanceEventStreamManager = new EventStreamManager();

    return __sharedInstanceEventStreamManager;
}

bool EventStreamManager::start(const string& fileName) {
    if (m_file)
        delete m_file;

    m_file = nullptr;

    if (fileName.empty())
        return false;

    m_fileName = fileName;

    m_file = new ofstream(m_fileName, ios::out | ios::binary);

    m_contentInBytes = 0;
    m_isArmed = true;

    return true;
}

void EventStreamManager::stopAndSave() {
    if (m_file && isArmed())
        m_file->close();

    if (m_file)
        delete m_file;

    m_file = nullptr;

    m_isArmed        = false;
    m_fileName       = "";
    m_contentInBytes = 0;
}

bool EventStreamManager::writeEvent(EventInfo  *event, int nBytes) {
    if (!event ||
            nBytes <= 0)
        return false;

    if (!m_file
            || !m_isArmed)
        return false;

    m_file->write((const char*)event, nBytes);

    m_contentInBytes += nBytes;

    return true;
}
