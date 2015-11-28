/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                logger.tpp
* Erstellt:                 22.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Definiert die Template-Klasse Logger.
*/

#ifndef __LOGGER_TPP_
#define __LOGGER_TPP_

// nur für intellisense
#include "logger.hpp"

namespace Icarus
{

/**********  Logger **********/

template<typename LogPolicyType>
Logger<LogPolicyType>::Logger(const std::string& name) :
    _pol(new LogPolicyType())
{
    _pol->openLogStream(name);
}

template<typename LogPolicyType>
Logger<LogPolicyType>::~Logger()
{
    if (_pol) _pol->closeLogStream();
}

template<typename LogPolicyType>
template<SeverityType sev, typename...Args>
void Logger<LogPolicyType>::print(unsigned line, std::string file, Args...args)
{
    _mx_log.lock();
    _s.str("");
    _s << " <" << file << ", line " << line << "> ";
    switch (sev)
    {
    case DEBUG:
        _s << "DEBUG: ";
        break;
    case INFO:
        _s << "INFO: ";
        break;
    case WARNING:
        _s << "WARNING: ";
        break;
    case ERROR:
        _s << "ERROR: ";
        break;
    }
    if(sev == ERROR) printRec(true, args...);
    else printRec(false, args...);
    _mx_log.unlock();
}


template<typename LogPolicyType>
template<typename First, typename...Rest>
void Logger<LogPolicyType>::printRec(bool critical, First first, Rest...rest)
{
    _s << first;
    printRec(critical, rest...);
}


template<typename LogPolicyType>
void Logger<LogPolicyType>::printRec(bool critical)
{
    if (critical)
    {
        _pol->write_err(getLogInfo() + _s.str());
        _mx_log.unlock();
        exit(EXIT_LOGFAIL);
    }
    else
    {
        _pol->write(getLogInfo() + _s.str());
    }
}


template <typename LogPolicyType>
std::string Logger<LogPolicyType>::getLogInfo()
{
    std::stringstream info;
    info.fill('0');
    info.width(FIELD_WIDTH);
    info << _lnum++;

    // zeitstempel der nachricht, menschenlesbar
    std::time_t tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    // ms compiler haben probleme mit localtime ("unsafe" trotz standard)
    // gcc < 5 hat put_time noch nicht implementiert
    std::stringstream tsstr;
#if _MSC_VER  > 1700
    std::tm tm;
    std::tm * ptm = &tm;
    localtime_s(ptm,&tt);
    tsstr << std::put_time(ptm, "%c");
#else
    std::tm * ptm = std::localtime(&tt);
    tsstr << ptm->tm_mday + 1 << '/' << ptm->tm_mon + 1 << '/'
          << (1900 + ptm->tm_year)%100 << ' '
          << ptm->tm_hour << ':' << ptm->tm_min << ':' << ptm->tm_sec;
#endif
    info << " [ " << tsstr.str() << " - ";

    // zeitstempel in clock ticks seit programmbeginn
    info.fill('0');
    info.width(FIELD_WIDTH);
    info << clock() << " ]";

    return info.str();
}

}

#endif // __LOGGER_TPP_
