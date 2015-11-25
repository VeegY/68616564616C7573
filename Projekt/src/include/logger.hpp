/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                logger.hpp
* Erstellt:                 17.11.15
* Autor / Ansprechpartner:  David 
*
* Kurzbeschreibung:
* - Definiert die Makros LOG_INFO, LOG_DEBUG, LOG_WARNING, LOG_ERROR zum einfachen, threadsicheren 
*   Logging von MPI- und Nicht-MPI-Anwendungen (Symbol USE_MPI).
* - Der Inhalt der Nachricht wird als Argumentkette übergeben. Es sind alle Typen zulässig, für
*   die eine entsprechende Überladung von std::stringstream::operator<< existiert.
* - Mit dem Symbol DEBUG_LEVEL (0 bis 3) kann eingestellt werden, welche Nachrichten aufgezeichnet 
*   werden. Wenn eine Nachricht nicht aufgezeichnet wird, wird der entsprechende Aufruf bereits zu 
*   Kompilierzeit entfernt, sodass kein Overhead entsteht.
*/

/// \file Logger.hpp
/// \author David
/// \brief
/// Stellt zentrale Loggingfunktionalitäten zur Verfügung.

//TODO: Dokumentation

#ifndef __LOGGER_HPP_
#define __LOGGER_HPP_

#define LOGGING_LEVEL 3

#include <string>
#include <iostream>
#include <sstream>
#include <memory>
#include <fstream>
#include <mutex>
#include <ctime>
#include <chrono>
#include <iomanip>

#include "proto.hpp"

#define EXIT_LOGFAIL -2

namespace Icarus
{
   // interface
    class LogPolicy : public NonCopyable
    {
    public:
        virtual void openLogStream(const std::string& name = "") = 0;

        virtual void closeLogStream() = 0;

        virtual void write(const std::string& msg) = 0;

        virtual void write_err(const std::string& msg) = 0;
    };

    // stdout, stderr
    class StdLogPolicy : public LogPolicy
    {
    public:
        virtual void openLogStream(const std::string &name = "");

        virtual void closeLogStream();

        virtual void write(const std::string& msg);

        virtual void write_err(const std::string& msg);
    };

    // single file logging
    class FileLogPolicy : public LogPolicy
    {
        std::unique_ptr<std::ofstream> _outfile;

    public:
        FileLogPolicy() : _outfile(new std::ofstream())
        { }

        virtual void openLogStream(const std::string &name = "");

        virtual void closeLogStream();

        virtual void write(const std::string& msg);

        virtual void write_err(const std::string& msg);
    };

    // dringlichkeit
    enum SeverityType
    {
        INFO, DEBUG, WARNING, ERROR
    };

    // der logger
    template<typename LogPolicyType>
    class Logger : public NonCopyable
    {
        unsigned _lnum;
        std::stringstream _s;
        std::mutex _mx_log;
        std::unique_ptr<LogPolicyType> _pol;        

    protected:
        static const unsigned FIELD_WIDTH = 7;

    public:
        explicit Logger(const std::string& name = "");

        virtual ~Logger();

        // zentrale schnittstelle des loggers
        template<SeverityType sev, typename...Args>
        void print(unsigned line, std::string file, Args...args);


    protected:
        // erzeuge den ersten teil des log-eintrags
        std::string getLogInfo();

    private:
        // rekursives herausschreiben der nachrichtenteile
        template<typename First, typename...Rest>
        void printRec(bool critical, First first, Rest...rest);

        // rekursionsabbruch
        void printRec(bool critical);
    };
}

#include "logger.tpp"

// gewünschten logger auswählen
extern Icarus::Logger<Icarus::StdLogPolicy> __log_inst;
//extern Icarus::Logger<Icarus::FileLogPolicy> __log_inst("test.log");

#define LOG_ERROR(...) __log_inst.print<Icarus::SeverityType::ERROR>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_ERROR_LF(L,F,...) __log_inst.print<Icarus::SeverityType::ERROR>(L,F,__VA_ARGS__)

#if LOGGING_LEVEL > 2

#define LOG_WARNING(...) __log_inst.print<Icarus::SeverityType::WARNING>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_INFO(...) __log_inst.print<Icarus::SeverityType::INFO>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_DEBUG(...) __log_inst.print<Icarus::SeverityType::DEBUG>(__LINE__,__FILE__,__VA_ARGS__)

#elif LOGGING_LEVEL == 2

#define LOG_WARNING(...) __log_inst.print<Icarus::SeverityType::WARNING>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_INFO(...) __log_inst.print<Icarus::SeverityType::INFO>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_DEBUG(...)

#elif LOGGING_LEVEL == 1

#define LOG_WARNING(...) __log_inst.print<Icarus::SeverityType::WARNING>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_INFO(...)
#define LOG_DEBUG(...)

#else

#define LOG_WARNING(...)
#define LOG_INFO(...)
#define LOG_DEBUG(...)

#endif

#endif // __LOGGER_HPP_
