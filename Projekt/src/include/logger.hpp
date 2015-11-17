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

#include <string>
#include <iostream>
#include <sstream>
#include <memory>
#include <fstream>
#include <mutex>
#include <ctime>
#include <chrono>
#include <iomanip>

namespace Icarus
{

    class NonCopyable
    {
    protected:
        NonCopyable() { }
    private:
        NonCopyable(const NonCopyable&) = delete;
        NonCopyable& operator=(const NonCopyable&) = delete;
    };

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
        virtual void openLogStream(const std::string &name = "")
        { 
            write(">>> Start logging on stdout/stderr.");
        }
        virtual void closeLogStream()
        {
            write(">>> Stop logging on stdout/stderr.");
        }
        virtual void write(const std::string& msg)
        {
            std::cout << msg << std::endl;
        }
        virtual void write_err(const std::string& msg)
        {
            std::cerr << msg << std::endl;
        }
    };

    // single file logging
    class FileLogPolicy : public LogPolicy
    {
        std::unique_ptr<std::ofstream> _outfile;

    public:
        FileLogPolicy() : _outfile(new std::ofstream())
        { }

        virtual void openLogStream(const std::string &name = "")
        {
            _outfile->open(name);
            if (!_outfile->is_open())
            {
                throw std::runtime_error("Failed to open specified logfile.");
            }
        }

        virtual void closeLogStream()
        {
            _outfile->close();
        }

        virtual void write(const std::string& msg)
        {
            *(_outfile) << msg << std::endl;
        }

        virtual void write_err(const std::string& msg)
        {
            write(msg);
        }    
    };

#ifdef USE_MPI
    // multi file logging
    class MpiFileLogPolicy : public FileLogPolicy
    {        
        int _my_rank;
    public:
        MpiFileLogPolicy(MpiHandler& mpi) : 
            FileLogPolicy(), 
            _my_rank(0)
        { 
            //TODO: wrapper benutzen
            MPI_Comm_rank(MPI_WORLD, &_my_rank);
        }

        virtual void openLogStream(const std::string &name = "")
        {
            // nummer des prozesses wird vor den dateinamen geschrieben
            std::stringstream name_s << "Proc_" << _my_rank << "_" name;
            FileLogPolicy::openLogStream(name_s.str());            
        }
    };
#endif

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
        explicit Logger(const std::string& name = "") :
            _pol(new LogPolicyType())
        {
            _pol->openLogStream(name);
        }

        virtual ~Logger()
        {
            if (_pol) _pol->closeLogStream();
        }        

        // zentrale schnittstelle des loggers
        template<SeverityType sev, typename...Args>
        void print(unsigned line, std::string file, Args...args)
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


    protected:
        // erzeuge den ersten teil des log-eintrags
        std::string getLogInfo()
        {
            std::stringstream info;
            info.fill('0');
            info.width(FIELD_WIDTH);
            info << _lnum++;  

            // zeitstempel der nachricht, menschenlesbar
            std::time_t tt = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());   
            
            // viele ms compiler haben probleme mit localtime
#ifdef _MSC_VER            
            std::tm tm;
            std::tm * ptm = &tm;
            localtime_s(ptm,&tt);
#else
            std::tm * ptm = std::localtime(&tt);
#endif
            info << " [ " << std::put_time(ptm, "%c")
                << " - ";

            // zeitstempel in clock ticks seit programmbeginn
            info.fill('0');
            info.width(FIELD_WIDTH);
            info << clock() << " ]";

            return info.str();
        }

    private:
        // rekursives herausschreiben der nachrichtenteile
        template<typename First, typename...Rest>
        void printRec(bool critical, First first, Rest...rest)
        {
            _s << first;
            printRec(critical, rest...);
        }

        // rekursionsabbruch
        void printRec(bool critical)
        {
            if (critical)
            {
                _pol->write_err(getLogInfo() + _s.str());
                throw std::runtime_error("A fatal error occured.");
            }
            else
            {
                _pol->write(getLogInfo() + _s.str());
            }
        }
    
    };
}

// gewünschten logger auswählen
static Icarus::Logger<Icarus::StdLogPolicy> __log_inst;
//static Icarus::Logger<Icarus::FileLogPolicy> __log_inst("test.log");

#if LOGGING_LEVEL > 2

#define LOG_ERROR(...) __log_inst.print<Icarus::SeverityType::ERROR>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_WARNING(...) __log_inst.print<Icarus::SeverityType::WARNING>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_INFO(...) __log_inst.print<Icarus::SeverityType::INFO>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_DEBUG(...) __log_inst.print<Icarus::SeverityType::DEBUG>(__LINE__,__FILE__,__VA_ARGS__)

#elif LOGGING_LEVEL == 2

#define LOG_ERROR(...) __log_inst.print<Icarus::SeverityType::ERROR>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_WARNING(...) __log_inst.print<Icarus::SeverityType::WARNING>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_INFO(...) __log_inst.print<Icarus::SeverityType::INFO>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_DEBUG(...)

#elif LOGGING_LEVEL == 1

#define LOG_ERROR(...) __log_inst.print<Icarus::SeverityType::ERROR>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_WARNING(...) __log_inst.print<Icarus::SeverityType::WARNING>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_INFO(...)
#define LOG_DEBUG(...)

#else

#define LOG_ERROR(...) __log_inst.print<Icarus::SeverityType::ERROR>(__LINE__,__FILE__,__VA_ARGS__)
#define LOG_WARNING(...)
#define LOG_INFO(...)
#define LOG_DEBUG(...)

#endif

#endif // __LOGGER_HPP_