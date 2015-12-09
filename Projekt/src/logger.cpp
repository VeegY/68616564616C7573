/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                logger.cpp
* Erstellt:                 22.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Definiert das globale Logger-Objekt.
* - Klassendefinition der LogPolicies.
*/

#include "include/logger.hpp"

namespace Icarus
{

/**********  StdLogPolicy **********/


 void StdLogPolicy::openLogStream(const std::string &name)
{
    write(">>> Start logging on stdout/stderr.");
}
 void StdLogPolicy::closeLogStream()
{
    write(">>> Stop logging on stdout/stderr.");
}
 void StdLogPolicy::write(const std::string& msg)
{
    std::cout << msg << std::endl;
}
 void StdLogPolicy::write_err(const std::string& msg)
{
    std::cerr << msg << std::endl;
}

/**********  FileLogPolicy **********/


void FileLogPolicy::openLogStream(const std::string &name)
{
    _outfile->open(name);
    // wenn der logger nicht startet, aufgeben
    if (!_outfile->is_open()) exit(EXIT_LOGFAIL);
}

void FileLogPolicy::closeLogStream()
{
    _outfile->close();
}

void FileLogPolicy::write(const std::string& msg)
{
    *(_outfile) << msg << std::endl;
}

void FileLogPolicy::write_err(const std::string& msg)
{
    write(msg);
}

}
