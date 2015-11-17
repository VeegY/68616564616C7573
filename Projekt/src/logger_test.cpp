/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                logger_test.cpp
* Erstellt:                 17.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Demonstration der Logger Klasse. (Noch) *kein* Unit-Test.
*/

#define LOGGING_LEVEL 3
#include "logger.hpp"
#include <iostream>

int main()
{    
    LOG_WARNING("Ich bin eine Warnung und eine echte ", 0, ".");
    LOG_INFO("Ich bin eine Falschinformation mit ", 2, " Zahlen darin!");
    LOG_DEBUG("Ich helfe beim Debugging in ", 2.2, " Prozent der ", "Faelle.");

    try
    {
        LOG_ERROR("Ich bin ein Fehler und der einzige, der Ausnahmen hochwerfen darf.");
    }
    catch (...)
    {
        std::cout << ">>> Ausnahme gefangen wie gewuenscht." << std::endl;
    }

	return EXIT_SUCCESS;
}