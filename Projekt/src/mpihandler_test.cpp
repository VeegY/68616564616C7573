/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                mpihandler_test.cpp
* Erstellt:                 23.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Test für den MpiHandler.
*/

#include <iostream>

#include "include/mpihandler.hpp"
#include "include/logger.hpp"

int mpihandler_test()
{
    std::cout << "Hello World from process "
         << Icarus::MpiHandler::inst().get_my_rank()
         << std::endl;
    LOG_INFO("Eine Information.");
    return 0;
}

