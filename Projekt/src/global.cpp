/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                global.hpp
* Erstellt:                 24.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Hier (und nur hier) werden die beiden (einzigen) globalen Variablen
*   des Programms definiert: Der Logger und der MpiHandler.
* - Beachte: Sowohl die Reihenfolge als auch die Tatsache, dass beide
*   Definitionen in einer Datei stehen, sind essentiell.
*/

#include "include/logger.hpp"
#include "include/mpihandler.hpp"

Icarus::Logger<Icarus::StdLogPolicy> __log_inst;
Icarus::MpiHandler __mpi_inst;
