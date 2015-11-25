/*
* Projekt:                  Studienprojekt TM 2015/16
* Dateiname:                mpihandler.hpp
* Erstellt:                 23.11.15
* Autor / Ansprechpartner:  David
*
* Kurzbeschreibung:
* - Diese Singleton-Klasse verwaltet Initialisierung und Deinitialisierung der MPI-Umgebung und
*   stellt die Makrofunktion MPI_SCALL(...) zur Verfügung, die als Wrapper für alle
*   Aufrufe an die MPI-API benutzt werden sollte.
* - Es wird garantiert, dass der Thread, der Init aufruft, auch Finalize aufruft. 
*   MpiSafeCall ist in dieser Form *nicht* Thread-sicher.
* - Das Symbol USE_MPI_ERROR_CHECKING steuert, ob MPI_SCALL wirklich etwas tut (leichter Overhead).
*/

/// \file mpihandler.hpp
/// \author David
/// \brief
/// MpiHandler.

//TODO: Dokumentation

#ifndef __MPIHANDLER_HPP_
#define __MPIHANDLER_HPP_

#define USE_MPI_ERROR_CHECKING

#include <string>
#include "mpi.h"
#include "proto.hpp"

namespace Icarus
{ 
   class MpiHandler : public NonCopyable
    {
        int _n_procs, _my_rank;

    public:       

        static MpiHandler& inst()
        {
            // seit c++11 sollte das threadsicher sein
            // aber es existieren wohl bugs in VS <= 2013
            static MpiHandler _the_instance;
            return _the_instance;
        }

        int get_n_procs() const { return _n_procs; }

        bool is_first() const {return _my_rank == 0;}

        bool is_last() const {return _my_rank == _n_procs-1;}

        int get_my_rank() const { return _my_rank; }

        void MpiSafeCall(int line, std::string file, int error) const;

    private:
        MpiHandler();
        ~MpiHandler();
    };
}

#ifdef USE_MPI_ERROR_CHECKING
#define MPI_SCALL(X) Icarus::MpiHandler::inst().MpiSafeCall(__LINE__,__FILE__,X)
#else
#define MPI_SCALL(...)
#endif

#endif // __MPIHANDLER_HPP_
