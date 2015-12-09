#include "include/mpihandler.hpp"
#include "include/logger.hpp"

namespace Icarus
{

void MpiHandler::MpiSafeCall(int line, std::string file, int error) const
{
    // nichts tun, wenncall erfolgreich
    if(error == MPI_SUCCESS) return;

    // sonst fehlertext aufbereiten
    static char errmsg[MPI_MAX_ERROR_STRING];
    static int eclass, elen;
    MPI_Error_class(error, &eclass);
    MPI_Error_string(error, errmsg, &elen);
    LOG_ERROR_LF(line, file, "Encountered Mpi Error: ", errmsg);
}

MpiHandler::MpiHandler() : _n_procs(0), _my_rank(0)
{
    MPI_Init(NULL,NULL);
    MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
    MpiSafeCall(__LINE__,__FILE__, MPI_Comm_rank(MPI_COMM_WORLD, &_my_rank));
    MpiSafeCall(__LINE__, __FILE__, MPI_Comm_size(MPI_COMM_WORLD, &_n_procs));
    LOG_INFO("Mpi successfully initialized.");
}

MpiHandler::~MpiHandler()
{
    MPI_Finalize();
    LOG_INFO("Mpi successfully finalized.");
}

}
