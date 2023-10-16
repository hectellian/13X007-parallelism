#include <mpi.h>
#include <iostream>
#include <vector>
#include <cstring>

#include "communications.hpp"

int main(int argc, char** argv) {
    // var declaration
    int rank, size;
    MPI_Status stat;

    // init MPI environment
    MPI_Init(&argc, &argv);
 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get number of processes

    if (argc < 2) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <hypercube|broadcast|ring>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Communicate between processes using MPI_Send and MPI_Recv
    // Choose communication pattern based on command-line argument
    if (strcmp(argv[1], "hypercube") == 0) {
        hypercube(rank, size, processor_name, stat);
    } else if (strcmp(argv[1], "broadcast") == 0) {
        sequential_broadcast(rank, size, processor_name, stat);
    } else if (strcmp(argv[1], "ring") == 0) {
        sequential_ring(rank, size, processor_name, stat);
    } else {
        if (rank == 0) {
            std::cerr << "Invalid argument. Choose <hypercube|broadcast|ring>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    MPI_Finalize();
    return 0;
}