#include "communications.hpp"
#include <mpi.h>

/// @brief Sends data from process 0 to all other processes
/// @param rank 
/// @param size 
/// @param data 
/// @param count 
/// @param tag1 
/// @param stat 
/// @param processor_name 
void sequential_broadcast(int rank, int size, double *data, int count, int tag1, MPI_Status stat, char processor_name[]) {
    if(rank == 0) {
        for(int i = 1; i < count; i++) {
            MPI_Send(data, size, MPI_DOUBLE, i, tag1, MPI_COMM_WORLD);
            std::cout << "Sent from " << processor_name << " to " << i << std::endl;
        }
    } else {
        MPI_Recv(data, size, MPI_DOUBLE, 0, tag1, MPI_COMM_WORLD, &stat);
        std::cout << rank << " Received at " << processor_name << " from 0" << std::endl;
    }
}

/// @brief Sends data from process 0 back to itself like a ring
/// @param rank 
/// @param size 
/// @param data 
/// @param count 
/// @param tag1 
/// @param stat 
/// @param processor_name 
void sequential_ring(int rank, int size, double *data, int count, int tag1, MPI_Status stat, char processor_name[]) {
    // modular arithmetic
    int next = (rank + 1) % count;
    int prev = (rank + count - 1) % count;

    // 0 -> 1 -> 2 -> 3 -> ... -> 0
    if(rank == 0) {
        MPI_Send(data, size, MPI_DOUBLE, next, tag1, MPI_COMM_WORLD);
        std::cout << "Sent from " << processor_name << " to " << next << std::endl;
    } else {
        MPI_Recv(data, size, MPI_DOUBLE, prev, tag1, MPI_COMM_WORLD, &stat);
        std::cout << rank << " Received at " << processor_name << " from " << prev << std::endl;
        MPI_Send(data, size, MPI_DOUBLE, next, tag1, MPI_COMM_WORLD);
        std::cout << "Sent from " << processor_name << " to " << next << std::endl;
    }
}

/// @brief 
/// @param rank 
/// @param size 
/// @param data 
/// @param count 
/// @param tag1 
/// @param stat 
/// @param processor_name 
void hypercube(int rank, int size, double *data, int count, int tag1, MPI_Status stat, char processor_name[]) {

}