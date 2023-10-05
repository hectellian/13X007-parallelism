#include <mpi.h>

#ifndef COMMUNICATIONS_HPP
#define COMMUNICATIONS_HPP

void sequential_broadcast(int rank, int size, double *data, int count, int tag1, MPI_Status stat, char processor_name[]);
void sequential_ring(int rank, int size, double *data, int count, int tag1, MPI_Status stat, char processor_name[]);
void hypercube(int rank, int size, double *data, int count, int tag1, MPI_Status stat, char processor_name[]);

#endif // COMMUNICATIONS_HPP
