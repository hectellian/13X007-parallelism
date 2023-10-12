#include <mpi.h>

#ifndef COMMUNICATIONS_HPP
#define COMMUNICATIONS_HPP

void sequential_broadcast(int rank, int size, char *processor_name, MPI_Status stat);
void sequential_ring(int rank, int size, char *processor_name, MPI_Status stat);
void hypercube(int rank, int size, char *processor_name, MPI_Status stat);

#endif // COMMUNICATIONS_HPP
