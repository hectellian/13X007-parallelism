#include "communications.hpp"
#include <mpi.h>
#include <cmath>
#include <iostream>

/// @brief Sequentially sends data from process 0 to all other processes
void sequential_broadcast(int rank, int size, char *processor_name, MPI_Status stat) {
    int data = rank;  // Initialize data with the rank of the process
    int received_data;
    
    if(rank == 0) {
        // Root process sends data to all other processes
        for(int i = 1; i < size; i++) {
            std::cout << "Broadcast: Rank " << rank << " sending data to rank " << i << std::endl;
            MPI_Send(&received_data, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    } else {
        // All other processes receive data from root process
        MPI_Recv(&data, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
        std::cout << "Broadcast: Rank " << rank << " received data from rank 0" << std::endl;
    }
}

/// @brief Sends data in a ring topology starting and ending at process 0
void sequential_ring(int rank, int size, char *processor_name, MPI_Status stat) {
    int data = rank;  // Initialize data with the rank of the process
    int received_data;
    
    int next = (rank + 1) % size;  // Calculate the next process in the ring
    int prev = (rank + size - 1) % size;  // Calculate the previous process in the ring

    if(rank == 0) {
        // Root process starts the ring
        std::cout << "Ring: Rank " << rank << " sending data to rank " << next << std::endl;
        MPI_Send(&data, 1, MPI_INT, next, 0, MPI_COMM_WORLD);
    } else {
        // Receive data from the previous process
        MPI_Recv(&received_data, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, &stat);
        std::cout << "Ring: Rank " << rank << " received data from rank " << prev << std::endl;
        
        // Send data to the next process
        std::cout << "Ring: Rank " << rank << " sending data to rank " << next << std::endl;
        MPI_Send(&data, 1, MPI_INT, next, 0, MPI_COMM_WORLD);
    }
}

/// @brief Executes hypercube communication among the processes
void hypercube(int rank, int size, char *processor_name, MPI_Status stat) {
    int dim = std::log2(size);  // Calculate the dimension of the hypercube
    if (std::pow(2, dim) != size) {
        if (rank == 0) {
            std::cerr << "Hypercube: The number of processes must be a power of 2." << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);  // Abort if number of processes is not a power of 2
    }
    int data = rank;  // Initialize data with the rank of the process

    for (int i = 0; i < dim; ++i) {
        int partner = rank ^ (1 << i);  // Compute partner rank by XOR-ing with 2^i
        int received_data;

        if (rank < partner) {
            // Lower-ranked process sends first, then receives
            std::cout << "Hypercube: Rank " << rank << " sending data to rank " << partner << std::endl;
            MPI_Send(&data, 1, MPI_INT, partner, 0, MPI_COMM_WORLD);

            MPI_Recv(&received_data, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, &stat);
            std::cout << "Hypercube: Rank " << rank << " received data from rank " << partner << std::endl;
        } else {
            // Higher-ranked process receives first, then sends
            MPI_Recv(&received_data, 1, MPI_INT, partner, 0, MPI_COMM_WORLD, &stat);
            std::cout << "Hypercube: Rank " << rank << " received data from rank " << partner << std::endl;

            std::cout << "Hypercube: Rank " << rank << " sending data to rank " << partner << std::endl;
            MPI_Send(&data, 1, MPI_INT, partner, 0, MPI_COMM_WORLD);
        }

        data += received_data;  // Update the data by adding the received_data
    }
}