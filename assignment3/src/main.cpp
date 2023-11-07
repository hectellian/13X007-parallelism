#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <mpi.h>
#include "utils.h"

// Define the domain size and simulation parameters
int rows = 8;
int cols = 8;
const double alpha = 1e-4; // Thermal diffusivity
const double dx = 0.1;
const double dy = 0.1;
const double dt = 0.01;
const int time_steps = 1e5; // Number of time steps to simulate

// Function to initialize the temperature grid
void initialize_temperature(std::vector<std::vector<double>> &grid, int rank,
                            int rows_per_rank, int total_rows) {
  // The 'total_rows' argument represents the total number of rows in the global grid.
  for (int i = 0; i < rows_per_rank; ++i) {
    for (int j = 0; j < cols; ++j) {
      // Set the top edge of the first rank to 1
      if (i == 0 && rank == 0) {
        grid[i][j] = 1.0;
      } 
      // Set the bottom edge of the last rank to 1
      else if (i == rows_per_rank - 1 && rank == (total_rows / rows_per_rank) - 1) {
        grid[i][j] = 1.0;
      }
      // Set the left and right edges to 1 for all ranks
      else if (j == 0 || j == cols - 1) {
        grid[i][j] = 1.0;
      } 
      // Interior points are set to 0
      else {
        grid[i][j] = 0.0;
      }
    }
  }
}

// Function to update the temperature grid using the FTCS scheme
void update_temperature(std::vector<std::vector<double>> &grid,
                        std::vector<double> &top_buffer,
                        std::vector<double> &bottom_buffer, int rank,
                        int rows_per_rank, int total_rows) {
  std::vector<std::vector<double>> new_grid = grid;
  // Only iterate over the interior cells
  for (int i = 1; i < rows_per_rank - 1; ++i) {
    for (int j = 1; j < cols - 1; ++j) {
      double temp_x =
          (grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx);
      double temp_y =
          (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy);
      new_grid[i][j] = grid[i][j] + alpha * dt * (temp_x + temp_y);
    }
  }

  // Update the top and bottom rows only if they are not the domain boundary
  if (rank != 0) {
    // Update the top interior row with the top buffer values
    int i = 0; // top row
    for (int j = 1; j < cols - 1; ++j) {
      double temp_x = (grid[i + 1][j] - 2.0 * grid[i][j] + top_buffer[j]) / (dx * dx);
      double temp_y = (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy);
      new_grid[i][j] = grid[i][j] + alpha * dt * (temp_x + temp_y);
    }
  }

  if (rank != total_rows / rows_per_rank - 1) {
    // Update the bottom interior row with the bottom buffer values
    int i = rows_per_rank - 1; // bottom row
    for (int j = 1; j < cols - 1; ++j) {
      double temp_x = (bottom_buffer[j] - 2.0 * grid[i][j] + grid[i - 1][j]) / (dx * dx);
      double temp_y = (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / (dy * dy);
      new_grid[i][j] = grid[i][j] + alpha * dt * (temp_x + temp_y);
    }
  }

  // Assign the updated grid back to the original grid reference
  grid = std::move(new_grid);
}


// Main program
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check command-line arguments
    if (argc == 2) {
        rows = cols = std::atoi(argv[1]);
        if (rows <= 0) {
            if (rank == 0) {
                std::cerr << "Invalid grid size. Please provide a positive integer." << std::endl;
            }
            MPI_Finalize();
            return 1;
        }
    } else if (argc > 2) {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <grid_size>" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }

    // Ensure the domain can be evenly distributed
    if (rows % size != 0) {
        if (rank == 0) {
            std::cerr << "Number of rows is not divisible by the size of MPI_COMM_WORLD." << std::endl;
        }
        MPI_Finalize();
        return -1;
    }

    const int rows_per_rank = rows / size;
    std::vector<std::vector<double>> local_grid(rows_per_rank, std::vector<double>(cols, 0.0));

    // Initialize the local grid with initial temperature values
    initialize_temperature(local_grid, rank, rows_per_rank, rows);

    // Convert local grid to a single vector for MPI communication
    std::vector<double> local_grid_linear(rows_per_rank * cols);
    for (int i = 0; i < rows_per_rank; ++i) {
        std::copy(local_grid[i].begin(), local_grid[i].end(), local_grid_linear.begin() + i * cols);
    }

    // Create buffers for non-blocking communications
    
    std::vector<double> top_buffer(cols, 0.0), bottom_buffer(cols, 0.0);
    MPI_Request top_send_request, top_recv_request;
    MPI_Request bottom_send_request, bottom_recv_request;
    
    double start_time = MPI_Wtime();

    // Simulation loop
    for (int step = 0; step < time_steps; ++step) {
        // Perform non-blocking sends and receives for the top and bottom rows
        if (rank > 0) {
            MPI_Isend(local_grid[0].data(), cols, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &top_send_request);
            MPI_Irecv(top_buffer.data(), cols, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &top_recv_request);
        }
        if (rank < size - 1) {
            MPI_Isend(local_grid[rows_per_rank - 1].data(), cols, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &bottom_send_request);
            MPI_Irecv(bottom_buffer.data(), cols, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &bottom_recv_request);
        }

        // Update the local grid
        update_temperature(local_grid, top_buffer, bottom_buffer, rank, rows_per_rank, rows);

        // Wait for non-blocking communications to complete
        if (rank > 0) {
            MPI_Wait(&top_send_request, MPI_STATUS_IGNORE);
            MPI_Wait(&top_recv_request, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Wait(&bottom_send_request, MPI_STATUS_IGNORE);
            MPI_Wait(&bottom_recv_request, MPI_STATUS_IGNORE);
        }
    }

    double end_time = MPI_Wtime();
    if (rank == 0) {
        std::cout << "Simulation took " << (end_time - start_time) << " seconds." << std::endl;
    }

    for (int i = 0; i < rows_per_rank; ++i) {
        std::copy(local_grid[i].begin(), local_grid[i].end(), local_grid_linear.begin() + i * cols);
    }

    // Gather the final temperature distribution at the root rank
    std::vector<double> final_grid;
    if (rank == 0) {
        final_grid.resize(rows * cols);
    }

    MPI_Gather(local_grid_linear.data(), rows_per_rank * cols, MPI_DOUBLE, final_grid.data(), rows_per_rank * cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Convert the final 1D grid to a 2D grid and save as BMP
    if (rank == 0) {
        std::vector<std::vector<double>> final_2d_grid(rows, std::vector<double>(cols));
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                final_2d_grid[i][j] = final_grid[i * cols + j];
            }
        }
        
        write_to_bmp(rows, final_2d_grid, time_steps, 0, 1);
    }

    MPI_Finalize();
    return 0;
}
