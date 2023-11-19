#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <mpi.h>
#include "utils.h"
#include "vector.h"

// Define the domain size and simulation parameters
int rows = 8;
int cols = 8;
int time_steps = 1e5; // Number of iter

// Function to initialize the temperature grid
void initialize_temperature(kt::vector2D<double> &grid, int rank,
                            int rows_per_rank, int total_rows) {
    // Set the top edge of the first rank to 1
    if (rank == 0) {
        for (int j = 0; j < cols; ++j) {
            grid[0][j] = 1.0;
        }
    }

    // Set the bottom edge of the last rank to 1
    if (rank == (total_rows / rows_per_rank) - 1) {
        for (int j = 0; j < cols; ++j) {
            grid[rows_per_rank - 1][j] = 1.0;
        }
    }

    // Set the left and right edges to 1 for all ranks
    for (int i = 0; i < rows_per_rank; ++i) {
        grid[i][0] = 1.0;       // Left edge
        grid[i][cols - 1] = 1.0; // Right edge
    }
}

// Function to update the temperature grid using the FTCS scheme
void update_temperature(kt::vector2D<double>& grid,
                        std::vector<double>& top_buffer,
                        std::vector<double>& bottom_buffer,
                        int rank, int rows_per_rank, int total_rows) {
  double hx = 1.0 / total_rows; 
  double hy = 1.0 / total_rows; 
  double C = 1.0; 
  double dt = 0.25 * hx * hx / C; 

  double diagx = -2.0 + hx * hx / (2 * C * dt);
  double diagy = -2.0 + hy * hy / (2 * C * dt);
  double weightx = C * dt / (hx * hx);
  double weighty = C * dt / (hy * hy);

  kt::vector2D<double> new_grid = grid; // This will store the new values

  // Update interior cells
  for (int i = 1; i < rows_per_rank - 1; ++i) {
    for (int j = 1; j < cols - 1; ++j) {
      new_grid[i][j] = weightx * (grid[i - 1][j] + grid[i + 1][j] + grid[i][j] * diagx) +
                       weighty * (grid[i][j - 1] + grid[i][j + 1] + grid[i][j] * diagy);
    }
  }

  // Update top interior row with the top buffer if not at domain boundary
  if (rank != 0) {
    int i = 0; // Top row
    for (int j = 1; j < cols - 1; ++j) {
      new_grid[i][j] = weightx * (top_buffer[j] + grid[i + 1][j] + grid[i][j] * diagx) +
                       weighty * (grid[i][j - 1] + grid[i][j + 1] + grid[i][j] * diagy);
    }
  }

  // Update bottom interior row with the bottom buffer if not at domain boundary
  if (rank != (total_rows / rows_per_rank) - 1) {
    int i = rows_per_rank - 1; // Bottom row
    for (int j = 1; j < cols - 1; ++j) {
      new_grid[i][j] = weightx * (grid[i - 1][j] + bottom_buffer[j] + grid[i][j] * diagx) +
                       weighty * (grid[i][j - 1] + grid[i][j + 1] + grid[i][j] * diagy);
    }
  }

  // Swap the new grid with the old grid to prepare for the next iteration
  grid = new_grid;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Check command-line arguments
    if (argc == 3) {
        rows = cols = std::atoi(argv[1]);
        time_steps = std::atoi(argv[2]);
        
        if (rows <= 0 || time_steps <= 0) {
            if (rank == 0) {
                std::cerr << "Invalid grid size or time steps. Please provide positive integers." << std::endl;
            }
            MPI_Finalize();
            return 1;
        }
    } else {
        if (rank == 0) {
            std::cerr << "Usage: " << argv[0] << " <grid_size> <time_steps>" << std::endl;
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
    kt::vector2D<double> local_grid(rows_per_rank, cols);

    // Initialize the local grid with initial temperature values
    initialize_temperature(local_grid, rank, rows_per_rank, rows);

    // Create buffers for non-blocking communications
    std::vector<double> top_buffer(cols, 0.0), bottom_buffer(cols, 0.0);
    MPI_Request top_request, bottom_request;
    
    double start_time = MPI_Wtime();

    // Simulation loop
    for (int step = 0; step < time_steps; ++step) {
        // Perform non-blocking sends and receives for the top and bottom rows
        if (rank > 0) {
            MPI_Isend(&local_grid[0][0], cols, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &top_request);
            MPI_Irecv(top_buffer.data(), cols, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &top_request);
        }
        if (rank < size - 1) {
            MPI_Isend(&local_grid[rows_per_rank - 1][0], cols, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &bottom_request);
            MPI_Irecv(bottom_buffer.data(), cols, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &bottom_request);
        }

        // Update the local grid
        update_temperature(local_grid, top_buffer, bottom_buffer, rank, rows_per_rank, rows);

        // Wait for non-blocking communications to complete
        if (rank > 0) {
            MPI_Wait(&top_request, MPI_STATUS_IGNORE);
        }
        if (rank < size - 1) {
            MPI_Wait(&bottom_request, MPI_STATUS_IGNORE);
        }
    }

    double end_time = MPI_Wtime();
    if (rank == 0) {
        std::cout << "Simulation took " << (end_time - start_time) << " seconds." << std::endl;
    }

    // Gather the final temperature distribution at the root rank
    kt::vector2D<double> final_grid(rows, cols);

    MPI_Gather(local_grid.data(), rows_per_rank * cols, MPI_DOUBLE, final_grid.data(), rows_per_rank * cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Convert the final 1D grid to a 2D grid to save as BMP
    if (rank == 0) {
        write_to_bmp(rows, final_grid, time_steps, 0, 1);
    }

    MPI_Finalize();
    return 0;
}
