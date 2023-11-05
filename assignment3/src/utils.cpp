#include "utils.h"
#include <fstream>
#include <iostream>
#include <mpi.h>
#
void save_grid(const std::vector<double> &grid, int rows, int cols, int rank,
               const std::string &grid_type) {
  // Construct the file name with the number of rows and columns
  std::string filename = grid_type + "_" + std::to_string(rows) + "x" +
                         std::to_string(cols) + "_grid.txt";

  // Open the file with the constructed name
  std::ofstream out_file(filename);
  if (!out_file) {
    std::cerr << "Error: Could not open file " << filename << " for writing.\n";
    // Use an MPI call to abort if file cannot be opened
    MPI_Abort(MPI_COMM_WORLD, 1);
    return;
  }

  // Write the grid to the file
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      out_file << grid[i * cols + j] << ' ';
    }
    out_file << '\n';
  }

  out_file.close();
  std::cout << "Output written to " << filename << std::endl;
}