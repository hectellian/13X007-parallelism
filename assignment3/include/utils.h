#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>

// Function to save the final grid to a file
void save_grid(const std::vector<double>& grid, int rows, int cols, int rank, const std::string& grid_type = "final");

#endif // UTILS_H
