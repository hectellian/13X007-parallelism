#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>

// Function to save the final grid to a file
void save_grid(const std::vector<double>& grid, int rows, int cols, int rank, const std::string& grid_type = "final");
void write_to_bmp(int N, std::vector<std::vector<double>>& data, int iter, double minval, double maxval);

#endif // UTILS_H
