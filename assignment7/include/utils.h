#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>

// Function to save the final grid to a fil
void write_to_bmp(int N, std::vector<double>& data, int iter, double minval, double maxval);

#endif // UTILS_H
