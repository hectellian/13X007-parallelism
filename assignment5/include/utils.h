#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <string>
#include "vector.h"

// Function to save the final grid to a file
void write_to_bmp(int N, kt::vector2D<double>& data, int iter, double minval, double maxval);

#endif // UTILS_H
