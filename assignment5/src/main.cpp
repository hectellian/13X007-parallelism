#include <complex>
#include <vector>
#include <iostream>
#include <omp.h>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "utils.h"
#include "vector.h"

// Function to check if a point is in the Mandelbrot set
int isInMandelbrotSet(std::complex<double> c, int maxIterations) {
    std::complex<double> z = 0;
    int n = 0;
    for (n = 0; n < maxIterations; ++n) {
        if (std::abs(z) > 2.0) {
            break;
        }
        z = z * z + c;
    }
    return n;
}

int main(int argc, char *argv[]) {
    int num_threads = 1; // Default number of threads
    int maxIterations = 256; // Default maximum number of iterations
    // Default coordinates for the fractal
    double frac_tl_x = -2.0;
    double frac_tl_y = -2.0;
    double frac_br_x = 2.0;
    double frac_br_y = 2.0;

    // Process command line arguments
    if (argc >= 2) {
        num_threads = std::atoi(argv[1]);
        num_threads = std::min(num_threads, 256); // Clamp to 256
    }
    if (argc >= 3) {
        maxIterations = std::atoi(argv[2]);
        maxIterations = std::min(maxIterations, 256); // Clamp to 256
    }

    const int width = 1000;
    const int height = 1000;
    kt::vector2D<double> grid(width, height);

    // Scale factors for x and y axis
    double x_scale = (frac_br_x - frac_tl_x) / width;
    double y_scale = (frac_br_y - frac_tl_y) / height;

    // Set the number of threads for OpenMP
    omp_set_num_threads(num_threads);

    double start_time = omp_get_wtime();

    // Parallelize the loop using OpenMP
    #pragma omp parallel for collapse(2)
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            std::complex<double> c(x * x_scale + frac_tl_x, y * y_scale + frac_tl_y);
            grid[y][x] = isInMandelbrotSet(c, maxIterations);
        }
    }

    double end_time = omp_get_wtime();

    // Find min and max values for color scaling
    double minval = std::numeric_limits<double>::max();
    double maxval = std::numeric_limits<double>::min();
    for (size_t y = 0; y < height; ++y) {
        for (size_t x = 0; x < width; ++x) {
            minval = std::min(minval, grid[y][x]);
            maxval = std::max(maxval, grid[y][x]);
        }
    }

    // Save the fractal to a BMP file
    write_to_bmp(width, grid, maxIterations, minval, maxval);

    std::cout << "Execution Time with " << num_threads << " threads: " << (end_time - start_time) << " seconds\n";
    std::cout << "Fractal image saved as BMP.\n";

    return 0;
}
