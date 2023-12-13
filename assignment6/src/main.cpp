#include "utils.h"
#include "vector.h"
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <sstream>
#include <vector>

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
  int num_threads = 256;   // Default number of threads
  int maxIterations = 256; // Default maximum number of iterations
  bool useStatic = false;  // Default scheduling method
  int chunk_size = 1;      // Default chunk size

  // Default coordinates for the fractal
  double frac_tl_x = -2.0;
  double frac_tl_y = -2.0;
  double frac_br_x = 2.0;
  double frac_br_y = 2.0;

  // 1. Check if no arguments are provided
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " [-t num_threads] [static|dynamic chunk_size]\n";
    return 1;
  }

  // Process command line arguments
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-t" && i + 1 < argc) {
      num_threads = std::atoi(argv[++i]);
      num_threads = std::min(num_threads, 256); // Clamp to 256
    } else if ((arg == "static" || arg == "dynamic") && i + 1 < argc) {
      useStatic = (arg == "static");
      chunk_size = std::atoi(argv[++i]);
    } else if (arg == "-i" && i + 1 < argc) {
      maxIterations = std::atoi(argv[++i]);
    } else if (arg == "-c" && i + 4 < argc) {
      frac_tl_x = std::atof(argv[++i]);
      frac_tl_y = std::atof(argv[++i]);
      frac_br_x = std::atof(argv[++i]);
      frac_br_y = std::atof(argv[++i]);
    } else {
      std::cerr << "Usage: " << argv[0]
                << " [-t num_threads] [static|dynamic chunk_size]\n";
      return 1;
    }
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

  omp_sched_t schedType = useStatic ? omp_sched_static : omp_sched_dynamic;
  omp_set_schedule(schedType, chunk_size); // Set the chunk size

#pragma omp parallel for schedule(runtime)
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

  // 2. Include number of iterations in the final print statement
  std::cout << "Execution Time with " << num_threads
            << " threads: " << (end_time - start_time) << " seconds, "
            << "Max Iterations: " << maxIterations << ", "
            << "Mode: " << (useStatic ? "static" : "dynamic") << ", "
            << "Chunk Size: " << chunk_size << std::endl;

  return 0;
}
