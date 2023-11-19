#include "utils.h"
#include "vector.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <mpi.h>

using namespace std;

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

void write_to_bmp(int N, kt::vector2D<double>& data, int iter, double minval, double maxval)
{
    unsigned char bmpfileheader[14] = {'B','M',0,0,0,0,0,0,0,0,54,0,0,0};
    unsigned char bmpinfoheader[40] = { 40, 0,  0, 0,
                                0, 0,  0, 0,
                                0, 0,  0, 0,
                                1, 0, 24, 0 };

    int width = N;
    int height = N;

    int padding = (4 - (width * 3) % 4) % 4;

    int datasize = (3 * width + padding) * height;
    int filesize = 54 + datasize;

    vector<unsigned char> img(datasize);

    // A linear interpolation function, used to create the color scheme.
    auto linear = [](double x, double x1, double x2, double y1, double y2)
       {
           return ( (y2-y1) * x + x2*y1 - x1*y2 ) / (x2-x1);
       };

    for(int iX = 0; iX < width; iX++){
        for(int iY = 0; iY < height; iY++){
            // Restrain the value to be plotted to [0, 1]
            double value = ((data[iY][iX] - minval) / (maxval - minval));
            double r = 0., g = 0., b = 0.;
            // For good visibility, use a color scheme that goes from black-blue to black-red.
            if (value <= 1./8.) {
                r = 0.;
                g = 0.;
                b = linear(value, -1./8., 1./8., 0., 1.);
            }
            else if (value <= 3./8.) {
                r = 0.;
                g = linear(value, 1./8., 3./8., 0., 1.);
                b = 1.;
            }
            else if (value <= 5./8.) {
                r = linear(value, 3./8., 5./8., 0., 1.);
                g = 1.;
                b = linear(value, 3./8., 5./8., 1., 0.);
            }
            else if (value <= 7./8.) {
                r = 1.;
                g = linear(value, 5./8., 7./8., 1., 0.);
                b = 0.;
            }
            else {
                r = linear(value, 7./8., 9./8., 1., 0.);
                g = 0.;
                b = 0.;
            }

            r = min(255. * r, 255.);
            g = min(255. * g, 255.);
            b = min(255. * b, 255.);

            img[(iX + iY*width)*3 + iY*padding + 2] = (unsigned char)(r);
            img[(iX + iY*width)*3 + iY*padding + 1] = (unsigned char)(g);
            img[(iX + iY*width)*3 + iY*padding + 0] = (unsigned char)(b);
        }
    }

    bmpfileheader[ 2] = (unsigned char)(filesize      );
    bmpfileheader[ 3] = (unsigned char)(filesize >>  8);
    bmpfileheader[ 4] = (unsigned char)(filesize >> 16);
    bmpfileheader[ 5] = (unsigned char)(filesize >> 24);

    bmpinfoheader[ 4] = (unsigned char)(       width      );
    bmpinfoheader[ 5] = (unsigned char)(       width >>  8);
    bmpinfoheader[ 6] = (unsigned char)(       width >> 16);
    bmpinfoheader[ 7] = (unsigned char)(       width >> 24);
    bmpinfoheader[ 8] = (unsigned char)(       height      );
    bmpinfoheader[ 9] = (unsigned char)(       height >>  8);
    bmpinfoheader[10] = (unsigned char)(       height >> 16);
    bmpinfoheader[11] = (unsigned char)(       height >> 24);
    bmpinfoheader[20] = (unsigned char)(datasize      );
    bmpinfoheader[21] = (unsigned char)(datasize >>  8);
    bmpinfoheader[22] = (unsigned char)(datasize >> 16);
    bmpinfoheader[23] = (unsigned char)(datasize >> 24);

    stringstream filename;
    filename << "T_";
    filename << setfill('0') << setw(4) << iter;
    filename << ".bmp";

    ofstream f;
    f.open(filename.str().c_str(), ios::binary | ios::out);

    f.write((const char*)bmpfileheader, 14);
    f.write((const char*)bmpinfoheader, 40);

    f.write((const char*)&img[0], datasize);
    f.close();
}