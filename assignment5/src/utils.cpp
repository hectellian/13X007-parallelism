#include "utils.h"
#include "vector.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <vector>

using namespace std;

tuple<unsigned char, unsigned char, unsigned char> iteration_to_color(double iteration, double max_iter) {
    double t = iteration / max_iter;
    // Linear interpolation between blue (min) and black (max)
    unsigned char blue = (1.0 - t) * 255;  // Blue fades to black
    return {0, 0, blue};  // RGB value
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

    // Write headers
    for (int i = 0; i < 14; i++) img[i] = bmpfileheader[i];
        for (int i = 0; i < 40; i++) img[14 + i] = bmpinfoheader[i];

    // Write pixel data
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            auto [r, g, b] = iteration_to_color(data[y][x], iter);
            int pos = 54 + (x * 3) + (width * 3 + padding) * (height - y - 1);
            img[pos + 0] = b; // Blue
            img[pos + 1] = g; // Green
            img[pos + 2] = r; // Red
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
    filename << "F_";
    filename << setfill('0') << setw(4) << iter;
    filename << ".bmp";

    ofstream f;
    f.open(filename.str().c_str(), ios::binary | ios::out);

    f.write((const char*)bmpfileheader, 14);
    f.write((const char*)bmpinfoheader, 40);

    f.write(reinterpret_cast<char*>(&img[0]), filesize);
    f.close();
}