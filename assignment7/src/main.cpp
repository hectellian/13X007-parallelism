#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <execution>
#include <chrono>
#include "utils.h"

using namespace std;
using namespace std::chrono;

const int max_iter = 1e5;
const int Nx = pow(10, 4);
const int Ny = pow(10, 4);

int main(int argc, char* argv[]) {
    std::vector<double> U(Nx * Ny, 0.0);
    std::vector<double> U0(Nx * Ny, 0.0);

    // Initialize boundaries
    for(int i = 0; i < Nx; ++i) {
        U[i * Ny] = 1;  // Left
        U0[i * Ny] = 1;
        U[(i + 1) * Ny - 1] = 1;  // Right
        U0[(i + 1) * Ny - 1] = 1;
    }
    for(int j = 0; j < Ny; ++j) {
        U[j] = 1;  // Top
        U0[j] = 1;
        U[(Nx - 1) * Ny + j] = 1;  // Bottom
        U0[(Nx - 1) * Ny + j] = 1;
    }

    auto hx = 1.0 / Nx;
    auto hy = 1.0 / Ny;
    auto C = 1.0;
    auto dt = 0.25 * hx * hx / C;
    auto diagx = -2.0 + hx * hx / (2 * C * dt);
    auto diagy = -2.0 + hy * hy / (2 * C * dt);
    auto weightx = C * dt / (hx * hx);
    auto weighty = C * dt / (hy * hy);

    auto start = high_resolution_clock::now(); // Start timing

    // Parallel computation
    for(int iT = 0; iT < max_iter; iT++) {
        std::for_each(std::execution::par_unseq, U.begin(), U.end(), [&](double& u_ij) {
            int i = &u_ij - &U[0];
            int x = i / Ny;
            int y = i % Ny;
            if(x > 0 && x < Nx - 1 && y > 0 && y < Ny - 1) {
                u_ij = weightx * (U0[(x-1) * Ny + y] + U0[(x+1) * Ny + y] + U0[x * Ny + y] * diagx)
                     + weighty * (U0[x * Ny + (y-1)] + U0[x * Ny + (y+1)] + U0[x * Ny + y] * diagy);
            }
        });

        swap(U, U0);
    }

    auto stop = high_resolution_clock::now(); // End timing
    auto duration = duration_cast<seconds>(stop - start);

    //Print time
    cout << "Time taken by function: "
         << duration.count() << " seconds" << endl;

    write_to_bmp(Nx, U0, max_iter, 0, 1);

    return 0;
}
