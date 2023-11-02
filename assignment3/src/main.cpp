
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <vector>
#include <iomanip>  // For setting the precision of floating-point values in display
#include <mpi.h>

using namespace std;

const int MAX_ITERATIONS = 1e5;
const int GRID_SIZE_X = 10;
const int GRID_SIZE_Y = 10;
const double KAPPA = 0.1;  // Thermal diffusivity
const double DT = 0.01;    // Time step
const double DX = 0.1;     // Space step in x direction
const double DY = 0.1;     // Space step in y direction

// Function to initialize the grid with initial temperatures
void initializeGrid(vector<vector<double>> &grid, int rows, int cols, double initialTemp) {
    grid.resize(rows, vector<double>(cols, initialTemp));
}

// Function to update the grid using the FTCS scheme
void updateGrid(vector<vector<double>> &grid, int rows, int cols) {
    vector<vector<double>> newGrid(rows, vector<double>(cols));

    for (int i = 1; i < rows - 1; ++i) {
        for (int j = 1; j < cols - 1; ++j) {
            newGrid[i][j] = grid[i][j] + KAPPA * DT * (
                (grid[i + 1][j] - 2 * grid[i][j] + grid[i - 1][j]) / (DX * DX) +
                (grid[i][j + 1] - 2 * grid[i][j] + grid[i][j - 1]) / (DY * DY)
            );
        }
    }
    grid = newGrid;
}

// Function to display grid state in the terminal
void displayGrid(const vector<vector<double>> &grid, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            cout << setprecision(2) << grid[i][j] << " ";
        }
        cout << endl;
    }
}

int main(int argc, char* argv[]) {
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Determine the number of rows each process will handle
    int rows_per_process = GRID_SIZE_X / size;

    vector<vector<double>> localTemperature(rows_per_process, vector<double>(GRID_SIZE_Y));

    // Initialize the grid only in the root process
    vector<vector<double>> globalTemperature;
    if (rank == 0) {
        initializeGrid(globalTemperature, GRID_SIZE_X, GRID_SIZE_Y, 0.0);
    }

    // Scatter initial grid data to all processes
    MPI_Scatter(&globalTemperature[0][0], rows_per_process * GRID_SIZE_Y, MPI_DOUBLE,
                &localTemperature[0][0], rows_per_process * GRID_SIZE_Y, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    double startTime = MPI_Wtime();  // Start time measurement

    // Simulation Loop
    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        // Exchange border data with neighboring processes
        if (rank > 0) {
            MPI_Send(&localTemperature[0][0], GRID_SIZE_Y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
        }
        if (rank < size - 1) {
            MPI_Recv(&localTemperature[rows_per_process - 1][0], GRID_SIZE_Y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Update local grid using FTCS scheme
        updateGrid(localTemperature, rows_per_process, GRID_SIZE_Y);

        // Exchange border data again, but this time in the opposite direction
        if (rank < size - 1) {
            MPI_Send(&localTemperature[rows_per_process - 1][0], GRID_SIZE_Y, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
        if (rank > 0) {
            MPI_Recv(&localTemperature[0][0], GRID_SIZE_Y, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (iter % 1000 == 0 && rank == 0) {  // Display every 1000 iterations for rank 0
          displayGrid(localTemperature, rows_per_process, GRID_SIZE_Y);
        }
    }

    double endTime = MPI_Wtime();  // End time measurement
    double elapsedTime = endTime - startTime;  // Compute elapsed time

    // Gather computed grid data back to the root process
    MPI_Gather(&localTemperature[0][0], rows_per_process * GRID_SIZE_Y, MPI_DOUBLE,
               &globalTemperature[0][0], rows_per_process * GRID_SIZE_Y, MPI_DOUBLE,
               0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Total time taken: " << elapsedTime << " seconds." << endl;
        displayGrid(globalTemperature, GRID_SIZE_X, GRID_SIZE_Y);  // Display the final grid state
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
