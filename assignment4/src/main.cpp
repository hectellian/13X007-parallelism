#include <iostream>
#include <cmath>
#include <omp.h>

// Function to approximate pi using Riemann sum and OpenMP
double approximatePiParallel(unsigned long long numRectangles, int numThreads) {
    double width = 1.0 / numRectangles;
    double sum = 0.0;

    omp_set_num_threads(numThreads); // Set the number of threads

    double startTime = omp_get_wtime(); // Start time measurement

    #pragma omp parallel for reduction(+:sum)
    for (unsigned long long i = 0; i < numRectangles; ++i) {
        double x = (i + 0.5) * width;
        double height = 4.0 / (1.0 + x * x);
        sum += height;
    }

    double endTime = omp_get_wtime(); // End time measurement

    std::cout << "Threads: " << numThreads << ", Time taken: " << (endTime - startTime) << " seconds" << std::endl;

    return sum * width;
}

int main() {
    unsigned long long numRectangles = 1e8; // 10^8 rectangles

    for (int numThreads = 2; numThreads <= pow(2, 8); numThreads *= 2) {
        double pi = approximatePiParallel(numRectangles, numThreads);
        std::cout << "Approximation of Pi with " << numThreads << " threads: " << pi << std::endl;
    }

    return 0;
}
