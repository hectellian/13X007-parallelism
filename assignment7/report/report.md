# Introduction

This report evaluates a heat simulation application, focusing on its parallelization approach and execution on different hardware configurations. The application is designed to solve a 2D heat equation, a typical problem in computational fluid dynamics and heat transfer simulations. The challenge lies in efficiently parallelizing the computations to leverage the capabilities of modern CPUs and GPUs. The parallelization strategy involves using standard C++ libraries and execution policies tailored to each hardware's strengths.

# Methodology

## Hardware Employed

- **Yggdrasil HPC:** A formidable configuration of 32 CPUs alongside a single GPU.
- **Personal Computer:** A setup comprising 12 CPU cores and an RTX 2070 Super GPU.

## Implementation Details

- **Data Structure:** The 2D domain of the heat simulation is stored in a linear vector, a common technique in high-performance computing to simplify memory access patterns and enhance data locality.
- **Parallelization Technique:** The application utilizes `std::for_each` combined with execution policies from the C++ Standard Library. This approach enables efficient parallel processing on multi-core CPUs.
- **CPU Execution:** On CPUs, the parallel execution policy (`std::execution::par_unseq`) is employed, allowing for unsequenced, concurrent execution across multiple threads.
- **GPU Execution:** The GPU implementation details were not provided in the initial code snippet. However, it typically involves using GPU-specific programming models like CUDA or OpenCL for parallel processing.

# Results

## Execution Time Ratio

The execution time ratio between CPU and GPU implementations was calculated based on the provided data:

- **Yggdrasil HPC:** CPU (32 cores) execution time was significantly higher than the GPU.
- **Personal Computer:** The execution times for the CPU (12 cores) and GPU (RTX 2070 Super) were almost identical.

## Performance Metrics

1. **On Yggdrasil HPC:**
   - **CPU (32 cores):** Marked an execution time of 41,391,736 microseconds.
   - **GPU:** Recorded an execution time of 53,800,258 microseconds.

2. **On Personal Computer:**
   - **CPU (12 cores):** Notched an execution time of 27,182,565 microseconds.
   - **GPU (RTX 2070 Super):** Clocking in at 25,766,272 microseconds.

## Graphical Representation

A bar graph vividly depicts the execution times, contrasting the performance across the varied hardware configurations.

![Execution Times](img/execution_time.png)

# Discussion

## Performance Analysis

- **CPU Performance:** The application demonstrates effective multi-threading on the CPU, particularly on the personal computer with fewer cores. The use of `std::for_each` with parallel execution policies contributes to this efficiency.
- **GPU Performance:** The relative underperformance on the GPU in the HPC environment suggests potential challenges, possibly related to memory transfer overhead or kernel optimization.

## Parallelization Challenges

- **Scalability:** The application's scalability with increasing CPU cores is a point of consideration, particularly on HPC systems.
- **Optimization:** Tuning the application for optimal performance on GPUs may require addressing specific challenges like memory management and execution efficiency.

# Conclusion

The analysis reveals that the heat simulation application performs comparably on CPU and GPU in a personal computer setup but faces challenges in the HPC environment, particularly in GPU optimization. The findings underscore the importance of tailored parallelization strategies and optimization techniques for different hardware configurations in high-performance computing applications. Further exploration and optimization, especially of the GPU implementation, could lead to significant performance improvements.