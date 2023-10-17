# 13X007 - Parallelism

Welcome to the repository of Parallelism Projects in C++. This collection contains multiple projects that showcase various techniques for achieving parallelism in C++.

## Table of Contents

- [13X007 - Parallelism](#13x007---parallelism)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
  - [Projects](#projects)
    - [Hello World](#hello-world)
    - [Communication](#communication)
  - [How to Run on Baobab](#how-to-run-on-baobab)
  - [License](#license)

## Introduction

Parallelism is the concept of doing multiple things at the same time, and this repository aims to explore this concept within the context of C++ programming. Examples range from simple multi-threaded applications to more complex projects using libraries and frameworks such as OpenMP, Boost.Asio, and more.

## Prerequisites

- OpenMPI
- MPI++ Compiler
- Make 

## Installation

Clone the repository:

```bash
git clone https://github.com/hectellian/13X007-parallelism.git
```

Navigate into the directory:

```bash
cd 13X007-parallelism
```

Compile the projects using Make:

```bash
make
```

## Projects

### Hello World

**Technologies Used:**
- OpenMPI

This project demonstrates basic hello world parallelism.

### Communication

This project showcases the use of MPI to communicate between processes.

**Technologies Used:**
- OpenMPI

## How to Run on Baobab

Baobab is a high-performance computing (HPC) cluster. To run your projects on Baobab, follow these steps:

1. **SSH into Baobab:**

    ```bash
    ssh username@login2.baobab.hpc.unige.ch
    ```

2. **Transfer Files:**

    Use `scp` or `rsync` to transfer your project files to Baobab.

    ```bash
    scp -r 13X007-parallelism username@login2.baobab.hpc.unige.ch:/path/to/destination
    ```

> [!NOTE]
> You can also directly clone the repository on Baobab.

3. **Load Modules:**

    Load the required modules using the `module load` command in your `.bashrc` file.

    ```bash
    # all your other bashrc commands

    module load CUDA
    module load foss
    ```

> - CUDA: This module is required for GPU programming.
> - foss: This module is required for OpenMPI, it contains GCC and other libraries.

4. **Compile and Run:**

    Follow the regular [Installation](#installation) steps to compile and then run your project by sending it to the queue with an `sbatch` script file that contains the following:

    ```bash
    #!/bin/sh
    #SBATCH --job-name jobname            # this is a parameter to help you sort your job when listing it
    #SBATCH --error jobname-error.e%j     # optional. By default a file slurm-{jobid}.out will be created
    #SBATCH --output jobname-out.o%j      # optional. By default the error and output files are merged
    #SBATCH --ntasks 1                    # number of tasks in your job. One by default
    #SBATCH --cpus-per-task 1             # number of cpus for each task. One by default
    #SBATCH --partition debug-cpu         # the partition to use. By default debug-cpu
    #SBATCH --time 15:00                  # maximum run time.

    module load CUDA
    module load foss                      # load a specific software using module, for example Python

    echo $SLURM_NODELIST

    srun --mpi=pmi2 ./your_program
    ```

    For example, to run the [Communication](#communication) project:

    ```bash
    #!/bin/sh
    #SBATCH --job-name broadcast          # this is a parameter to help you sort your job when listing it
    #SBATCH --error broadcast-error.e%j   # optional. By default a file slurm-{jobid}.out will be created
    #SBATCH --output broadcast-out.o%j    # optional. By default the error and output files are merged
    #SBATCH --ntasks 16                   # number of tasks in your job. One by default
    #SBATCH --cpus-per-task 1             # number of cpus for each task. One by default
    #SBATCH --partition debug-cpu         # the partition to use. By default debug-cpu
    #SBATCH --time 15:00                  # maximum run time.

    module load CUDA
    module load foss                      # load a specific software using module, for example Python

    echo $SLURM_NODELIST

    srun --mpi=pmi2 ./build/communication broadcast
    srun --mpi=pmi2 ./build/communication ring
    srun --mpi=pmi2 ./build/communication hypercube
    ```


    Compile:

    ```bash
    make
    ```

    Run:

    ```bash
    sbatch run.sh
    ```

    You should obtain 2 files: `jobname-error.e{jobid}` and `jobname-out.o{jobid}`. The first one contains the error output and the second one contains the standard output.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

---

For more information or questions, feel free to contact the repository owner.