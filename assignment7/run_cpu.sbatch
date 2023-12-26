#!/bin/sh
#SBATCH --job-name heat_cpu          # this is a parameter to help you sort your job when listing it
#SBATCH --error heat_cpu-error.e%j   # optional. By default a file slurm-{jobid}.out will be created
#SBATCH --output heat_cpu-out.o%j    # optional. By default the error and output files are merged
#SBATCH --ntasks 1                   # number of tasks in your job. One by default
#SBATCH --cpus-per-task 32           # number of cpus for each task. One by default
#SBATCH --partition debug-cpu        # the partition to use. By default debug-cpu
#SBATCH --time 15:00                 # maximum run time.

module load foss/2020b
module load tbb

echo $SLURM_NODELIST

make cpu
srun ./build/heat_cpu