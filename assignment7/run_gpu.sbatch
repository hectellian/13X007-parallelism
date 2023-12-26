#!/bin/sh
#SBATCH --job-name heat_gpu          # this is a parameter to help you sort your job when listing it
#SBATCH --error heat_gpu-error.e%j   # optional. By default a file slurm-{jobid}.out will be created
#SBATCH --output heat_gpu-out.o%j    # optional. By default the error and output files are merged
#SBATCH --ntasks 1                   # number of tasks in your job. One by default
#SBATCH --gpus 2                       # number of cpus for each task. One by default
#SBATCH --partition public-gpu        # the partition to use. By default debug-cpu
#SBATCH --time 15:00                 # maximum run time.

module load foss/2018b
module load NVHPC/21.9
module load CUDA

echo $SLURM_NODELIST

srun make gpu
srun ./build/heat_gpu