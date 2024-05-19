#!/bin/bash

# Request the number of CPUs and compute nodes needed by your job.
# You must instruct the number of nodes and the number of cores.
#
# Request one processor core:
#   --nodes=1 --ntasks-per-node=1
#
# Request two compute nodes, each with ten processor cores:
#   --nodes=2 --ntasks-per-node=10
#
# Request all twenty processor cores on one compute node:
#   --nodes=1 --ntasks-per-node=20
#
#SBATCH --nodes=1 --ntasks-per-node=1


# Estimate how much memory you expect your job to use (in megabytes).
# Common values:
#     4GB    4096
#     8GB    8192
#   ~16GB   16000
#   ~32GB   32000
#   ~64GB   64000
#
#SBATCH --mem=32000


# Specify how long you expect your job to run. By default SLURM will kill jobs
# that over-run their reservation, so you need to make a realistic estimate.
# PLAY NICE!
#
# Request 1 day:
#   --time=24:00:00
#
# Request 1 hour:
#   --time=1:00:00
#
# Request 15 minutes:
#   --time=15:00
#
#SBATCH --time=1-0:00:00 


# Set the job name
#SBATCH --job-name=hair


# If desired, you may set the filenames for your program's output.
# The %j variable will be replaced with the Job ID # of the job.
#SBATCH --error=hair.%j.output.errors
#SBATCH --output=hair.%j.output.log

# Request GPUs
#SBATCH -p gpu
#SBATCH --gpus-per-node=1 -C 2080ti

# Also, change to the working directory
cd .

################################################################################
# THIS IS WHERE YOU SPECIFY THE JOBS/TASKS TO RUN
################################################################################

ulimit -s hard
#module load uri OpenMPI/4.1.1-GCC-11.2.0 CUDA/11.4.1   
module load openmpi/4.1.3+cuda11.6.2-mpirun  
#module load nvhpc

if [ "x$SLURM_NPROCS" = "x" ]
then
  if [ "x$SLURM_NTASKS_PER_NODE" = "x" ]
  then
    SLURM_NTASKS_PER_NODE=1
  fi
  SLURM_NPROCS=`expr $SLURM_JOB_NUM_NODES \* $SLURM_NTASKS_PER_NODE`
else
  if [ "x$SLURM_NTASKS_PER_NODE" = "x" ]
  then
    SLURM_NTASKS_PER_NODE=`expr $SLURM_NPROCS / $SLURM_JOB_NUM_NODES`
  fi
fi

srun hostname -s | sort -u > /tmp/hosts.$SLURM_JOB_ID
awk "{ print \$0 \" slots=$SLURM_NTASKS_PER_NODE\"; }" /tmp/hosts.$SLURM_JOB_ID >/tmp/tmp.$SLURM_JOB_ID
mv /tmp/tmp.$SLURM_JOB_ID /tmp/hosts.$SLURM_JOB_ID
cat /tmp/hosts.$SLURM_JOB_ID 
echo $SLURM_NPROCS

#mpirun -hostfile /tmp/hosts.$SLURM_JOB_ID --mca plm_slurm_args "--mpi=none -G 0 --gpus-per-node=0" -np $SLURM_NPROCS  ./runme.bin
#mpirun -hostfile /tmp/hosts.$SLURM_JOB_ID --mca btl_openib_allow_ib "true" -np $SLURM_NPROCS  ./runme.bin
#date
mpirun -hostfile /tmp/hosts.$SLURM_JOB_ID -np $SLURM_NPROCS  ./runme.bin
#date

