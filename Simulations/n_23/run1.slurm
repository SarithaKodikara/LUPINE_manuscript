#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Tue Aug 03 2021 11:37:53 GMT+1000 (AEST)

# Partition for the job:
#SBATCH -p mig,cascade,long,physicaltest
##SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# The name of the job:
#SBATCH --job-name="1.BEME_FULL"

# The project ID which this job should run under:
#SBATCH --account="punim0613"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=90G


# Use this email address:
#SBATCH --mail-user=saritha.kodikara@unimelb.edu.au

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=5-1:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)
##Job monitor command to list the resource usage
my-job-stats -a -n -s
echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "$SLURM_NTASKS_PER_NODE tasks per node"
echo "Running on nodes: $SLURM_NODELIST"
echo "------------------------------------------------------------"

# The modules to load:
module load GCC/11.3.0  OpenMPI/4.1.4
module load R/4.2.2 


# The job command(s):
export R_LIBS_USER='~/R/'
mpirun -n 1 R --slave -f run_1.R 

my-job-stats -a -n -s
touch DONE
