#!/bin/bash
#
#SBATCH --mail-user=vhaghani@ucdavis.edu                       # User email to receive updates
#SBATCH --mail-type=ALL                                        # Get an email when the job begins, ends, or if it fails
#SBATCH -p production                                          # Partition, or queue, to assign to
#SBATCH -J test                                                # Name for job
#SBATCH -o test_slurm.j%j.out                                  # File to write STDOUT to
#SBATCH -e test_slurm.j%j.err                                  # File to write error output to
#SBATCH -N 1                                                   # Number of nodes/computers
#SBATCH -n 1                                                   # Number of cores
#SBATCH -c 1                                                   # Eight cores per task
#SBATCH -t 00:05:00                                            # Ask for no more than 5 minutes
#SBATCH --mem=1gb                                              # Ask for no more than 1 GB of memory
#SBATCH --chdir=/share/korflab/home/viki/mechip2               # Directory I want the job to run in

# Source .profile so conda can be used
source /share/korflab/home/viki/profile

# Initialize conda
. /software/anaconda3/4.8.3/lssc0-linux/etc/profile.d/conda.sh

# Activate your desired conda environment
conda activate /home/vhaghani/.conda/envs/mechip2

# Fail on weird errors
set -o nounset
set -o errexit
set -x

# Run the snakemake!
echo Hello World > helloworld.txt

# Print out various information about the job
env | grep SLURM                                               # Print out values of the current jobs SLURM environment variables
  
scontrol show job ${SLURM_JOB_ID}                              # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch

# Note: Run dos2unix {filename} if sbatch DOS line break error occurs