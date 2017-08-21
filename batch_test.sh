#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=225G
#SBATCH -t 0-48:00:00
#SBATCH -o rise_1_%A_%a.out
#SBATCH -e rise_1_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=john.graves@vanderbilt.edu
module load Intel/2016.3.210-GCC-5.4.0-2.26  IntelMPI/5.1.3.181
module load GCC/5.4.0-2.26  OpenMPI/1.10.3 R/3.4.0-X11-20160819
R --version

echo "SLURM_JOBID: " $SLURM_JOBID
Rscript test0.R 
