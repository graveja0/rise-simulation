#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH -t 0-0:30:00
#SBATCH -o rise_1_%A_%a.out
#SBATCH -e rise_1_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=zilu.zhou@vanderbilt.edu

module load GCC/5.4.0-2.26  OpenMPI/1.10.3 R/3.3.3-X11-20160819
R --version

echo "SLURM_JOBID: " $SLURM_JOBID
Rscript test0.R 
