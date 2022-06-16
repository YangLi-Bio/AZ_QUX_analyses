#!/bin/bash
#SBATCH --job-name=Convert_qsave_to_RDS
#SBATCH --time=2:20:59
#SBATCH --output=Convert_qsave_to_RDS.out
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=500GB
#SBATCH --gpus-per-node=1


set -e


cd /fs/ess/PCON0022/liyang/astrazeneca/QUX/


module load R/4.1.0-gnu9.1
Rscript Convert_qsave_to_rds.R
