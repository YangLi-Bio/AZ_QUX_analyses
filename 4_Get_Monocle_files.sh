#!/bin/bash
#SBATCH --job-name=Get_Monocle_files
#SBATCH --time=24:20:59
#SBATCH --output=Get_Monocle_files.out
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=500GB
#SBATCH --gpus-per-node=1


set -e


cd /fs/ess/PCON0022/liyang/astrazeneca/QUX/


module load R/4.1.0-gnu9.1
Rscript Codes/4_Get_Monocle_files.R
