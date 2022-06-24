#!/bin/bash
#SBATCH --job-name=Integrate_Inserm_and_Jose
#SBATCH --time=12:20:59
#SBATCH --output=Integrate_Inserm_and_Jose.out
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=300GB
#SBATCH --gpus-per-node=1


set -e


cd /fs/ess/PCON0022/liyang/astrazeneca/QUX/Codes/


module load R/4.1.0-gnu9.1
Rscript 6_Integrate_Inserm_Jose_data.R
