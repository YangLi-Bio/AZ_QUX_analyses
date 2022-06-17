#!/bin/bash
#SBATCH --job-name=Run_Generate_annData
#SBATCH --time=5:20:59
#SBATCH --output=Run_Generate_annData.out
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=500GB
#SBATCH --gpus-per-node=1


set -e


cd /fs/ess/PCON0022/liyang/astrazeneca/QUX/


module load python/3.7-2019.10
source activate GLUE_env
python Codes/3_Generate_annData.py
conda deactivate
