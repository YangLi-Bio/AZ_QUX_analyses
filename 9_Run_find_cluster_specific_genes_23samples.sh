#!/bin/bash
#SBATCH --job-name=Find_cluster_specific_genes_23samples
#SBATCH --time=22:20:59
#SBATCH --output=Find_cluster_specific_genes_23samples.out
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=300GB
#SBATCH --gpus-per-node=1


set -e


cd /fs/ess/PCON0022/liyang/astrazeneca/QUX/Codes/


module load R/4.1.0-gnu9.1
Rscript 9_Find_cluster_specific_genes_23samples.R
