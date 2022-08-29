#!/bin/bash
#SBATCH -p biostat
##SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=sim-03-assoc-%a
#SBATCH --output=sim-03-assoc-%a.out
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=1
##SBATCH --mail-user=alejandro.ochoa@duke.edu
##SBATCH --mail-type=END,FAIL

module load R/4.0.0
module load Plink/2.00a3LM

name='tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01'
rep=$SLURM_ARRAY_TASK_ID
# NOTE: use a single thread per run!
time Rscript sim-03-assoc.R --bfile $name -p 10 -r $rep -t 1

module unload R/4.0.0
module unload Plink/2.00a3LM
