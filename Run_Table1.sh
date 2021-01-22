#!/bin/bash
#! What is this job called?
#SBATCH --job-name=R_table
#! What should I call the reports
##SBATCH --output=table_%A_%a
##SBATCH --error=table_%A_%a
#! Which project should be charged:
#SBATCH -A BURGESS-SL3-CPU
#! Which partition/cluster am I using?
#SBATCH -p skylake
#! How many nodes should be allocated? 
#SBATCH --nodes=1
#! How many tasks will there be in total? By default SLURM
#SBATCH --ntasks=1
#! How much memory in MB is required _per node_? 
##SBATCH --mem=5Gb
#! How many cpus per task
#SBATCH --cpus-per-task=8
#! How much wallclock time will be required?
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --array=1

#FILENAME: snptest_run
#AUTHOR : Amy Mason
#PURPOSE: run 5_gen_save_output on range of inputs 
#INPUT: degree, beta1, beta2, par1, par2
#OUTPUT: csv file of associations with snps for each outcome called Table 1

module purge
module load R/4.0.3

Rscript "/rds/project/jmmh2/rds-jmmh2-projects/zz_mr/nlmr/Non-linear MR/5_gen_save_output.R" "1" "2"