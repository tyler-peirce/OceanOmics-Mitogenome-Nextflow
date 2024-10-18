#!/bin/bash --login
#SBATCH --account=pawsey0812
#SBATCH --job-name=fastqc
#SBATCH --partition=work
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=12G

#SBATCH --export=NONE

#-----------------
#Loading the required modules

module load nextflow/23.10.0

. ../configfile.txt

unset SBATCH_EXPORT

nextflow run mitogenomes.nf -c nextflow.config -profile setonix -resume -disable-jobs-cancellation -with-report mitogenomes.html
