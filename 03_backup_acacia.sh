#!/bin/bash

#SBATCH --account=pawsey0812
#SBATCH --job-name=mito_backup
#SBATCH --partition=copy
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00
#SBATCH --export=ALL
#SBATCH --output=%x-%j.out #SBATCH --error=%x-%j.err


rclone copy /scratch/pawsey0964/tpeirce/_MITOGENOMES/mitogenomes/ pawsey0964:oceanomics-mitochondrial-genomes --checksum --progress