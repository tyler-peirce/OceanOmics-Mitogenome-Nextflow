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


#rclone move /scratch/pawsey0812/tpeirce/MITOGENOMES/ilmn/ pawsey0812:oceanomics-mitogenomes --checksum --progress
rclone copy /scratch/pawsey0812/tpeirce/MITOGENOMES/ilmn/ pawsey0812:oceanomics-mitogenomes --checksum --progress