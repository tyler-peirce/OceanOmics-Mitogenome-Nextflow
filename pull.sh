#!/bin/bash --login
#SBATCH --account=pawsey0812
#SBATCH --job-name=aws-raw-backup
#SBATCH --partition=work
#SBATCH --mem=15GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --export=NONE

. ~/.bashrc
#this script is to pull the fastp files if pawsey has deleted them

list=list.txt


while IFS= read -r line; do
    og=$(echo "$line")
    echo "og = $og"
    
    mkdir -p /scratch/pawsey0812/tpeirce/MITOGENOMES/TREND/$og
    rclone copy pawsey0812:oceanomics-mitogenomes/$og /scratch/pawsey0812/tpeirce/MITOGENOMES/TREND/$og --checksum --progress
    
done < "$list"
