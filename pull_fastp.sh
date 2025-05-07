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
    echo "Processing: $og"

    dest_dir="/scratch/pawsey0964/tpeirce/_MITOGENOMES/mitogenomes/$og"
    src_path="pawsey0964:oceanomics-mitogenomes/$og"

    mkdir -p "$dest_dir"

    # Check if source path exists using rclone ls (silent if missing)
    if rclone ls "$src_path" >/dev/null 2>&1; then
        rclone copy "$src_path" "$dest_dir" --checksum --progress
    else
        echo "❌ ERROR: Source not found for OG: $og" >&2
    fi

done < "$list"
