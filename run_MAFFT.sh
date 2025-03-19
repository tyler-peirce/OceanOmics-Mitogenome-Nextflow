#!/bin/bash

# Directories
ALIGNED_PROTEINS_DIR="../aligned_proteins"
#SINGULARITY_PATH="/path/to/singularity"  # Update this with the correct path
MAFFT_CONTAINER="$SING/mafft:7.515.sif"

# Loop through all aligned .fa files
#for aligned_file in "$ALIGNED_PROTEINS_DIR"/*.aligned.fa; do # for running on multiple files in a dir
for aligned_file in CO3.cat.fa; do # for running on a single file
    if [[ -f "$aligned_file" ]]; then
        input_path="$aligned_file"
        output_path="${input_path%.fa}.out.fa"
        
        # Run MAFFT with Singularity
        echo "Running: singularity run $MAFFT_CONTAINER mafft --auto $input_path > $output_path"
        singularity run $MAFFT_CONTAINER mafft --auto "$input_path" > "$output_path"
        
        echo "✅ MAFFT aligned $aligned_file -> $output_path"
    fi
done

echo "🚀 All files processed and aligned with MAFFT!"
