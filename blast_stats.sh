#!/bin/bash
#!/bin/bash
#SBATCH -J stats_acacia.slurm
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --partition=work
#SBATCH --clusters=setonix
#SBATCH --account=pawsey0812
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tyler.peirce@uwa.edu.au
#SBATCH --output=%x-%j.out  #SBATCH --error=%x-%j.err


#create list of all the samples that have mitogenomes
file_list=$(rclone lsf pawsey0812:oceanomics-mitogenomes)
touch blast.12s.tsv blast.16s.tsv blast.CO1.tsv

for og in $file_list; do
    echo "Processing OG: $og"

    # Get the contents of each directory
    name=$(rclone lsf "pawsey0812:oceanomics-mitogenomes/$og")
    echo "Contents of $og: $name"

    # Check if there are files in the directory
    if [ -n "$name" ]; then
        echo "$name" | while read -r line1; do
            echo "Processing directory: $line1"

            # Path to the LCA directory
            lca_path="pawsey0812:oceanomics-mitogenomes/${og}${line1}lca"
            lcas=$(rclone lsf "$lca_path")
            echo "LCA directory contents: $lcas"

            # Check if the LCA directory has files
            if [ -n "$lcas" ]; then
                
                echo "$lcas" | while read -r line; do
                    echo "Processing LCA file: $line"
                    # Pull out just the LCA files and append to appropriate file
                    if [[ $line == blast.12s* ]]; then
                        path="$lca_path/$line"
                        rclone cat "$path" >> blast.12s.tsv
                    elif [[ $line == blast.16s* ]]; then
                        path="$lca_path/$line"
                        rclone cat "$path" >> blast.16s.tsv
                    elif [[ $line == blast.CO1* ]]; then
                        path="$lca_path/$line"
                        rclone cat "$path" >> blast.CO1.tsv
                    else
                        echo "Unknown file pattern: $line"
                    fi
                done
            else
                echo "No LCA files found in: ${og}${line1}lca"
            fi
        done
    else
        echo "No files found in OG: $og"
    fi
done


