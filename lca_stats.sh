#!/bin/bash
#!/bin/bash
#SBATCH -J lca_sats.slurm
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
touch lca.12s.tsv lca.16s.tsv lca.CO1.tsv

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
                    if [[ $line == lca.12s* ]]; then
                        path="$lca_path/$line"
                        rclone cat "$path" >> lca.12s.tsv
                    elif [[ $line == lca.16s* ]]; then
                        path="$lca_path/$line"
                        rclone cat "$path" >> lca.16s.tsv
                    elif [[ $line == lca.CO1* ]]; then
                        path="$lca_path/$line"
                        rclone cat "$path" >> lca.CO1.tsv
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


