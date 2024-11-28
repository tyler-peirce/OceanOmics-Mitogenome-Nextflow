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

# this code pulls all of the stats for the mitogenomes on Acaia

# Create output file
echo -e "OG_id\ttech\tdate\tcode\tstats\tlength\tlength_emma\tcds\ttrna\trrna\t12s\t16s\tATP6\tATP8\tCOX1\tCOX2\tCOX3\tCYTB\tND1\tND2\tND3\tND4\tND4L\tND5\tND6" >> mtdnastat."$(date '+%y%m%d')".tsv

#create list of all the samples that have mitogenomes
file_list=$(rclone lsf pawsey0812:oceanomics-mitogenomes)

for og in $file_list; do
    echo "og = $og"

    #Get the contents of each directory
    name=$(rclone lsf pawsey0812:oceanomics-mitogenomes/"$og"/)
    echo "name = $name"

    #check if there are any files and if there is more then 1 in the directory
        if [ $(echo "$name" | wc -l) -ge 1 ]; then
            echo "$name" | while read -r line; do
            echo "line = $line"

            # Pull out the stats of the mitogenome for data assembled through GetOrganelle
            if [[ $line != *hifi* ]]; then
                # Details on whether the mitogenome is circular
                path=pawsey0812:oceanomics-mitogenomes/$og"$line"
                stats=$(rclone cat "$path"mtdna/get_org.log.txt | grep "Result status of animal_mt: " | tail -n 1 | awk -F 'animal_mt: ' '{print $NF}')
                
                # Count the length of the sequence
                length=$(rclone cat "$path"mtdna/${line%/}.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                length_emma=$(rclone cat "$path"emma/${line%/}.rotated.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ###need to change this as it doesnt do it for the scaffolds###
                gff=""$path"emma/${line%/}.rotated.fasta.gff"
                cds=$(rclone cat $gff | awk '$3 == "CDS" { count++ } END { print count }')
                trna=$(rclone cat $gff | awk '$3 == "tRNA" { count++ } END { print count }')
                rrna=$(rclone cat $gff | awk '$3 == "rRNA" { count++ } END { print count }')

                codes=""$path"emma/cds/${line%/}.rotated"
                s12=$(rclone cat "$codes".12srna.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                s16=$(rclone cat "$codes".16srna.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ATP6=$(rclone cat "$codes".ATP6.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ATP8=$(rclone cat "$codes".ATP8.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                COX1=$(rclone cat "$codes".COX1.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                COX2=$(rclone cat "$codes".COX2.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                COX3=$(rclone cat "$codes".COX3.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                CYTB=$(rclone cat "$codes".CYTB.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND1=$(rclone cat "$codes".ND1.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND2=$(rclone cat "$codes".ND2.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND3=$(rclone cat "$codes".ND3.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND4=$(rclone cat "$codes".ND4.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND4L=$(rclone cat "$codes".ND4L.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND5=$(rclone cat "$codes".ND5.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND6=$(rclone cat "$codes".ND6.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                
                # Split up the assmbly name for the database
                key=$(echo ${line%/} | tr '.' '\t')
                echo -e "$key\t$stats\t$length\t$length_emma\t$cds\t$trna\t$rrna\t$s12\t$s16\t$ATP6\t$ATP8\t$COX1\t$COX2\t$COX3\t$CYTB\t$ND1\t$ND2\t$ND3\t$ND4\t$ND4L\t$ND5\t$ND6" >> mtdnastat."$(date '+%y%m%d')".tsv
            else
            
            
            
            # Pull out the stats of the mitogenome for hifi data assembled by mitohifi
                # Details on whether the mitogenome is circular
                path=pawsey0812:oceanomics-mitogenomes/$og"$line"
                stats=$(rclone cat "$path"mtdna/contigs_stats.tsv | grep "final_mitogenome" | awk -F "\t" {'print $6'})
                    if [ "$stats" = True ]; then
                        stats="circular genome"
                    else
                        stats="not circular"
                    fi
                # Count the length of the sequence
                length=$(rclone cat "$path"mtdna/${line%/}.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                length_emma=$(rclone cat "$path"emma/${line%/}.rotated.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                
                gff=""$path"emma/${line%/}.rotated.fasta.gff"
                cds=$(rclone cat $gff | awk '$3 == "CDS" { count++ } END { print count }')
                trna=$(rclone cat $gff | awk '$3 == "tRNA" { count++ } END { print count }')
                rrna=$(rclone cat $gff | awk '$3 == "rRNA" { count++ } END { print count }')

                codes=""$path"emma/cds/${line%/}.rotated"
                s12=$(rclone cat "$codes".12srna.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                s16=$(rclone cat "$codes".16srna.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ATP6=$(rclone cat "$codes".ATP6.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ATP8=$(rclone cat "$codes".ATP8.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                COX1=$(rclone cat "$codes".COX1.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                COX2=$(rclone cat "$codes".COX2.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                COX3=$(rclone cat "$codes".COX3.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                CYTB=$(rclone cat "$codes".CYTB.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND1=$(rclone cat "$codes".ND1.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND2=$(rclone cat "$codes".ND2.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND3=$(rclone cat "$codes".ND3.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND4=$(rclone cat "$codes".ND4.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND4L=$(rclone cat "$codes".ND4L.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND5=$(rclone cat "$codes".ND5.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                ND6=$(rclone cat "$codes".ND6.fasta | grep -v "^>" | tr -d '\n' | awk '{print length}')
                
                # Split up the assmbly name for the database
                key=$(echo ${line%/} | tr '.' '\t')
                echo -e "$key\t$stats\t$length\t$length_emma\t$cds\t$trna\t$rrna\t$s12\t$s16\t$ATP6\t$ATP8\t$COX1\t$COX2\t$COX3\t$CYTB\t$ND1\t$ND2\t$ND3\t$ND4\t$ND4L\t$ND5\t$ND6" >> mtdnastat."$(date '+%y%m%d')".tsv
        
                    
            fi
        done
    else
        echo "no files for $name"
        fi

done

