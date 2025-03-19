#!/bin/bash

rundir=/scratch/pawsey0812/tpeirce/MITOGENOMES/ilmn
rundate=241101
file="/scratch/pawsey0812/tpeirce/OceanOmics-Mitogenome-Nextflow/nomspecies.txt"
echo -e "OG\tSample\tNomSpecies\tCO1\t16s\t12s" > lcas.$rundate.tsv

for O in $rundir/*; do
    #echo "O = $O"
    samp=$(basename $O/*)
    #echo "sampe = $samp"
    OG=$(basename $O)
    #echo "OG = $OG"
    lca=$O/OG*/lca
    #echo "lca=$lca"
    #ls $lca
    CO1=$(cat $lca/lca.CO1.* | awk -F '\t' '{print $3}')
    sixteen=$(cat $lca/lca.16s.* | awk -F '\t' '{print $3}')
    twelve=$(cat $lca/lca.12s.* | awk -F '\t' '{print $3}')
    species=$(cat $file | grep -w -m 1 "$OG" | awk -F '\t' '{print $2}')
    
    #echo "$species"
    #echo "CO1 = $CO1"
    #echo "sixteen = $sixteen"
    #echo "twelve = $twelve"

    echo -e "$OG\t$samp\t$species\t$CO1\t$sixteen\t$twelve" >> lcas.$rundate.tsv
done