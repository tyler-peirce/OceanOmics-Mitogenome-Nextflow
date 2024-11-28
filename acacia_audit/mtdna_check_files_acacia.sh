#!/bin/bash
# Load in the configfile
. ../configfile.txt
mkdir FileCheck
while IFS= read -r line; do
# Set OG and DATE variables (replace with actual values)
    OG=$(echo $line | awk '{print $1}')
    DATE=$(echo $line | awk '{print $4}')
    RUN=$(echo $line | awk '{print $3}')
    TAX=$(cat taxon.txt | grep -w $OG | awk -F'\t' '{print substr($2, 1, 4)}')
    TAX_no=$(cat taxon.txt | grep -w $OG | awk -F'\t' '{print $4}')
    OGdir="pawsey0812:oceanomics-genomes/genomes.v2/$OG"
    #echo $TAX
    # List of expected files and directories
    file_list=(
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.busco_sequences.tar.gz
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.busco_sequences.tar.gz.md5
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.full_table.tsv
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.logs/bbtools_err.log
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.logs/bbtools_out.log
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.logs/busco.log
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.logs/hmmsearch_err.log
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.logs/hmmsearch_out.log
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.logs/metaeuk_run1_err.log
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.logs/metaeuk_run1_out.log
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.logs/metaeuk_run2_err.log
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.logs/metaeuk_run2_out.log
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.logs
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.missing_busco_list.tsv
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.short_summary.json
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.short_summary.json.tsv
        $OGdir/assemblies/genome/busco_$TAX/$OG.ilmn.$DATE.v129mh.fna.busco.$TAX.short_summary.txt
        $OGdir/assemblies/genome/busco_$TAX
        $OGdir/assemblies/genome/bwa/$OG.sorted.bam
        $OGdir/assemblies/genome/bwa/$OG-sn_results.tsv
        $OGdir/assemblies/genome/bwa
        $OGdir/assemblies/genome/depthsizer/$OG.ilmn.$DATE.v129mh.fna.depthsizer.busco.dupcnv.tsv
        $OGdir/assemblies/genome/depthsizer/$OG.ilmn.$DATE.v129mh.fna.depthsizer.fulltable.tsv
        $OGdir/assemblies/genome/depthsizer/$OG.ilmn.$DATE.v129mh.fna.depthsizer.gensize.tdt
        $OGdir/assemblies/genome/depthsizer/$OG.sorted.bam.busco.fastmp.scdepth
        $OGdir/assemblies/genome/depthsizer
        $OGdir/assemblies/genome/NCBI/adaptor/cleaned_sequences/$OG.ilmn.$DATE.v129mh.fa
        $OGdir/assemblies/genome/NCBI/adaptor/cleaned_sequences
        $OGdir/assemblies/genome/NCBI/adaptor/combined.calls.jsonl
        $OGdir/assemblies/genome/NCBI/adaptor/fcs.log
        $OGdir/assemblies/genome/NCBI/adaptor/$OG.ilmn.$DATE.fcs_adaptor.log
        $OGdir/assemblies/genome/NCBI/adaptor/$OG.ilmn.$DATE.fcs_adaptor_report.txt
        $OGdir/assemblies/genome/NCBI/adaptor/logs.jsonl
        $OGdir/assemblies/genome/NCBI/adaptor/$OG.ilmn.$DATE.pipeline_args.yaml
        $OGdir/assemblies/genome/NCBI/adaptor/$OG.ilmn.$DATE.skipped_trims.jsonl
        $OGdir/assemblies/genome/NCBI/adaptor/validate_fasta.txt
        $OGdir/assemblies/genome/NCBI/adaptor
        $OGdir/assemblies/genome/NCBI/$OG.contam.fasta
        $OGdir/assemblies/genome/NCBI/$OG.ilmn.$DATE.contig_count_500bp.txt
        $OGdir/assemblies/genome/NCBI/$OG.ilmn.$DATE.filter_report.txt
        $OGdir/assemblies/genome/NCBI/$OG.ilmn.$DATE.review_scaffolds_1kb.txt
        $OGdir/assemblies/genome/NCBI/$OG.ilmn.$DATE.v129mh.fasta.$TAX_no.fcs_gx_report.txt
        $OGdir/assemblies/genome/NCBI/$OG.ilmn.$DATE.v129mh.fasta.$TAX_no.taxonomy.rpt
        $OGdir/assemblies/genome/NCBI/$OG.ilmn.$DATE.v129mh.rc.fasta
        $OGdir/assemblies/genome/NCBI/$OG.ilmn.$DATE.v129mh.rf.fa
        $OGdir/assemblies/genome/NCBI
        $OGdir/assemblies/genome/$OG.ilmn.$DATE.adaptor-contam.fasta
        $OGdir/assemblies/genome/$OG.ilmn.$DATE.rmadapt.fasta
        $OGdir/assemblies/genome/$OG.ilmn.$DATE.v129mh.fa
        $OGdir/assemblies/genome/$OG.ilmn.$DATE.v129mh.fasta.gz
        $OGdir/assemblies/genome/$OG.ilmn.$DATE.v129mh.fna
        $OGdir/assemblies/genome/tiara/log_$OG.ilmn.$DATE.tiara.txt
        $OGdir/assemblies/genome/tiara/mitochondrion_$OG.ilmn.$DATE.rmadapt.fasta
        $OGdir/assemblies/genome/tiara/$OG.ilmn.$DATE.tiara.contig_removal.txt
        $OGdir/assemblies/genome/tiara/$OG.ilmn.$DATE.tiara.txt
        $OGdir/assemblies/genome/tiara/$OG.ilmn.$DATE.tiara_filter_summary.txt
        $OGdir/assemblies/genome/tiara/plastid_$OG.ilmn.$DATE.rmadapt.fasta
        $OGdir/assemblies/genome/tiara/prokarya_$OG.ilmn.$DATE.rmadapt.fasta
        $OGdir/assemblies/genome/tiara
        $OGdir/assemblies/genome
        $OGdir/assemblies
        $OGdir/fastp/$OG.ilmn.$RUN.fastp.json.tsv
        $OGdir/fastp/$OG.ilmn.$RUN.fastp.json
        $OGdir/fastp/$OG.ilmn.$DATE.fastp.html
        $OGdir/fastp/$OG.ilmn.$DATE.R1.fastq.gz
        $OGdir/fastp/$OG.ilmn.$DATE.R2.fastq.gz
        $OGdir/fastp/fastqc/$OG.ilmn.${RUN}_logs/$OG.ilmn.$RUN.R1_fastqc.html
        $OGdir/fastp/fastqc/$OG.ilmn.${RUN}_logs/$OG.ilmn.$RUN.R1_fastqc.zip
        $OGdir/fastp/fastqc/$OG.ilmn.${RUN}_logs/$OG.ilmn.$RUN.R2_fastqc.html
        $OGdir/fastp/fastqc/$OG.ilmn.${RUN}_logs/$OG.ilmn.$RUN.R2_fastqc.zip
        $OGdir/fastp/fastqc/$OG.ilmn.${RUN}_logs
        $OGdir/fastp/fastqc
        $OGdir/fastp
        $OGdir/kmers/$OG.ilmn.$DATE/$OG.ilmn.$DATE-genomescope_linear_plot.png
        $OGdir/kmers/$OG.ilmn.$DATE/$OG.ilmn.$DATE-genomescope_log_plot.png
        $OGdir/kmers/$OG.ilmn.$DATE/$OG.ilmn.$DATE-genomescope_model.txt
        $OGdir/kmers/$OG.ilmn.$DATE/$OG.ilmn.$DATE-genomescope_progress.txt
        $OGdir/kmers/$OG.ilmn.$DATE/$OG.ilmn.$DATE-genomescope_summary.txt
        $OGdir/kmers/$OG.ilmn.$DATE/$OG.ilmn.$DATE-genomescope_transformed_linear_plot.png
        $OGdir/kmers/$OG.ilmn.$DATE/$OG.ilmn.$DATE-genomescope_transformed_log_plot.png
        $OGdir/kmers/$OG.ilmn.$DATE
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000000.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000000.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000001.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000001.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000010.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000010.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000011.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000011.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000100.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000100.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000101.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000101.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000110.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000110.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000111.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x000111.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001000.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001000.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001001.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001001.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001010.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001010.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001011.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001011.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001100.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001100.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001101.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001101.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001110.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001110.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001111.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x001111.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010000.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010000.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010001.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010001.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010010.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010010.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010011.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010011.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010100.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010100.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010101.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010101.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010110.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010110.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010111.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x010111.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011000.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011000.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011001.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011001.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011010.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011010.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011011.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011011.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011100.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011100.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011101.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011101.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011110.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011110.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011111.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x011111.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100000.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100000.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100001.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100001.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100010.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100010.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100011.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100011.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100100.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100100.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100101.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100101.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100110.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100110.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100111.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x100111.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101000.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101000.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101001.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101001.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101010.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101010.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101011.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101011.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101100.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101100.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101101.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101101.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101110.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101110.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101111.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x101111.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110000.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110000.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110001.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110001.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110010.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110010.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110011.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110011.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110100.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110100.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110101.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110101.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110110.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110110.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110111.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x110111.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111000.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111000.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111001.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111001.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111010.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111010.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111011.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111011.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111100.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111100.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111101.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111101.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111110.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111110.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111111.merylData
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/0x111111.merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl/merylIndex
        $OGdir/kmers/$OG.ilmn.$DATE.meryl
        $OGdir/kmers/$OG.ilmn.$DATE.meryl.tar.gz
        $OGdir/kmers/$OG.ilmn.$DATE.meryl.hist
        $OGdir/kmers/$OG.ilmn.$DATE.v129mh.merqury.completeness.stats
        $OGdir/kmers/$OG.ilmn.$DATE.v129mh.merqury.$OG.ilmn.$DATE.v129mh.fna.qv
        $OGdir/kmers/$OG.ilmn.$DATE.v129mh.merqury.$OG.ilmn.$DATE.v129mh.fna.spectra-cn.fl.png
        $OGdir/kmers/$OG.ilmn.$DATE.v129mh.merqury.$OG.ilmn.$DATE.v129mh.fna.spectra-cn.ln.png
        $OGdir/kmers/$OG.ilmn.$DATE.v129mh.merqury.$OG.ilmn.$DATE.v129mh.fna.spectra-cn.st.png
        $OGdir/kmers/$OG.ilmn.$DATE.v129mh.merqury.qv
        $OGdir/kmers/$OG.ilmn.$DATE.v129mh.merqury.spectra-asm.fl.png
        $OGdir/kmers/$OG.ilmn.$DATE.v129mh.merqury.spectra-asm.ln.png
        $OGdir/kmers/$OG.ilmn.$DATE.v129mh.merqury.spectra-asm.st.png
        $OGdir/kmers
    )

    # Iterate over the list and check if each file or directory exists
    echo -e $OG > FileCheck/$OG.files.tsv
    for file in "${file_list[@]}"; do
        
        if [ -e "$file" ]; then
            echo -e "Found"
        else
            echo -e "Missing"
        fi
    done >> FileCheck/$OG.files.tsv
done



paste 00a_file_list.tsv FileCheck/* > FileCheck.tsv