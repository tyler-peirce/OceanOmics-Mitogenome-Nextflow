#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


//_________________________________________________________________________________________________________
// |||| Pipeline input parameters ||||
//_________________________________________________________________________________________________________

params.rundir = "/scratch/pawsey0812/tpeirce/MITOGENOMES/OG115" // Put the path to the parent directory for the OG dirs to follow file path for params.fastq
params.mitodir = "/scratch/pawsey0812/tpeirce/MITOGENOMES/ilmn" // The output parent directory
params.fastq="$params.rundir/*/*.{R1,R2}.fastq.gz"
//params.fastq="$params.rundir/OG*/fastp/*.{R1,R2}.fastq.gz" // This is connected to the Draft Genome pipeline output dir
//params.getorg_db = "/scratch/pawsey0812/tpeirce/.GetOrganelle" // this is now redundant because of the new process, delete if it works fine.
params.organelle_type = "animal_mt"
params.lca = "/scratch/pawsey0812/pbayer/OceanGenomes.CuratedNT.NBDLTranche1.CuratedBOLD.fasta" // The curated OG database, curated by Philipp
params.taxdb = "/scratch/pawsey0812/tpeirce/MITOGENOMES/blast_database/*" // The directory that you have "wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz" and then "tar xzvf taxdb.tar.gz"
params.taxonkit="/scratch/pawsey0812/tpeirce/MITOGENOMES/blast_database/" // The direcotry that you have "wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz" and then "tar xzvf taxdump.tar.gz"
//params.EMMA = "/scratch/pawsey0812/tpeirce/MITOGENOMES/ACACIA/OG*/OG*/mtdna/OG*.fasta" // For test running Emma process
//params.annotation = "${params.mitodir}/OG*/*/mtdna/*1.1*.fasta" // This params is to be used when you dont need to do the assembly
//_________________________________________________________________________________________________________
// |||| Processes ||||
//_________________________________________________________________________________________________________

    //_________________________________________________________________________________________________________
    // MITOHIFI_FINDMITOREFERENCE
    //_________________________________________________________________________________________________________

    process MITOHIFI_FINDMITOREFERENCE {
        tag "$species"
        label 'process_single'

        // Docker image available at the project github repository
        container 'ghcr.io/marcelauliano/mitohifi:master'

        input:
        tuple val(meta), val(species)

        output:
        tuple val(meta), path("*.fasta"), emit: fasta
        tuple val(meta), path("*.gb")   , emit: gb
        path "versions.yml"             , emit: versions

        when:
        task.ext.when == null || task.ext.when

        script:
        def args = task.ext.args ?: ''
        """
        findMitoReference.py \\
            --species "$species" \\
            --outfolder . \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mitohifi: \$( mitohifi.py -v | sed 's/.* //' )
        END_VERSIONS
        """

        stub:
        """
        touch accession.fasta
        touch accession.gb

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mitohifi: \$( mitohifi.py -v | sed 's/.* //' )
        END_VERSIONS
        """
    }
    
    
    //_________________________________________________________________________________________________________
    // MITOHIFI_MITOHIFI - Extracting the mitogenome from fasta files
    //_________________________________________________________________________________________________________

    process MITOHIFI_MITOHIFI {
        tag "$meta.id"
        label 'process_high'


        // Docker image available at the project github repository
        container 'ghcr.io/marcelauliano/mitohifi:master'

        input:
        tuple val(meta), path(input)
        path ref_fa
        path ref_gb
        val input_mode
        val mito_code

        output:
        tuple val(meta), path("*fasta")                          , emit: fasta
        tuple val(meta), path("*contigs_stats.tsv")              , emit: stats
        tuple val(meta), path("*gb")                             , emit: gb, optional: true
        tuple val(meta), path("*gff")                            , emit: gff, optional: true
        tuple val(meta), path("*all_potential_contigs.fa")       , emit: all_potential_contigs, optional: true
        tuple val(meta), path("*contigs_annotations.png")        , emit: contigs_annotations, optional: true
        tuple val(meta), path("*contigs_circularization")        , emit: contigs_circularization, optional: true
        tuple val(meta), path("*contigs_filtering")              , emit: contigs_filtering, optional: true
        tuple val(meta), path("*coverage_mapping")               , emit: coverage_mapping, optional: true
        tuple val(meta), path("*coverage_plot.png")              , emit: coverage_plot, optional: true
        tuple val(meta), path("*final_mitogenome.annotation.png"), emit: final_mitogenome_annotation, optional: true
        tuple val(meta), path("*final_mitogenome_choice")        , emit: final_mitogenome_choice, optional: true
        tuple val(meta), path("*final_mitogenome.coverage.png")  , emit: final_mitogenome_coverage, optional: true
        tuple val(meta), path("*potential_contigs")              , emit: potential_contigs, optional: true
        tuple val(meta), path("*reads_mapping_and_assembly")     , emit: reads_mapping_and_assembly, optional: true
        tuple val(meta), path("*shared_genes.tsv")               , emit: shared_genes, optional: true
        path  "versions.yml"                                     , emit: versions

        when:
        task.ext.when == null || task.ext.when

        script:
        // Exit if running this module with -profile conda / -profile mamba
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            error "MitoHiFi module does not support Conda. Please use Docker / Singularity instead."
        }

        def args = task.ext.args ?: ''
        if (! ["c", "r"].contains(input_mode)) {
            error "r for reads or c for contigs must be specified"
        }
        """
        mitohifi.py -${input_mode} ${input} \\
            -f ${ref_fa} \\
            -g ${ref_gb} \\
            -o ${mito_code} \\
            -t $task.cpus ${args}
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mitohifi: \$( mitohifi.py --version 2>&1 | head -n1 | sed 's/^.*MitoHiFi //; s/ .*\$//' )
        END_VERSIONS
        """
    }


    //_________________________________________________________________________________________________________
    // EMMA - Annotation of the mitogenomes
    //_________________________________________________________________________________________________________

    process EMMA { 
        tag "EMMA annotation on $og_num"

        publishDir "${params.mitodir}/${og_num}/${prefix}", mode:'copy'
        
        input:

            tuple val(og_num), val(prefix), path(fasta)
            

        output:
            tuple val(og_num), val(prefix), path("emma_prefix.txt"), path("emma"), emit: blast
            path "emma/versions_emma.yml"
          
            
        script:
            def args   = task.ext.args ?: '--rotate MT-TF'
            
            """ 
            version=\$(julia -v | sed 's/^julia version //g' | sed 's/\\.//g')
            emma_version=\$(cat /opt/Emma/Project.toml | grep version | sed -E 's/version = "([0-9]+)\\.([0-9]+)\\.([0-9]+)"/\\1\\2\\3/')
                     
            emma_prefix="${prefix}.emma"\${emma_version}""
            
            mkdir -p tempdir cds emma
            
            /opt/julia-1.10.5/bin/julia \\
                /opt/Emma/src/command.jl \\
                $args \\
                --fa emma/\${emma_prefix}.fa \\
                --gff emma/\${emma_prefix}.fa.gff \\
                --tbl emma/\${emma_prefix}.tbl \\
                --svg emma/\${emma_prefix}.svg \\
                --tempdir tempdir/ \\
                --loglevel debug \\
                ${fasta} 
            
            /opt/julia-1.10.5/bin/julia /opt/extract_proteins.jl \\
                emma/ \\
                emma/
            
            mv cds emma/  

            # Add in the sample to the file names in cds
            for file in emma/cds/*; do
                new_file="\${file%.fa}.\${emma_prefix}.fa"     
                mv "\$file" "\$new_file"
            done
                       
            cat <<-END_VERSIONS > emma/versions_emma.yml
            "${task.process}":
                julia: \$(julia -v | sed 's/^julia version //g' )
                emma: \$(cat /opt/Emma/Project.toml | grep version)
            END_VERSIONS

            # Return the calculated emma_prefix so it can be passed to the next process
            echo \$emma_prefix > emma_prefix.txt

            """
            
    }

    //_________________________________________________________________________________________________________
    // BLAST - this process blasts the COX1, 16s and 12s genes against the OceanGenomes curated databse
    //_________________________________________________________________________________________________________

    process BLAST { 
        tag "BLAST on $og_num"

        publishDir "${params.mitodir}/${og_num}/${prefix}", mode:'copy'

        input:
            tuple val(og_num), val(prefix), val(emma_prefix), path(emma)
            path "*"
        

        output:
            tuple val(og_num), val(prefix), val(emma_prefix), path("lca"), emit: lca
            path "lca/versions_BLAST.yml"

        script:
            """ 
            version=\$(blastn -version)
            cds=emma/cds
            mkdir -p lca

            blastn \\
                -query \${cds}/MT-CO1.*.fa \\
                -db ${params.lca}  \\
                -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \\
                -evalue 1 \\
                -task blastn  \\
                -num_alignments 999999 \\
                > lca/blast.CO1.${emma_prefix}.tsv
            
            blastn \\
                -query \${cds}/MT-RNR1.*.fa \\
                -db ${params.lca}  \\
                -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \\
                -evalue 1 \\
                -task blastn  \\
                -num_alignments 999999 \\
                > lca/blast.12s.${emma_prefix}.tsv

            blastn \\
                -query \${cds}/MT-RNR2.*.fa \\
                -db ${params.lca}  \\
                -outfmt "6 qseqid sseqid staxids sscinames scomnames sskingdoms pident length qlen slen mismatch gapopen gaps qstart qend sstart send stitle evalue bitscore qcovs qcovhsp" \\
                -evalue 1 \\
                -task blastn  \\
                -num_alignments 999999 \\
                > lca/blast.16s.${emma_prefix}.tsv
            
            wait

            awk -F '\t' '{if ((\$15 - \$14 > 200) && (\$7 > 98)) print}' lca/blast.CO1.${emma_prefix}.tsv > lca/blast.CO1.${emma_prefix}.filtered.tsv
            awk -F '\t' '{if ((\$15 - \$14 > 200) && (\$7 > 98)) print}' lca/blast.16s.${emma_prefix}.tsv > lca/blast.16s.${emma_prefix}.filtered.tsv
            awk -F '\t' '{if ((\$15 - \$14 > 200) && (\$7 > 98)) print}' lca/blast.12s.${emma_prefix}.tsv > lca/blast.12s.${emma_prefix}.filtered.tsv


            
            cat <<-END_VERSIONS > lca/versions_BLAST.yml
            "${task.process}":
                \$(blastn -version)
            END_VERSIONS
            """
    }

    //_________________________________________________________________________________________________________
    // LCA - this process takes the filtered blast results and calculates the LCA
    //_________________________________________________________________________________________________________

    process LCA { 
        tag "LCA on $og_num"

        publishDir "${params.mitodir}/${og_num}/${prefix}", mode:'copy'

        input:
            tuple val(og_num), val(prefix), val(emma_prefix), path(lca)
        

        output:        
            path("lca")
            

        script:
            """ 
            version=\$(python -V | sed 's/Python //g')
            
            env TAXONKIT_DB=$params.taxonkit

            python /computeLCA.py \\
                lca/blast.CO1*filtered.tsv \\
                > lca/lca.CO1.${emma_prefix}.tsv
            
            python /computeLCA.py \\
                lca/blast.16s*filtered.tsv \\
                > lca/lca.16s.${emma_prefix}.tsv
            
            python /computeLCA.py \\
                lca/blast.12s*filtered.tsv \\
                > lca/lca.12s.${emma_prefix}.tsv
            
            sed -i "s/\$/\\t\$(date +%y%m%d)/" lca/lca.*.tsv
            wait
            sed -i "s/\$/\\tCO1/" lca/lca.CO1s*.tsv
            sed -i "s/\$/\\t16s/" lca/lca.16s*.tsv
            sed -i "s/\$/\\t12s/" lca/lca.12s*.tsv

            cat <<-END_VERSIONS > lca/versions_LCA.yml
            "${task.process}":
                Python: \$(python -V | sed 's/Python //g')
                TaxonKit: \$(/taxonkit version | sed 's/taxonkit //g')
            END_VERSIONS
            """
    }


//_________________________________________________________________________________________________________
// |||| Workflow ||||
//_________________________________________________________________________________________________________

workflow {

    //_________________________________________________________________________________________________________________
    // This is the section for running the whole pipeline, starting with using get organelle to assemble the mitogenome.
    // Comment this section out to the solid lines if you just want to run the annotation and LCA.
    //_________________________________________________________________________________________________________________

        organelle_ch = params.organelle_type

        GETORGANELLE_CONFIG(organelle_ch)
                
        read_pairs_ch = Channel
            .fromFilePairs(params.fastq, checkIfExists: true)
            .map { pair ->
                def tokens = pair[0].tokenize('.')
                def og_num = tokens[0]              // Get the first token (og_num)
                def sample = tokens.take(3).join('.') // Join the first 3 tokens to form the sample name
                return tuple(og_num, sample, pair[1]) 
            }
        read_pairs_ch.subscribe { item -> println "read_pairs_ch: $item"}
        

        GETORGANELLE_FROMREADS(read_pairs_ch, GETORGANELLE_CONFIG.out.db)
        
        prefix_ch = GETORGANELLE_FROMREADS.out.fasta  // Change to Channel.fromPath(params.annotation) if you already have the assemblies and "//" out everything above this channel in the workflow
            .map {
                // Use the baseName method to extract the filename without the directo
                def fileName = it.getFileName().toString()
                // Extract the first token (e.g., "OG33")
                def og_num = fileName.tokenize('.')[0] // the [0] takes the first element of the list thats created and provides just the value
                // Extract the first four tokens and join them (e.g., "OG33.ilmn.240716.getorg1770")
                def prefix = fileName.tokenize('.').take(4).join('.')
                return tuple(og_num, prefix, it) 
            }
        prefix_ch.subscribe { item -> println "prefix_ch: $item"}
        
        EMMA(prefix_ch)

    //_________________________________________________________________________________________________________________
    // Comment out to here to just run the annotation and LCA.
    // The next section needs to be un commented to run the annotation and LCA.
    // The next section needs to be commented out to run the whole pipeline. Until the next solid line
    //_________________________________________________________________________________________________________________

        //Comment out everything above to just run the nextflow from EMMA and unclomment out this section
    //    emma_ch = Channel
    //        .fromPath(params.EMMA, checkIfExists: true)
    //        .map {
                // Use the baseName method to extract the filename without the directo
    //            def fileName = it.getFileName().toString()
                // Extract the first token (e.g., "OG33")
    //            def og_num = fileName.tokenize('.')[0] // the [0] takes the first element of the list thats created and provides just the value
                // Extract the first four tokens and join them (e.g., "OG33.ilmn.240716.getorg1770")
    //            def prefix = fileName.tokenize('.').take(4).join('.')
    //            return tuple(og_num, prefix, it) 
    //        }   
        //emma_ch.subscribe { item -> println "emma_ch: $item"}  // Uncomment if you want to check the channel output
        
    //    EMMA(emma_ch)

        //EMMA.out.blast.subscribe { item -> println "Output from EMMA: $item" } // Uncomment if you want to check the channel output
    
    //_________________________________________________________________________________________________________________
    // This remaining section remains uncommented for either workflow.
    //_________________________________________________________________________________________________________________
    

    emma_prefix_ch = EMMA.out.blast
        .map { 
            og_num, prefix, emma_prefix, emma -> 
            def new_emma_prefix = file(emma_prefix).text.trim()
            return tuple(og_num, prefix, new_emma_prefix, emma)
        }
    //emma_prefix_ch.subscribe { item -> println "emma_prefix_ch: $item" } // Uncomment if you want to check the channel output
    

    taxdb_ch = Channel
        .fromPath(params.taxdb, checkIfExists: true)
        
    BLAST(emma_prefix_ch, taxdb_ch.collect())
    LCA(BLAST.out.lca)
    

}