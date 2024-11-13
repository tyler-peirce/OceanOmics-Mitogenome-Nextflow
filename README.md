# Mitogenome Nextflow

This is the Mitogenome Nextflow for assembling, annotating and determining the lowest common ancestor for species ID, using ilumina short read data. There are four steps to this pipeline, Getorganelle (assembly), EMMA (annotation), Blast and LCA.

This pipeline can be run by submitting the 01_slurm.sh, or run in the workflow node terminal in a TMUX session.

Make sure to update the params at the start of the nextflow to match your directory sturctures.

**Running only the annotation and LCA*
This pipeline can also be altered to just run the annotation and LCA. The params.EMMA needs to be changed to the path of the mitogenomes and in the workflow you need to comment out the first part, there are instructions in the nextflow script.

**Getorganelle**
The first process uses GetOrgnanelle to assemble the mitogenome using filtered and trimmed reads. The resulting mitogenome will be output as ${sample}.getorg1770.fasta.

**EMMA**
EMMA is an annotation program that is accessible at https://github.com/ian-small/Emma. This program uses the julia coding language and has been containerised for this step.

This process will give you the rotated fasta file starting at tRNA-Phe, the *.gff annotation file which is readable in genious prime for manual curation and visualisation, the *.tbl annotation file which is the file used when uploading to NCBI Genbank and a *.svg file for a nice visualisation of the circuralised mitogenome. In the cds directory is the coding sequences for all the proteins and the 12s and 16s RNA.

**BLAST**
The BLAST process blasts the mitogenome against the OceanOmics curated database which includes NCBI, BOLD and our own data. This is periodically updated so the run date of the BLAST is what will allow you to know if it has been updated since last using the database.

The 16s, 12s and COX1 sequences are used for this processes. The results are then filtered to remove any matches that are less than 200bp and less than 98% match.

**LCA**
This process uses the computeLCA.py script to find the lowest common ancestor using the fitered blast results.

# After the Nextflow is complete
When the nextflow is complete there is the 02_tidy_ilmn.sh which will remove the unnessarcery files for backing up. Some of these files may be used to improve the mitogenome so have been left there should any further interigation is needed. However, we do not want to back all these files up onto Acacia, so remove them before backing up.

After the tidy up script has been run you can then back up the mitogenomes onto Acacia. You can run the 03_backup_acacia script for this. This process moves all the files onto Acacia so after it has been run there will be no files left on your scratch.