import os
import glob
import subprocess
from Bio import SeqIO

# Directories
mitogenomes_dir = "../mitogenomes"
references_dir = "references"
output_dir = "../aligned_proteins"
aligned_proteins_dir = "../aligned_proteins"
mafft_container = "/software/projects/pawsey0812/singularity/mafft:7.515.sif"  # Singularity MAFFT container

singularity_path = "/software/setonix/2024.05/software/linux-sles15-zen3/gcc-12.2.0/singularityce-4.1.0-2gadr2xoc2nb4prnnyq2vvztjh6x4wzl/bin/singularity"


# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Loop through sample directories in mitogenomes/
for sample_folder in sorted(os.listdir(mitogenomes_dir)):
    sample_path = os.path.join(mitogenomes_dir, sample_folder)

    if not os.path.isdir(sample_path):
        continue  # Skip if not a directory

    # Extract the sample ID (e.g., OG228 from OG228.ilmn.230324.v177getorg.emma102_proteins)
    sample_id = sample_folder.split('.')[0]
    # Locate the reference file dynamically
    reference_files = glob.glob(os.path.join(references_dir, sample_id, "*.fasta"))

    if not reference_files:
        print(f"⚠️ No reference found for {sample_id}, skipping...")
        continue

    # Use the first found reference file
    reference_path = reference_files[0]

    # Load reference protein sequences into a dictionary
    reference_sequences = {}
    with open(reference_path, "r") as ref_file:
        for record in SeqIO.parse(ref_file, "fasta"):
            gene_name = record.description.split('|')[0]  # Extract gene name (e.g., "ND1")
            reference_sequences[gene_name] = record

    # Loop through translated protein files in the sample directory
    for protein_file in sorted(os.listdir(sample_path)):
        if not protein_file.endswith(".fa"):
            continue  # Skip non-FASTA files

        # Extract gene name from filename (e.g., "MT-ATP6" from "MT-ATP6.OG228...")
        gene_name = protein_file.split('.')[0].replace("MT-", "")  # Normalize gene name

        # Ensure the reference has this protein
        if gene_name not in reference_sequences:
            print(f"⚠️ No reference match found for {gene_name} in {sample_id}")
            continue

        # Load sample protein sequence
        sample_protein_path = os.path.join(sample_path, protein_file)
        with open(sample_protein_path, "r") as sample_file:
            sample_record = next(SeqIO.parse(sample_file, "fasta"))

        # Get corresponding reference sequence
        reference_record = reference_sequences[gene_name]

        # Create concatenated output
        output_filename = f"{gene_name}.{sample_id}.aligned.fa"
        output_path = os.path.join(output_dir, output_filename)

        with open(output_path, "w") as out_file:
            SeqIO.write(sample_record, out_file, "fasta")
            SeqIO.write(reference_record, out_file, "fasta")

        print(f"✅ Created {output_filename}")

print("🚀 All files processed! Ready for MAFFT.")

