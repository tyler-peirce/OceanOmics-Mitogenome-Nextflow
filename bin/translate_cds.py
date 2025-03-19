import os
from Bio import SeqIO
from Bio.Seq import Seq

# Define base directory
basedir = "../mitogenomes"

# Walk through all matching CDS files
for root, _, files in os.walk(basedir):
    #if "OG" in root:  # Ensure we are in the cds directory. This lines is for normal file structure
    if os.path.basename(root).startswith("OG"):  # Only process folders starting with "OG". For trend stuff
        for file in files:
            if file.endswith(".fasta") or file.endswith(".fa"):  # Ensure it's a FASTA file
                cds_path = os.path.join(root, file)
                
                # Define output directory
                #proteins_dir = root.replace("cds", "proteins") # this line is for normal file structure
                proteins_dir = root + "_proteins" # for trend stuff
                os.makedirs(proteins_dir, exist_ok=True)
                
                # Define output file path
                protein_file = os.path.join(proteins_dir, file)

                # Translate and save proteins
                with open(protein_file, "w") as out_fasta:
                    for record in SeqIO.parse(cds_path, "fasta"):
                        protein_seq = Seq(record.seq).translate(table=2)  # Table 2 = Vertebrate mitochondrial code
                        out_fasta.write(f">{record.description}\n{protein_seq}\n")

                print(f"Translated {cds_path} -> {protein_file}")
