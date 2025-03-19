import os
from Bio import SeqIO

# Define the input directory
input_dir = "../mitogenomes"  # Change this to your actual directory

# Process all .fa files in the directory
for fa_file in sorted(os.listdir(input_dir)):
    if fa_file.endswith(".fa"):  # Only process FASTA files
        base_name = os.path.splitext(fa_file)[0]  # Get the basename (e.g., OG174)
        tbl_file = os.path.join(input_dir, f"{base_name}.tbl")
        fa_path = os.path.join(input_dir, fa_file)

        # Ensure the matching .tbl file exists
        if not os.path.exists(tbl_file):
            print(f"Skipping {fa_file} - No matching .tbl file found.")
            continue

        # Create an output folder named after the file's basename
        output_dir = os.path.join(input_dir, base_name)
        os.makedirs(output_dir, exist_ok=True)

        # Load genome sequence from .fa file
        genome_record = SeqIO.read(fa_path, "fasta")

        # Read the .tbl file and extract CDS information
        gene_name = "unknown_gene"  # Default gene name
        product_name = "unknown_product"  # Default product name

        with open(tbl_file, "r") as tbl:
            lines = tbl.readlines()  # Read all lines for proper multi-line parsing

            for i, line in enumerate(lines):
                parts = line.strip().split("\t")  # Split by tab
                
                # Check for gene name (stored before the CDS entry)
                if len(parts) == 2 and parts[0] == "gene":
                    gene_name = parts[1]  # Capture the gene name

                # If it's a CDS entry, extract the sequence
                elif len(parts) >= 3 and parts[2] == "CDS":
                    start, end = int(parts[0]), int(parts[1])  # Extract positions

                    # Extract product name from the next few lines (if available)
                    product_name = "unknown_product"
                    for j in range(i + 1, min(i + 4, len(lines))):  # Check the next few lines
                        if "product" in lines[j]:
                            product_name = lines[j].split("\t")[-1].strip()
                            break  # Stop once product is found

                    # Extract CDS sequence (handle reverse strand)
                    if start < end:
                        cds_seq = genome_record.seq[start - 1:end]  # 1-based indexing
                    else:
                        cds_seq = genome_record.seq[end - 1:start].reverse_complement()

                    # Define the output filename (e.g., MT-ND4.OG174.fa)
                    output_fasta = os.path.join(output_dir, f"{gene_name}.{base_name}.fa")

                    # Write to separate FASTA file
                    with open(output_fasta, "w") as out_fasta:
                        fasta_header = f">{gene_name}|{product_name}|{base_name}"
                        out_fasta.write(f"{fasta_header}\n{cds_seq}\n")

                    print(f"Extracted {gene_name} ({product_name}) -> {output_fasta}")

print("CDS extraction complete!")
