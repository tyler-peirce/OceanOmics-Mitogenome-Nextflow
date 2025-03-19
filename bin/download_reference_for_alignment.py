import os
import ssl
from Bio import Entrez, SeqIO
from urllib.request import urlopen

# Disable SSL verification globally
ssl._create_default_https_context = ssl._create_unverified_context

# Set your NCBI email (required for API access)
Entrez.email = "tyler.peirce@uwa.edu.au"  # Replace with your email

# Input file with OG ID and species name (tab-separated)
species_file = "species_list.txt"

# Base directory for reference sequences
base_output_dir = "references"
os.makedirs(base_output_dir, exist_ok=True)

def search_ncbi_for_accession(species):
    """Search for the mitochondrial genome accession number of a species."""
    search_term = f"{species} mitochondrion complete genome"
    handle = Entrez.esearch(db="nuccore", term=search_term, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    
    if record["IdList"]:
        return record["IdList"][0]  # Return first accession number
    else:
        print(f"❌ No results found for {species}")
        return None

def fetch_protein_sequences(og_id, species, accession):
    """Fetch mitochondrial protein sequences from NCBI and save them."""
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    # Create the OG-specific references folder
    output_dir = os.path.join(base_output_dir, og_id)
    os.makedirs(output_dir, exist_ok=True)

    # Define output FASTA file
    output_fasta = os.path.join(output_dir, f"{species.replace(' ', '_')}.fasta")

    with open(output_fasta, "w") as out_fasta:
        for feature in record.features:
            if feature.type == "CDS":
                gene_name = feature.qualifiers.get("gene", ["unknown_gene"])[0]
                product_name = feature.qualifiers.get("product", ["unknown_product"])[0]

                # Get protein sequence if available
                if "translation" in feature.qualifiers:
                    protein_seq = feature.qualifiers["translation"][0]
                    out_fasta.write(f">{gene_name}|{product_name}|{accession}\n{protein_seq}\n")

    print(f"✅ {species} proteins saved to {output_fasta}")

# Read OG and species names from file
with open(species_file, "r") as f:
    species_list = [line.strip().split("\t") for line in f.readlines() if line.strip()]

# Process each OG-species pair
for og_id, species_name in species_list:
    print(f"🔍 Searching for {species_name} (OG: {og_id})...")
    reference_accession = search_ncbi_for_accession(species_name)

    if reference_accession:
        print(f"📌 Found accession: {reference_accession} for {species_name} (OG: {og_id})")
        fetch_protein_sequences(og_id, species_name, reference_accession)
    else:
        print(f"⚠️ Skipping {species_name} (OG: {og_id}) - No mitochondrial genome found.")

print("✅ Protein extraction complete!")
