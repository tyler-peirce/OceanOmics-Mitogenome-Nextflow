import psycopg2
import pandas as pd
import os
import glob
import numpy as np  # Required for handling infinity values

# run using: singularity run $SING/psycopg2:0.1.sif python 04_push_mitogenome_results_to_sqldb.py
# PostgreSQL connection parameters
db_params = {
    'dbname': 'oceanomics',
    'user': 'postgres',
    'password': 'oceanomics',
    'host': '203.101.227.69',
    'port': 5432
}


# Access file paths


# File containing mitogenome data
mito_path = "/scratch/pawsey0964/tpeirce/_MITOGENOMES/OceanOmics-Mitogenome-Nextflow/mtdnastat.250521.tsv"  # Update this with the correct file path

# Import Mitogenome Data
print(f"Importing data from {mito_path}")

# Load and preprocess mito data
mito = pd.read_csv(mito_path, sep="\t")

# Replace NaN with None
mito = mito.replace({np.nan: None})

# Normalize column names (remove spaces & make lowercase)
mito.columns = mito.columns.str.strip().str.lower()

try:
    # Connect to PostgreSQL
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    row_count = 0  # Track number of processed rows

    for index, row in mito.iterrows():
        values = tuple(row.values)
        # Extract primary key values
        og_id, tech, seq_date, code = row['og_id'], row['tech'], row['seq_date'], row['code']

        # UPSERT: Insert if not exists, otherwise update the row
        upsert_query = """
        INSERT INTO mitogenome_data (
            og_id, tech, seq_date, code, stats, length, length_emma, seqlength_12s,
            seqlength_16s, seqlength_co1, cds_no, trna_no, rrna_no, status, genbank, rrna12s,
            rrna16s, atp6, atp8, cox1, cox2, cox3, cytb, nad1, nad2, nad3, nad4, nad4l, 
            nad5, mad6, trna_phe, trna_val, trna_leuuag, trna_leuuaa, trna_ile, trna_met, 
            trna_thr, trna_pro, trna_lys, trna_asp, trna_glu, trna_sergcu, trna_seruga, 
            trna_tyr, trna_cys, trna_trp, trna_ala, trna_asn, trna_gly, trna_arg, trna_his, 
            trna_gln
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, 
                %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT (og_id, tech, seq_date, code) 
        DO UPDATE SET
            stats = EXCLUDED.stats,
            length = EXCLUDED.length,
            length_emma = EXCLUDED.length_emma,
            seqlength_12s = EXCLUDED.seqlength_12s,
            seqlength_16s = EXCLUDED.seqlength_16s,
            seqlength_co1 = EXCLUDED.seqlength_co1,
            cds_no = EXCLUDED.cds_no,
            trna_no = EXCLUDED.trna_no,
            rrna_no = EXCLUDED.rrna_no,
            status = EXCLUDED.status,
            genbank = EXCLUDED.genbank,
            rrna12s = EXCLUDED.rrna12s,
            rrna16s = EXCLUDED.rrna16s,
            atp6 = EXCLUDED.atp6,
            atp8 = EXCLUDED.atp8,
            cox1 = EXCLUDED.cox1,
            cox2 = EXCLUDED.cox2,
            cox3 = EXCLUDED.cox3,
            cytb = EXCLUDED.cytb,
            nad1 = EXCLUDED.nad1,
            nad2 = EXCLUDED.nad2,
            nad3 = EXCLUDED.nad3,
            nad4 = EXCLUDED.nad4,
            nad4l = EXCLUDED.nad4l,
            nad5 = EXCLUDED.nad5,
            mad6 = EXCLUDED.mad6,
            trna_phe = EXCLUDED.trna_phe,
            trna_val = EXCLUDED.trna_val,
            trna_leuuag = EXCLUDED.trna_leuuag,
            trna_leuuaa = EXCLUDED.trna_leuuaa,
            trna_ile = EXCLUDED.trna_ile,
            trna_met = EXCLUDED.trna_met,
            trna_thr = EXCLUDED.trna_thr,
            trna_pro = EXCLUDED.trna_pro,
            trna_lys = EXCLUDED.trna_lys,
            trna_asp = EXCLUDED.trna_asp,
            trna_glu = EXCLUDED.trna_glu,
            trna_sergcu = EXCLUDED.trna_sergcu,
            trna_seruga = EXCLUDED.trna_seruga,
            trna_tyr = EXCLUDED.trna_tyr,
            trna_cys = EXCLUDED.trna_cys,
            trna_trp = EXCLUDED.trna_trp,
            trna_ala = EXCLUDED.trna_ala,
            trna_asn = EXCLUDED.trna_asn,
            trna_gly = EXCLUDED.trna_gly,
            trna_arg = EXCLUDED.trna_arg,
            trna_his = EXCLUDED.trna_his,
            trna_gln = EXCLUDED.trna_gln;
        """

        # Debugging Check
        print("\n🔹 DEBUG: Checking SQL Query and Values")
        print("Number of %s placeholders:", upsert_query.count("%s"))
        print("Number of values in row:", len(values))
        print("Values:", values)
        # Ensure placeholders match values before executing
        if len(values) != upsert_query.count("%s"):
            print("❌ ERROR: Mismatch in number of placeholders and values!")
            

        cursor.execute(upsert_query, values)
        row_count += 1  # Increment successful row counter

        print(f"✅ Successfully processed {row_count} rows. Committing changes...")
    
    

    conn.commit()
    print("Database update complete.")

except Exception as e:
    print(f"Error: {e}")

finally:
    if cursor:
        cursor.close()
    if conn:
        conn.close()
    print("🔹 Connection closed.")