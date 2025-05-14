import psycopg2
import pandas as pd
import os
import glob
import numpy as np  # Required for handling infinity values

# run using: singularity run $SING/psycopg2:0.1.sif python 04_push_results_to_sqldb.py

# PostgreSQL connection parameters
db_params = {
    'dbname': 'oceanomics',
    'user': 'postgres',
    'password': 'oceanomics',
    'host': '115.146.85.41',
    'port': 5432
}


# Access file paths
#DATE = config.get("DATE")
DATE = "250131"


# Define directory where files are stored
base_dir_pattern = "../mitogenomes/*/*/lca"
base_dirs = glob.glob(base_dir_pattern)
print(base_dirs)


## Inserting the filtered BLAST results to SQL database
# Define file patterns for the three regions
file_patterns = {
    "12s": "blast.12s.*.filtered.tsv",
    "16s": "blast.16s.*.filtered.tsv",
    "CO1": "blast.CO1.*.filtered.tsv"
}

# Define correct column names
column_headers = [
    "query_id", "match_sequence_id", "taxon_id", "scientific_name",
    "common_name", "kingdoms", "percent_identity", "alignment_length",
    "query_length", "subject_length", "mismatch", "gap_open", "gaps",
    "query_start", "query_end", "subject_start", "subject_end",
    "subject_title", "evalue", "bit_score", "query_coverage", "subject_coverage",
    "region"  # New column
]

# Create an empty list to store DataFrames
merged_data = []

for base_dir in base_dirs:
    # Loop through each file pattern
    for region, pattern in file_patterns.items():
        # Get all matching files
        files = glob.glob(os.path.join(base_dir, pattern))
        print(files)
        print(region)
        print(pattern)
        
        for file in files:
            # Skip empty files
            if os.stat(file).st_size == 0:
                    print(f"Skipping empty file: {file}")
                    continue
            try:
            # Read file with tab separator
                df = pd.read_csv(file, sep='\t', header=None)
            
                # Ensure correct number of columns
                if df.shape[1] == len(column_headers) - 1:  # Excluding 'region'
                    df.columns = column_headers[:-1]  # Assign headers (excluding 'region')
                    df["region"] = region  # Add 'region' column
                    merged_data.append(df)
                else:
                    print(f"Skipping {file} due to column mismatch")
            
            except Exception as e:
                    print(f"Error reading {file}: {e}")

# Concatenate all DataFrames if any valid data exists
if merged_data:
    merged_df = pd.concat(merged_data, ignore_index=True)

    if 'query_id' in merged_df.columns:
        # Split 'query_id' into three new columns
        merged_df[['og_id', 'tech', 'seq_date', 'code', 'annotation']] = merged_df['query_id'].str.split('.', expand=True)
        
        print("File successfully processed! New columns added.")
    else:
        print("Error: 'query_id' column not found in the input file.")

    # Save the merged results
    output_file = "merged_blast_results.tsv"
    merged_df.to_csv(output_file, sep="\t", index=False)

    print(f"Merged file saved as {output_file}")
else:
    print("No valid files found for merging.")






try:
    # Connect to PostgreSQL
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    row_count = 0  # Track number of processed rows

    for index, row in merged_df.iterrows():
        row_dict = row.to_dict()

        # Extract primary key values
        og_id, tech, seq_date, code, annotation, match_sequence_id, region = row_dict["og_id"], row_dict["tech"], row_dict["seq_date"], row_dict["code"], row_dict["annotation"], row_dict["match_sequence_id"], row_dict["region"]

        # UPSERT: Insert if not exists, otherwise update
        upsert_query = """
        INSERT INTO blast_filtered_lca (
            og_id, tech, seq_date, code, annotation, match_sequence_id, taxon_id,
            scientific_name, common_name, kingdoms, percent_identity, alignment_length,
            query_length, subject_length, mismatch, gap_open, gaps, query_start,
            query_end, subject_start, subject_end, subject_title, evalue, bit_score,
            query_coverage, subject_coverage, region
        )
        VALUES (
            %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(annotation)s, %(match_sequence_id)s, %(taxon_id)s,
            %(scientific_name)s, %(common_name)s, %(kingdoms)s, %(percent_identity)s, %(alignment_length)s,
            %(query_length)s, %(subject_length)s, %(mismatch)s, %(gap_open)s, %(gaps)s, %(query_start)s,
            %(query_end)s, %(subject_start)s, %(subject_end)s, %(subject_title)s, %(evalue)s, %(bit_score)s,
            %(query_coverage)s, %(subject_coverage)s, %(region)s
        )   
        ON CONFLICT (og_id, tech, seq_date, code, annotation, match_sequence_id, region)
        DO UPDATE SET
            taxon_id = EXCLUDED.taxon_id,
            scientific_name = EXCLUDED.scientific_name,
            common_name = EXCLUDED.common_name,
            kingdoms = EXCLUDED.kingdoms,
            percent_identity = EXCLUDED.percent_identity,
            alignment_length = EXCLUDED.alignment_length,
            query_length = EXCLUDED.query_length,
            subject_length = EXCLUDED.subject_length,
            mismatch = EXCLUDED.mismatch,
            gap_open = EXCLUDED.gap_open,
            gaps = EXCLUDED.gaps,
            query_start = EXCLUDED.query_start,
            query_end = EXCLUDED.query_end,
            subject_start = EXCLUDED.subject_start,
            subject_end = EXCLUDED.subject_end,
            subject_title = EXCLUDED.subject_title,
            evalue = EXCLUDED.evalue,
            bit_score = EXCLUDED.bit_score,
            query_coverage = EXCLUDED.query_coverage,
            subject_coverage = EXCLUDED.subject_coverage;
        """
        
        # Params as a dictionary as the columns arent in the right order
        params = {
            "og_id": row_dict["og_id"],
            "tech": row_dict["tech"],
            "seq_date": row_dict["seq_date"],
            "code": row_dict["code"],
            "annotation": row_dict["annotation"],
            "match_sequence_id": row_dict["match_sequence_id"],
            "taxon_id": row_dict["taxon_id"] if row_dict["taxon_id"] else None,
            "scientific_name": row_dict["scientific_name"],
            "common_name": row_dict["common_name"],
            "kingdoms": row_dict["kingdoms"],
            "percent_identity": row_dict["percent_identity"] if row_dict["percent_identity"] else None,
            "alignment_length": row_dict["alignment_length"] if row_dict["alignment_length"] else None,
            "query_length": row_dict["query_length"] if row_dict["query_length"] else None,
            "subject_length": row_dict["subject_length"] if row_dict["subject_length"] else None,
            "mismatch": row_dict["mismatch"] if row_dict["mismatch"] else 0,
            "gap_open": row_dict["gap_open"] if row_dict["gap_open"] else 0,
            "gaps": row_dict["gaps"] if row_dict["gaps"] else 0,
            "query_start": row_dict["query_start"] if row_dict["query_start"] else None,
            "query_end": row_dict["query_end"] if row_dict["query_end"] else None,
            "subject_start": row_dict["subject_start"] if row_dict["subject_start"] else None,
            "subject_end": row_dict["subject_end"] if row_dict["subject_end"] else None,
            "subject_title": row_dict["subject_title"],
            "evalue": row_dict["evalue"] if row_dict["evalue"] else 0,
            "bit_score": row_dict["bit_score"] if row_dict["bit_score"] else None,
            "query_coverage": row_dict["query_coverage"] if row_dict["query_coverage"] else None,
            "subject_coverage": row_dict["subject_coverage"] if row_dict["subject_coverage"] else None,
            "region": row_dict["region"],
        }

        # Debugging Check
        #print(f"Number of columns in query: {upsert_query.count('%s')}")
        print(f"Number of rows being passed: {len(row)}")
        print(f"Column names in DataFrame: {merged_df.columns.tolist()}")
        print("row:", row_dict)
        print("params:", params)

        cursor.execute(upsert_query, params)
        row_count += 1  

        conn.commit()
        print(f"✅ Successfully processed {row_count} rows!")

except Exception as e:
    conn.rollback()
    print(f"❌ Error: {e}")

finally:
    cursor.close()
    conn.close()

print("Connection Closed")



## Inserting the LCA results to SQL database
# Define file patterns for the three regions
file_patterns = {"lca.*.tsv"}

# Define correct column names
column_headers = [
    "query_id", "taxonomy", "lca", "percent_match",
    "length", "lca_run_date", "region"
]

# Create an empty list to store DataFrames
merged_data = []

for base_dir in base_dirs:
    # Loop through each file pattern
    for pattern in file_patterns:
        # Get all matching files
        files = glob.glob(os.path.join(base_dir, pattern))
        print(f"files: {files}")
        print(f"pattern: {pattern}")
        
        for file in files:
            # Skip empty files
            if os.stat(file).st_size == 0:
                    print(f"Skipping empty file: {file}")
                    continue
            try:
            # Read file with tab separator
                df = pd.read_csv(file, sep='\t', header=None)
            
                # Ensure correct number of columns
                if df.shape[1] == len(column_headers):
                    df.columns = column_headers  # Assign headers
                    merged_data.append(df)
                else:
                    print(f"Skipping {file} due to column mismatch")
            
            except Exception as e:
                    print(f"Error reading {file}: {e}")

# Concatenate all DataFrames if any valid data exists
if merged_data:
    merged_df = pd.concat(merged_data, ignore_index=True)

    if 'query_id' in merged_df.columns:
        # Split 'query_id' into three new columns
        merged_df[['og_id', 'tech', 'seq_date', 'code', 'annotation']] = merged_df['query_id'].str.split('.', expand=True)
        
        print("File successfully processed! New columns added.")
    else:
        print("Error: 'query_id' column not found in the input file.")

    # Save the merged results
    output_file = "merged_lca_results.tsv"
    merged_df.to_csv(output_file, sep="\t", index=False)

    print(f"Merged file saved as {output_file}")
else:
    print("No valid files found for merging.")

try:
    # Connect to PostgreSQL
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    row_count = 0  # Track number of processed rows

    for index, row in merged_df.iterrows():
        row_dict = row.to_dict()

        # Extract primary key values
        og_id, tech, seq_date, code, annotation, region = row_dict["og_id"], row_dict["tech"], row_dict["seq_date"], row_dict["code"], row_dict["annotation"], row_dict["region"]

        # UPSERT: Insert if not exists, otherwise update
        upsert_query = """
        INSERT INTO lca (
            og_id, tech, seq_date, code, annotation, taxonomy, lca, percent_match,
            length, lca_run_date, region
        )
        VALUES (
            %(og_id)s, %(tech)s, %(seq_date)s, %(code)s, %(annotation)s, %(taxonomy)s, %(lca)s,
            %(percent_match)s, %(length)s, %(lca_run_date)s, %(region)s
        )   
        ON CONFLICT (og_id, tech, seq_date, code, annotation, region)
        DO UPDATE SET
            taxonomy = EXCLUDED.taxonomy,
            lca = EXCLUDED.lca,
            percent_match = EXCLUDED.percent_match,
            length = EXCLUDED.length,
            lca_run_date = EXCLUDED.lca_run_date;
        """
        
        # Params as a dictionary as the columns arent in the right order
        params = {
            "og_id": row_dict["og_id"],
            "tech": row_dict["tech"],
            "seq_date": row_dict["seq_date"],
            "code": row_dict["code"],
            "annotation": row_dict["annotation"],
            "taxonomy": row_dict["taxonomy"],
            "lca": row_dict["lca"],
            "percent_match": row_dict["percent_match"],
            "length": row_dict["length"],
            "lca_run_date": row_dict["lca_run_date"],
            "region": row_dict["region"],
        }

        # Debugging Check
        #print(f"Number of columns in query: {upsert_query.count('%s')}")
        print(f"Number of rows being passed: {len(row)}")
        print(f"Column names in DataFrame: {merged_df.columns.tolist()}")
        print("row:", row_dict)
        print("params:", params)

        cursor.execute(upsert_query, params)
        row_count += 1  

        conn.commit()
        print(f"✅ Successfully processed {row_count} rows!")

except Exception as e:
    conn.rollback()
    print(f"❌ Error: {e}")

finally:
    cursor.close()
    conn.close()

print("Connection Closed")