import pandas as pd
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import os 

background_pop = Path("background_pop")
outdir = Path("final_background_pop")

print("Creating output directory...")
outdir.mkdir(parents=True, exist_ok=True)

print("Reading UniProt metadata...")
uniprot_info = pd.read_csv(
    "processed_README",
    sep='\t'
).drop(columns=["Proteome_ID", "OSCODE", "#(1)", "#(2)", "#(3)"]).rename(columns={
    "Species Name": "Species_Name"
})

uniprot_cel_orgs = uniprot_info[uniprot_info["SUPERREGNUM"] != "viruses"]

print("Reading GOA metadata...")
goa_info = pd.read_csv(
    "proteome2taxid",
    sep='\t',
    header=None,
    names=[
        "Species_Name",
        "Tax_ID",
        "GOA_file"
    ])

cellular_organisms = goa_info[goa_info['Tax_ID'].isin(uniprot_cel_orgs['Tax_ID'])]
cellular_organisms['index'] = cellular_organisms['GOA_file'].str.split('.').str[0]
index_df = cellular_organisms[['index', 'Tax_ID']]

print("Creating Tax ID dictionary...")
tax_id_dict = index_df.set_index('index')['Tax_ID'].to_dict()

column_names = [
    "DB", "DB_Object_ID", "DB_Object_Symbol", "Relation", "GO_ID",
    "DB:Reference_(|DB:Reference)", "Evidence_Code", "With_(or)_From",
    "Aspect", "DB_Object_Name", "DB_Object_Synonym_(|Synonym)",
    "DB_Object_Type", "Taxon_(|taxon)", "Date", "Assigned_By",
    "Annotation_Extension", "Gene_Product_Form_ID"
]

dtypes = {
    "DB": "category", "DB_Object_ID": "category", "DB_Object_Symbol": "category",
    "Relation": "object", "GO_ID": "category", "DB:Reference_(|DB:Reference)": "category",
    "Evidence_Code": "category", "With_(or)_From": "object", "Aspect": "category",
    "DB_Object_Name": "category", "DB_Object_Synonym_(|Synonym)": "category",
    "DB_Object_Type": "category", "Taxon_(|taxon)": "category", "Date": "int64",
    "Assigned_By": "category", "Annotation_Extension": "object",
    "Gene_Product_Form_ID": "object",
}

def create_background(file, outdir, taxon_id, column_names, dtypes):
    """
    Reads a GOA file, filters it, and saves the background population.
    Returns the taxon_id upon success or raises an exception on failure.
    """
    try:
        print(f"Processing {taxon_id} from {file.name}...")
        gaf = pd.read_csv(
            file,
            sep='\t',
            header=None,
            comment="!",
            names=column_names,
            dtype=dtypes,
            low_memory=False
        )

        filtered_gaf = gaf[
            (gaf["DB"] == "UniProtKB") &
            (~gaf["Relation"].str.contains('NOT', case=False, na=False)) &
            (gaf["DB_Object_Type"] == "protein")
        ]

        filtered_gaf = filtered_gaf[['DB_Object_ID', 'GO_ID', 'Evidence_Code']].drop_duplicates()

        mod_gaf = outdir.joinpath(f"{taxon_id}_background.txt")

        filtered_gaf.to_csv(mod_gaf, sep='\t', index=False, header=False)
        return taxon_id 
        
    except Exception as e:
        print(f"!!! ERROR processing {taxon_id} from {file.name}: {e}")
        raise Exception(f"Error in create_background for {taxon_id}") from e

def main():
    tasks_to_submit = []
    print(f"Scanning directory {background_pop} for .goa files...")
    for file in background_pop.iterdir():
        if file.suffix == ".goa":
            index = file.name.split(".")[0]
            if index in tax_id_dict:
                tax_id = tax_id_dict[index]
                tasks_to_submit.append({
                    "file": file,
                    "tax_id": tax_id
                })
            else:
                print(f"Warning: Index '{index}' from file {file.name} not found in tax_id_dict. Skipping.")

    if not tasks_to_submit:
        print("No .goa files found matching the tax_id dictionary. Exiting.")
        return

    print(f"Found {len(tasks_to_submit)} files to process.")

    num_workers = max(1, os.cpu_count() - 2) if os.cpu_count() else 4
    print(f"Starting processing with up to {num_workers} parallel workers...")

    futures = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for task in tasks_to_submit:
            
            future = executor.submit(
                create_background,
                task["file"],
                outdir,
                task["tax_id"],
                column_names, 
                dtypes        
            )
            futures.append(future)

        processed_count = 0
        total_tasks = len(futures)
        for future in as_completed(futures):
            try:
                completed_tax_id = future.result()
                processed_count += 1
                print(f"({processed_count}/{total_tasks}) Successfully processed Taxon ID: {completed_tax_id}")
            except Exception as exc:
                 processed_count += 1
                 print(f"({processed_count}/{total_tasks}) A task failed: {exc}")


    print(f"\n--- All processing finished. Processed {processed_count} files. ---")

if __name__ == "__main__":
    main()
