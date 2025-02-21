import argparse
from pathlib import Path
import pandas as pd

def load_goea_results(path, dictionary):
    for file in path.glob("*_GOEA_results.txt"):
            file_name = file.stem.replace("_GOEA_results", "")
            file_name = file_name.replace("_", " ")
            goea_results_df = pd.read_csv(file, sep='\t')  
            dictionary[file_name] = goea_results_df

def main():
    parser = argparse.ArgumentParser(description="Example script")
    parser.add_argument("-i", help="Directory containing the GOEA results.", required=True)
    
    args = parser.parse_args()

    goea_results_dir = Path(args.i)
    dirs = [entry.name for entry in goea_results_dir.iterdir() if entry.is_dir()]
    goea_results_dict = dict()

    if "combined_taxonomy_results" in dirs:
        combined_results_dir = goea_results_dir / "combined_taxonomy_results"
        load_goea_results(combined_results_dir, goea_results_dict)
        plots_dir = combined_results_dir / "plots"
    if "single_taxon_results" in dirs:
        single_taxon_results_dir = goea_results_dir / "single_taxon_results"
        load_goea_results(single_taxon_results_dir, goea_results_dict)
        plots_dir = single_taxon_results_dir / "plots"

if __name__ == "__main__":
    main()