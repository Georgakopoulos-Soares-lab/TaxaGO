# <img src="/logo.svg" alt="TaxaGO Logo" width="200"/>

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Rust](https://img.shields.io/badge/Rust-1.65+-orange.svg)](https://www.rust-lang.org/)

</div>

## Features

- **GOEA at Taxonomic level**: Perform enrichment analysis across different taxonomic levels (e.g. Kingdom, Family etc.)
- **Single-taxon Analysis**: Conduct GOEA for individual species 
- **Common Ancestor Analysis**: Identify shared ancestors between specified GO terms
- **High Performance**: Built with Rust for optimal speed and memory efficiency
- **Wide Species Support**: Compatible with 19,077 UniProtKB reference proteomes
- **Flexible Input**: Accepts either a directory with FASTA files containing UniProtKB protein accessions for each input species or a .csv where each column represent a different species.
- **Statistical Rigor**: Implements either Fisher's exact test of Hypergeometric test for enrichment analysis
- **Easy Integration**: Can be easily integrated into bioinformatics pipelines

## Prerequisites

Before installing TaxaGO, ensure you have the following:

- Rust
- Cargo package manager
- Git (for cloning the repository)
- 4GB RAM minimum (8GB recommended)
- Unix-like operating system (Linux, macOS) or Windows
- NPM (UPDATE)
- Mermaid (UPDATE)
- 

## Installation

Follow these steps to install TaxaGO:

```bash
# Clone the repository
git clone https://github.com/Georgakopoulos-Soares-lab/TaxaGO

# Navigate to the repository directory
cd TaxaGO

# Install using Cargo
cargo install --path .

# Verify installation
taxago --help
```

## Example usage

Simple GOEA for multiple species
```bash
taxago --study-pop example_study_pop.csv --out-dir results
```

Taxonomic GOEA at the Phylum level
```bash
taxago --group-results --taxonomic-level phylum --study-pop example_study_pop.csv --out-dir results 
```

### Common Ancestor Analysis
```bash
common_ancestor --terms GO:0016070,GO:0140187 --graph results
```

### Available Commands

#### Main Commands
- `taxago`: Perform GOEA analysis
- `common_ancestor`: Identify common ancestors between GO terms

#### Options for taxago 
- `--obo`: Path to the Gene Ontology file in OBO format. (default: $HOME/.cargo/taxago_assets/go.obo)
- `--study-pop`: Directory containing study popultion for each taxon in FASTA format or CSV file with the study population for each species. (required)
- `--background-pop`: Directory containing background populations. (default: $HOME/.cargo/taxago_assets/background_pop)
- `--out-dir`: Directory to write results for each taxon and the combined results for the taxonomic level (if specified). (required)
- `--propagate-counts`: Propagates GO term counts upwards the Ontology graph (from child to parent). (must be specified to propagate the counts)
- `--statistical-test`: Statistical test to use (fisher or hypergeometric). (default: fisher)
- `--min-protein-count`: Minimum protein count a GO Term must have to be processed. (default: 5)
- `--min-score`: Minimum score (log(odds ratio)) a GO Term must have to be written in the results. Keeps GO terms with score ≥ threshold. (default: 2.0)
- `--significance-threshold`: P-value / Q-value threshold to determine significant results. (default: 0.05)
- `--adjustment-method`: Method to adjust p-values for multiple test correction (bonferroni or bh or by or no). (default: bonferroni)
- `--group-results`: Combine results from all taxa into a single output. (must be specified to group the results)
- `--taxonomic-level`: Desired taxonomic level for result combination (superkingdom, kingdom, phylum, class, order, family, genus). (default: kingdom)
- `--lineage-percentage`: Percentage of species inside the desired taxonomic level in which the GO term must be found in (from 0.0 to 1.0). (default: 0.25)
- `--pm-iterations`: Number of maximum iterations the Paule-Mandel estimator can reach when calculating the τ² estimate. (default: 1000)
- `--pm-tolerance`: Minimum acceptable tolerance between two τ² estimates. (default: 1e-6)

#### Options for common_ancestor
- `--obo`: Path to the Gene Ontology file in OBO format. (default: $HOME/.cargo/taxago_assets/go.obo)
- `--terms`: Comma-separated list of GO terms (e.g., GO:0016070,GO:0140187). (required)
- `--graph`: Path to write the ancestor graph between the input terms. (default: current path)

## Documentation

### Input File Formats
1. **CSV file with the study population (proteins meant for the analysis for each unique species) (`study_pop.csv`):**

| Taxon ID 1 | Taxon ID 2 | Taxon ID 3 | Taxon ID 4 | Taxon ID 5 | ... | Taxon ID N |
|------------|------------|------------|------------|------------|-----|------------|
| PROTEIN 1  | PROTEIN 1  | PROTEIN 1  | PROTEIN 1  | PROTEIN 1  | ... | PROTEIN 1  |
| PROTEIN 2  | PROTEIN 2  | PROTEIN 2  | PROTEIN 2  | PROTEIN 2  | ... | PROTEIN 2  |
| PROTEIN 3  | PROTEIN 3  | PROTEIN 3  | PROTEIN 3  | PROTEIN 3  | ... | PROTEIN 3  |
| ...        | ...        | ...        | ...        | ...        | ... | ...        |
| PROTEIN K  | PROTEIN K  | PROTEIN K  | PROTEIN K  | PROTEIN K  | ... | PROTEIN K  |

2. **Directory containing FASTA files for each species used in the GOEA, where each FASTA is of this format:**
```text
>TaxonID
PROTEIN 1
PROTEIN 2
PROTEIN 3
PROTEIN 4
...
PROTEIN K
```

### Output Format

The tool generates results in CSV format with the following columns:

Single taxon results:

| GO Term ID | Name | Namespace | log(Odds Ratio) | Statistical significance | N Study with term | N Study without term | N Background with term | N Background without term |
|------------|------|-----------|---------|------------|------------|-------------|-----------|-----------|
| GO:0042391 | Regulation of membrane potential | Biological Process | 4.328 | 5.53201e-23 | 21 | 60 | 51 | 11045 |
| GO:0015075 | Monoatomic ion transmembrane transporter activity | Molecular Function | 2.005 | 2.39047e-4 | 20 | 61 | 469 | 10627 |
| GO:0098916 | Anterograde trans-synaptic signaling | Biological Process | 4.180 | 1.90007e-25 | 24 | 57 | 71 | 11025 |
| GO:0004888 | Transmembrane signaling receptor activity | Molecular Function | 2.798 | 2.10769e-20 | 37 | 44 | 541 | 10555 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |

Taxonomy results:

| GO Term ID | Name | Namespace | log(Odds Ratio) | Statistical significance | Heterogeneity | Species Percentage | N with GO term| N in taxonomy |
|------------|------|-----------|---------|------------|----------------|------------|----------|------------|
| GO:0015770 | Sucrose transport | Biological Process | 2.205 | 9.01783e-8 | 3.59790e-1 | 25.160 | 157 | 624 |
| GO:0043652 | Engulfment of apoptotic cell | Biological Process | 2.541 | 3.27587e-8 | 4.67309e-1 | 27.244 | 170 | 624 |
| GO:0004703 | G protein-coupled receptor kinase activity | Molecular Function | 2.056 | 3.74576e-8 | 1.20615e0 | 76.603 | 478 | 624 |
| GO:0004993 | G protein-coupled serotonin receptor activity | Molecular Function | 2.295 | 4.17601e-9 | 0.00000e0 | 98.237 | 613 | 624 |
| ... | ... | ... | ... | ... | ... | ... | ... | ... |


### Example Workflow

1. **Prepare Input Data**
   - Create a CSV file with protein UniProtKB identifiers
   - Column name must be the TaxonID
   - One protein per line
   - Supported formats: UniProt IDs

2. **Run Analysis**
   ```bash
   taxago \
   --group-results \
   --taxonomic-level phylum \
   --study-pop input_study_pop.csv \
   --out-dir results
   ```

3. **Interpret Results**
   - Review significant GO terms

## Contact

For support or questions:
- izg5139@psu.edu
- left.bochalis@gmail.com
- antonpapg@gmail.com 

