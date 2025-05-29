<p align="center">
  <img src="logo.svg" alt="TaxaGO Logo" width="500"/>
</p>

---

# TaxaGO: a novel, high-performance, multi-taxonomic Gene Ontology Enrichment Analysis tool

<p align="left">
  <a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" alt="License: GPL v3"></a>
  <a href="https://www.rust-lang.org"><img src="https://img.shields.io/badge/Rust-1.87+-orange.svg" alt="Rust Version"></a>
  <img src="https://img.shields.io/badge/Version-v1.0.0-green.svg" alt="TaxaGO Version"></a>
  </p>


## Table of Contents

1.  [Graphical Abstract](#1-graphical-abstract)
2.  [Key Features](#2-key-features)
3.  [Installation](#3-installation)
    * [Prerequisites](#prerequisites)
    * [Required Assets](#required-assets)
    * [From Source](#from-source)
4.  [Usage](#4-usage)
    * [Phylogenetic GO Enrichment Analysis](#taxago)
    * [Semantic Similarity](#semantic-similarity)
    * [Common Ancestor Analysis](#common-ancestors)
    * [Interactive Interface](#taxago-interactive)
5.  [Input File Formats](#4-input-file-formats)
    * [OBO File](#obo-file)
    * [Study Population](#study-population)
    * [Background Population](#background-population)
    * [Lineage File](#lineage-file)
    * [Variance-Covariance (VCV) Matrix](#variance-covariance-vcv-matrix)
6.  [Output File Formats](#6-output-file-formats)
7.  [Interpreting Results](#7-interpreting-results)
8.  [Contributing](#8-contributing)
9.  [License](#9-license)
10. [Citation](#10-citation)
11. [Contact](#11-contact)

## 1. Graphical Abstract

<p align="center">
  <img src="abstract.png" alt="TaxaGO Logo" width="1000"/>
</p>

## 2. Key Features

* **Single-species GOEA :** Performs standard Gene Ontology Enrichment Analysis for a single species or for multiple species at once.
* **Phylogenetically-aware meta-analysis:** Unifying enrichment scores across different taxonomic levels, considering the evolutionary relationships between species.
* **Advanced GO hierarchy handling:** Implementing various count propagation algorithms (Classic, Elim, Weight) to refine enrichment signals.
* **Comprehensive GO toolkit:** Including semantic similarity calculations and common ancestor analysis to further explore GO term relationships.
* **User-friendly interfaces:** Offers a locally hosted interactive user interface for easier data input, parameter tuning and results exploration.
* **High Performance:** Designed for speed and can leverage multiple CPU cores to efficiently handle GOEA across multiple species and/or taxonomic levels simultaneously.

## 3. Installation

### Prerequisites

* **Rust Toolchain:** Version 1.87.0 or later is recommended. Install from the original [Rust website](https://www.rust-lang.org/tools/install).
* **Mermaid CLI (`mmdc`):** Required **only** for the `common-ancestors` tool to generate PDF outputs from Mermaid diagrams. Installation from the official [Mermaid CLI repository](https://github.com/mermaid-js/mermaid-cli).
* **Jemalloc:** TaxaGO uses `jemallocator` for potentially better memory allocation performance.

### Required Assets

TaxaGO uses several data files for its operations. We provide a pre-compiled `taxagp_assets.tar.gz` containing the ontology information, pre-processed background populations for 12,131 species with their corresponding taxonomic and phylogenetic information. The pre-compiled file can be downloaded from [Zenodo](test).

If you want to use your own data, here is a brief description of each one:
* **`go.obo`**: The Gene Ontology OBO file, can be downloaded from the [Gene Ontology Consortium](http://geneontology.org/docs/download-ontology/).
* **`background_pop/`**: A directory containing pre-processed background population files. Each file should be named `{taxon_id}_background.txt` (e.g., `9606_background.txt`).
* **`lineage.txt`**: A tab-separated file mapping NCBI Taxon IDs to their full taxonomic lineage.
* **`vcv.dmat`**: A Variance-Covariance matrix displaying the evolutionary relationship between species. 
See [Input File Formats](#input-file-formats) for additional details.

### From Source

1.  **Install prerequisites.**

2.  **Clone the repository:**
    ```bash
    git clone https://github.com/Georgakopoulos-Soares-lab/TaxaGO
    cd TaxaGO
    ```
3.  **Download `taxago_assets.tar.gz` from Zenodo.**

4.  **Move `taxago_assets.tar.gz` inside the cloned repository.**

5.  **Install TaxaGO:**
    ```bash
    cargo install --path .
    ```
    After installation a `taxago_assets` directory will be created in the `$CARGO_HOME` (specified during Rust installation).

6.  **Once installed, the following executables will be available in your system's PATH:**
    * `taxago`: Main executable for the GOEA analyses.
    * `semantic-similarity`: Tool for calculating GO term semantic similarity.
    * `common-ancestors`: Tool for finding and visualizing common GO ancestors.
    * `taxago-interactive`: Interactive user-interface executable.

## 4. Usage

TaxaGO provides a suite of tools for Gene Ontology Enrichment Analysis. The main executables are `taxago`, `semantic-similarity`, `common-ancestors`, and `taxago-interactive`.

### Gene Ontology Enrichment Analysis: taxago

The `taxago` executable is the primary tool for performing Gene Ontology Enrichment Analysis across single species or taxonomic levels.

### Synopsis:

```bash
taxago [OPTIONS] --study <FILE_OR_DIR> --dir <DIRECTORY>
```
### Options:

**Input Files**
- `-o, --obo <FILE>`: Gene Ontology file in OBO format  
  **Default:** `$CARGO_HOME/taxago_assets/go.obo`

- `-s, --study <FILE_OR_DIRECTORY>`: **Required.** Study population data. Accepts FASTA format (single file for one species, or directory of files for multi-species analysis) or a CSV file containing study populations for one or multiple species

- `-b, --background <DIRECTORY>`: Background population data. Either a single file for custom background or a directory containing background population files for multiple species. Background files must be pre-processed 
  **Default:** `$CARGO_HOME/taxago_assets/background_pop`

**Analysis Parameters**
- `-e, --evidence <CATEGORY>`: Evidence code categories to include from background associations  
  **Options:** `all`, `experimental`, `phylogenetic`, `computational`, `author`, `curator`, `automatic`  
  **Default:** `all`

- `-p, --propagate-counts <METHOD>`: Method for propagating GO term counts up the ontology hierarchy  
  **Options:** `none`, `classic`, `elim`, `weight`  
  **Default:** `none`

- `-t, --test <TEST>`: Statistical test for enrichment analysis  
  **Options:** `fishers`, `hypergeometric`  
  **Default:** `fishers`

**Filtering Thresholds**
- `-m, --min-prot <COUNT>`: Minimum number of proteins required for a GO term to be analyzed. GO terms with associations less than this number will be excluded
  **Default:** `5`

- `-r, --min-score <SCORE>`: Minimum log(Odds Ratio) threshold for GO terms to be reported or further analyzed. GO terms with observed log(Odds Ratio) less than this will be excluded.
  **Default:** `0.2`

- `-a, --alpha <THRESHOLD>`: Statistical significance threshold. Refers to either the corrected or uncorrected p-value
  **Default:** `0.05`

- `-c, --correction-method <METHOD>`: Multiple testing correction method  
  **Options:** `none`, `bonferroni`, `benjamini-hochberg`, `benjamini-yekutieli`  
  **Default:** `bonferroni`

**Meta-Analysis Options**
- `-g, --group-results <LEVEL>`: Group results by taxonomic level to be subjected to  phylogenetic meta-analysis. 
   **Requires** `--vcv-matrix`

- `-l, --lineage-percentage <PERCENTAGE>`: Minimum percentage (range 0.0 to 1.0) of species within a taxonomic group where a GO term must be found enriched 
  **Default:** `0.25` (25%)

- `--vcv-matrix <FILE>`: Variance-covariance matrix file for phylogenetic meta-analysis

- `--permutations <COUNT>`: Number of permutations for phylogenetic meta-analysis  
  **Default:** `1000`

**Output Options**
- `-d, --dir <DIRECTORY>`: **Required.** Output directory for results (individual taxon results and combined analysis). Previous results will be overwritten

- `--save-plots <FORMAT>`: Format for saving enrichment plots. `interactive`: HTML format, `static`: PDF format
  **Options:** `none`, `interactive`, `static`, `both`  
  **Default:** `interactive`

**System Options**
- `--cores <NUMBER>`: Number of CPU cores to use for parallel processing  
  **Default:** All available cores

- `-h, --help`: Display help information
- `-V, --version`: Display version information

### Example:

```bash
taxago -s ./my_study_data/ -d ./taxago_results/ -g kingdom --vcv-matrix ./assets/vcv_matrix.dmat -p classic -c benjamini-hochberg -a 0.01 --save-plots both
```

This command runs TaxaGO using study data from `./my_study_data/`, outputs results to `./taxago_results/`, combines results at the kingdom level using the VCV matrix from `./assets/vcv_matrix.dmat`, uses the classic count propagation, Benjamini-Hochberg for p-value correction with an alpha of 0.01, and saves both interactive HTML and static plots.

### Semantic Similarity: semantic-similarity

Calculates semantic similarity between GO terms.

### Synopsis:

```bash
semantic-similarity [OPTIONS] --terms <GO_TERMS_OR_FILE>
```

### Options:

**Input Files**
- `-o, --obo <OBO_FILE>`: Gene Ontology file in OBO format  
  **Default:** `$CARGO_HOME/taxago_assets/go.obo`

- `-t, --terms <GO_TERMS_OR_FILE>`: **Required.** GO terms to analyze. Either comma-separated terms (e.g., `GO:0016070,GO:0140187`) or path to a file containing one term per line

- `-b, --background <BACKGROUND_DIR>`: Directory containing background population files  
  **Default:** `$CARGO_HOME/taxago_assets/background_pop`

**Analysis Parameters**
- `-i, --ids <TAXON_IDS>`: Comma-separated list of taxonomic IDs to analyze (e.g., `9606,10090`)  
  **Default:** `9606` (Homo Sapiens)

- `-m, --method <METHOD>`: Semantic similarity calculation method  
  **Options:** `resnik`, `lin`, `jiang-conrath`, `wang`  
  **Default:** `resnik`

- `-p, --propagate-counts`: Propagate GO term counts up the ontology hierarchy.
  **Default:** Disabled

**Output Options**
- `-d, --dir <RESULTS_DIR>`: Directory for output files  
  **Default:** `./` (current directory)

- `-h, --help`: Display help information

### Example:

```bash
semantic-similarity -t "GO:0008150,GO:0005575" -i 9606 -m lin -d ./similarity_results/ --propagate-counts
```

This calculates Lin semantic similarity for GO:0008150 and GO:0005575 in taxon 9606, using propagated counts, and saves results to `./similarity_results/`.

### Common Ancestor Analysis: common-ancestors

Finds and visualizes common ancestors of specified GO terms.

### Synopsis:

```bash
common-ancestors [OPTIONS] --terms <GO_TERMS>
```

### Options:

**Input Options**
- `-o, --obo <OBO_FILE>`: Gene Ontology file in OBO format  
  **Default:** `$CARGO_HOME/taxago_assets/go.obo`

- `-t, --terms <GO_TERMS>`: **Required.** Comma-separated list of GO terms (e.g., `GO:0016070,GO:0140187`)

**Output Options**
- `-d, --dir <RESULTS_DIR>`: Output directory for results (generates Mermaid .mmd and .pdf files)  
  **Default:** `./` (current directory)

- `-h, --help`: Display help information

### Example:

```bash
common-ancestors -t "GO:0044237,GO:0006412" -d ./ancestor_analysis/
```

This command analyzes common ancestors for GO:0044237 and GO:0006412 and outputs the Mermaid graph and PDF to the `./ancestor_analysis/` directory.

### Interactive Interface: taxago-interactive

Launches a web-based interactive user interface for TaxaGO.

### Synopsis:

```bash
taxago-interactive
```

### Usage:

1. Run the executable:
   ```bash
   taxago-interactive
   ```

2. A browser window will be opened with the interactive interface. If not, please open your web browser and navigate to `http://127.0.0.1:8080`.

The interactive interface allows you to:

- Upload study population data (CSV or FASTA).
- Upload custom background population files.
- Upload a custom Gene Ontology (OBO) file.
- Configure analysis parameters (statistical test, count propagation, multiple testing correction, thresholds).
- Configure and run Phylogenetic Meta-Analysis (PMA) including VCV matrix upload.
- Filter by evidence codes.
- Initiate the analysis and view results, including interactive plots.
- Download results in various formats.
- Perform Common Ancestor and Semantic Similarity analyses through the UI.


## 5. Input File Formats

## 6. Output File Formats

## 7. Interpreting Results

## 8. Contributing

We warmly welcome contributions to TaxaGO! Whether it's reporting a bug, suggesting a new feature, improving documentation, or submitting code changes, your help is greatly appreciated and valued.

To ensure a smooth and effective collaboration process, please consider the following guidelines:

**Ways to Contribute:**

* **Reporting Bugs:** If you encounter a bug, please open an issue on our [GitHub Issues page](https://github.com/Georgakopoulos-Soares-lab/TaxaGO/issues). Describe the bug in detail, including steps to reproduce it, the expected behavior, and the actual behavior. Include your operating system, Rust version.
* **Suggesting Enhancements or New Features:** We are always open to new ideas! Please open an issue on GitHub to suggest an enhancement or new feature. Provide a clear and detailed explanation of the feature and why it would be beneficial to TaxaGO.
* **Improving Documentation:** Good documentation is key. If you find areas that are unclear, incorrect, or could be improved, please let us know by opening an issue or submitting a pull request with your suggested changes.
* **Submitting Code Changes (Pull Requests):**
    1.  **Fork the repository** on GitHub.
    2.  **Create a new branch** for your feature or bug fix: `git checkout -b feature/your-feature-name` or `git checkout -b fix/your-bug-fix-name`.
    3.  **Make your changes.** Ensure your code adheres to the existing style and that you add relevant tests.
    4.  **Test your changes thoroughly.**
    5.  **Commit your changes** with a clear and descriptive commit message: `git commit -m "feat: Add new feature X"`.
    6.  **Push your branch** to your forked repository: `git push origin feature/your-feature-name`.
    7.  **Open a Pull Request (PR)** against the `main` branch of the `Georgakopoulos-Soares-lab/TaxaGO`repository.
    8.  Clearly describe the changes in your PR, why they were made, and reference any related issues.

We look forward to your contributions and to making TaxaGO a better tool together!

## 9. License

This project is licensed under the **GNU GPL v3**.

See the [LICENSE.txt](LICENSE.txt) file for further details.

## 10. Citation

The citation will be placed here after publication.

## 11. Contact
For any questions or support, please contact:
* izg5139@psu.edu
* left.bochalis@gmail.com
* antonpapg@gmail.com