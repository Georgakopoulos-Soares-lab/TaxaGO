use clap::Parser;
use std::error::Error;
use std::fs::create_dir_all;
use std::env::var;
use std::collections::HashSet;
use TaxaGO::parsers::{
    obo_parser::*,
    study_parser::read_study_pop,
    background_parser::{BackgroundPop,read_background_pop},
};
use lazy_static::lazy_static;
use std::path::PathBuf;
use dirs::home_dir;

use TaxaGO::analysis::{
    enrichment_analysis::*, 
    multiple_testing_correction::*, 
    write_results::*,
    handle_lineage::*,
    result_combination::*,

};

lazy_static! {
    static ref DEFAULT_OBO_PATH: String = {
        let cargo_home = var("CARGO_HOME")
            .unwrap_or_else(|_| {
                home_dir()
                    .expect("Could not determine home directory")
                    .join(".cargo")
                    .to_string_lossy()
                    .into_owned()
            });
        PathBuf::from(cargo_home)
            .join("taxago_assets")
            .join("go.obo")
            .to_string_lossy()
            .into_owned()
    };
}

lazy_static! {
    static ref DEFAULT_BACKGROUND: String = {
        let cargo_home = var("CARGO_HOME")
            .unwrap_or_else(|_| {
                home_dir()
                    .expect("Could not determine home directory")
                    .join(".cargo")
                    .to_string_lossy()
                    .into_owned()
            });
        PathBuf::from(cargo_home)
            .join("taxago_assets")
            .join("background_pop")
            .to_string_lossy()
            .into_owned()
    };
}

lazy_static! {
    static ref DEFAULT_LINEAGE: String = {
        let cargo_home = var("CARGO_HOME")
            .unwrap_or_else(|_| {
                home_dir()
                    .expect("Could not determine home directory")
                    .join(".cargo")
                    .to_string_lossy()
                    .into_owned()
            });
        PathBuf::from(cargo_home)
            .join("taxago_assets")
            .join("full_lineage.txt")
            .to_string_lossy()
            .into_owned()
    };
}

#[derive(Parser, Debug)]
#[command(name = "taxago")]
struct CliArgs {
    #[arg(
        long = "obo",
        value_name = "OBO_FILE",
        help = "Path to the Gene Ontology file in OBO format.",
        default_value_t = DEFAULT_OBO_PATH.to_string(),
    
    )]
    obo_file: String,
    
    #[arg(
        long = "study-pop",
        value_name = "STUDY_DIR",
        help = "Directory containing study popultion for each taxon in FASTA format or CSV file with the study population for each species.",
        required = true
    )]
    study_dir: String,
    
    #[arg(
        long = "background-pop",
        value_name = "BACKGROUND_DIR",
        help = "Directory containing background populations.",
        default_value_t = DEFAULT_BACKGROUND.to_string(),
    )]
    background_dir: String,
    
    #[arg(
        long = "out-dir",
        value_name = "RESULTS_DIR",
        help = "Directory to write results for each taxon and the combined results for the taxonomic level [if specified].",
        required = true
    )]
    output_dir: String,
    
    #[arg(
        long = "propagate-counts",
        help = "Propagates GO term counts upwards the Ontology graph (from child to parent). [Must be specified to propagate the counts]",
        default_value_t = false
    )]
    propagate_counts: bool,

    #[arg(
        long = "statistical-test",
        value_name = "TYPE",
        help = "Statistical test to use. [available: fisher, hypergeometric]",
        default_value = "fisher"
    )]
    statistical_test: String,

    #[arg(
        long = "min-protein-count",
        value_name = "MIN_COUNT",
        help = "Minimum protein count a GO Term must have to be processed.",
        default_value = "5",
        required = false
    )]
    min_protein_count: usize,
    
    #[arg(
        long = "min-score",
        value_name = "MIN_SCORE",
        help = "Minimum score (log(odds ratio)) a GO Term must have to be written in the results. Keeps GO terms with score ≥ threshold.",
        default_value = "0.5",
        required = false
    )]
    min_odds_ratio: f64,

    #[arg(
        long = "significance-threshold",
        value_name = "THRESHOLD",
        help = "P-value / Adjusted P-value threshold to determine significant results.",
        default_value = "0.05",
        required = false
    )]
    significance_threshold: f64,

    #[arg(
        long = "correction-method",
        value_name = "METHOD",
        help = "Method to adjust p-values for multiple test correction. [available: bonferroni, bh, by, none]",
        default_value = "none",
        required = false
    )]
    correction_method: String,

    #[arg(
        long = "combine-results",
        value_name = "TAXONOMIC_LEVEL",
        help = "Combine results from all taxa into a single output based on the taxonomic level specified. [Must be specified to group the results]",
    )]
    combine_results: Option<String>,
    
    #[arg(
        long = "lineage-percentage",
        value_name = "PERCENTAGE",
        help = "Percentage of species inside the desired taxonomic level in which the GO term must be found in. [From 0.0 to 1.0]",
        default_value = "0.25"
    )]
    lineage_percentage: f64,

    #[arg(
        long = "pm-iterations",
        value_name = "ITERATIONS_NUMBER",
        help = "Number of maximum iterations the Paule-Mandel estimator can reach when calculating the τ² estimate.",
        default_value = "1000",
        required = false
    )]
    pm_iterations: usize,

    #[arg(
        long = "pm-tolerance",
        value_name = "TOLERANCE",
        help = "Minimum acceptable tolerance between two τ² estimates.",
        default_value = "1e-6",
        required = false
    )]
    pm_tolerance: f64,
}

fn main() -> Result<(), Box<dyn Error>> {

    let cli_args: CliArgs = CliArgs::parse();  
    create_dir_all(&cli_args.output_dir)?;
    
    println!("\nReading ontology information from: {}\n\nBuilding ontology graph\n", &cli_args.obo_file);

    let ontology = parse_obo_file(&cli_args.obo_file)?;
    let (ontology_graph, node_index) = build_ontology_graph(&ontology)?;
    
    println!("Reading study populations from: {}\n", &cli_args.study_dir);

    let mut study_population  = read_study_pop(&cli_args.study_dir)?;
    
    let taxon_ids = study_population.taxon_map.keys().copied().collect::<HashSet<u32>>();

    println!("Reading background populations from: {}\n", &cli_args.background_dir);
    
    let mut background_data: BackgroundPop = read_background_pop(
        &taxon_ids, 
        &cli_args.background_dir);

    study_population.map_go_terms(&background_data.protein_to_go)?;

    if cli_args.propagate_counts{
        println!("Propagating counts up the Ontology graph\n");
        study_population.propagate_counts(&ontology_graph, &node_index);
        background_data.propagate_counts(&ontology_graph, &node_index);   
    }

    println!("Performing Gene Ontology (GO) term enrichment analysis\n");
    let analysis = EnrichmentAnalysis::new(
        cli_args.min_protein_count,
        match cli_args.statistical_test.to_lowercase().as_str() {
            "fisher" => StatisticalTest::Fishers,
            "hypergeometric" => StatisticalTest::Hypergeometric,
            _ => return Err("Invalid statistical test type".into()),
        }
    );
    
    let enrichment_results = analysis.analyze(
        &taxon_ids,          
        &background_data.go_term_count,
        &study_population.go_term_count,
        &background_data.taxon_protein_count,
        &study_population.taxon_protein_count,
    );
    
    let significant_fishers_results = adjust_species_p_values(&enrichment_results, &cli_args.correction_method, Some(cli_args.significance_threshold));
    
    write_single_taxon_results(&significant_fishers_results, &ontology, cli_args.min_odds_ratio, &cli_args.output_dir)?;
    
    if let Some(level_to_combine) = &cli_args.combine_results {

        println!("Grouping species based on {}\n", level_to_combine);

        println!("Reading taxonomic lineage information from: {}\n", DEFAULT_LINEAGE.to_string());
        let lineage = read_lineage(DEFAULT_LINEAGE.to_string());
        
        let grouped_species = taxid_to_level(
            &enrichment_results,
            &lineage?,
            level_to_combine
        );

        let lineage_organized_results = group_results_by_taxonomy(&grouped_species, &enrichment_results, cli_args.lineage_percentage);
        
        println!("Applying Paule-Mandel estimator with {} tolerance and {} iterations \n", &cli_args.pm_tolerance, &cli_args.pm_iterations);

        let complex_map = combine_taxonomic_results(&lineage_organized_results, cli_args.pm_tolerance, cli_args.pm_iterations);

        let significant_taxonomy_results = adjust_taxonomy_p_values(&complex_map, &cli_args.correction_method, Some(cli_args.significance_threshold), level_to_combine);
        
        write_taxonomy_results(&significant_taxonomy_results, &ontology, cli_args.min_odds_ratio, &cli_args.output_dir, level_to_combine)?;
    }

    println!("Finished analysis!\n");
    Ok(())
}