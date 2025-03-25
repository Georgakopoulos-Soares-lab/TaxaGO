use clap::Parser;
use std::error::Error;
use std::fs::create_dir_all;
use std::env::var;
use std::collections::HashSet;
use TaxaGO::parsers::{
    background_parser::*, obo_parser::*, study_parser::*
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
    count_propagation::*

};
use std::process::Command;

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
// MUST ADD IT TO TAXAGO ASSETS !!!!!!
lazy_static! {
    static ref ENRICHMENT_PLOTS_SCRIPT: String = {
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
            .join("enrichment_plots.py")
            .to_string_lossy()
            .into_owned()
    };
}

#[derive(Parser, Debug)]
#[command(name = "taxago")]
struct CliArgs {
    #[arg(
        short='o',
        long = "obo",
        value_name = "OBO_FILE",
        help = "Path to the Gene Ontology file in OBO format.",
        default_value_t = DEFAULT_OBO_PATH.to_string(),
    
    )]
    obo_file: String,
    
    #[arg(
        short='s',
        long = "study",
        value_name = "STUDY_POP",
        help = "Directory containing study popultion for each taxon in FASTA format or CSV file with the study population for each species.",
        required = true
    )]
    study_pop: String,
    
    #[arg(
        short='b',
        long = "background",
        value_name = "BACKGROUND_DIR",
        help = "Directory containing background populations.",
        default_value_t = DEFAULT_BACKGROUND.to_string(),
    )]
    background_dir: String,
    
    #[arg(
        short='d',
        long = "dir",
        value_name = "RESULTS_DIR",
        help = "Directory to write results for each taxon and the combined results for the taxonomic level [if specified].",
        required = true
    )]
    output_dir: String,
    
    #[arg(
        short='p',
        long = "propagate-counts",
        help = "Propagates GO term counts upwards the Ontology graph (from child to parent). [Must be specified to propagate the counts]",
        default_value_t = false
    )]
    propagate_counts: bool,

    #[arg(
        short='t',
        long = "test",
        value_name = "TYPE",
        help = "Statistical test to use. [available: fisher, hypergeometric]",
        default_value = "fisher"
    )]
    statistical_test: String,

    #[arg(
        short='m',
        long = "min-prot",
        value_name = "MIN_COUNT",
        help = "Minimum protein count a GO Term must have to be processed.",
        default_value = "5",
        required = false
    )]
    min_protein_count: usize,
    
    #[arg(
        short='r',
        long = "min-score",
        value_name = "MIN_SCORE",
        help = "Minimum score (log(odds ratio)) a GO Term must have to be written in the results. Keeps GO terms with score ≥ threshold.",
        default_value = "0.5",
        required = false
    )]
    min_odds_ratio: f64,

    #[arg(
        short='a',
        long = "alpha",
        value_name = "THRESHOLD",
        help = "P-value / Adjusted P-value threshold to determine significant results.",
        default_value = "0.05",
        required = false
    )]
    significance_threshold: f64,

    #[arg(
        short='c',
        long = "correction-method",
        value_name = "METHOD",
        help = "Method to adjust p-values for multiple test correction. [available: bonferroni, bh, by, none]",
        default_value = "none",
        required = false
    )]
    correction_method: String,

    #[arg(
        short='g',
        long = "group-results",
        value_name = "TAXONOMIC_LEVEL",
        help = "Combine results from all taxa into a single output based on the taxonomic level specified. [Must be specified to group the results]",
    )]
    combine_results: Option<String>,
    
    #[arg(
        short = 'l',
        long = "lineage-percentage",
        value_name = "PERCENTAGE",
        help = "Percentage of species inside the desired taxonomic level in which the GO term must be found in. [From 0.0 to 1.0]",
        default_value = "0.25"
    )]
    lineage_percentage: f64,

    #[arg(
        short = 'i',
        long = "pm-iterations",
        value_name = "ITERATIONS_NUMBER",
        help = "Number of maximum iterations the Paule-Mandel estimator can reach when calculating the τ² estimate.",
        default_value = "1000",
        required = false
    )]
    pm_iterations: usize,

    #[arg(
        short = 'e',
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
    println!("Cleaning previous results");

    clean_directory(&cli_args.output_dir)?;

    create_dir_all(&cli_args.output_dir)?;

    println!("\nReading ontology information from: {}\n\nBuilding ontology graph\n", &cli_args.obo_file);

    let ontology = parse_obo_file(&cli_args.obo_file)?;
    let (ontology_graph, go_id_to_node_index) = build_ontology_graph(&ontology)?;
    
    println!("Reading background populations from: {}\n", &cli_args.background_dir);

    let taxon_ids: HashSet<TaxonID> = collect_taxon_ids(&cli_args.study_pop)?;
    
    let mut background_population = match BackgroundPop::read_background_pop(&taxon_ids, &cli_args.background_dir)? {
        Some(background_pop) => {
            println!("Successfully loaded background population for {} taxa\n", &taxon_ids.len());
            background_pop
        },
        None => {
            return Err(Box::new(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "No background population data could be loaded\n"
            )));
        }
    };

    println!("Reading study populations from: {}\n", &cli_args.study_pop);
    

    let mut study_population = match StudyPop::read_study_pop(&cli_args.study_pop, &background_population.protein_to_go)? {
        Some(study_pop) => {
            println!("Successfully loaded study population for {} taxa\n", &taxon_ids.len());
            study_pop
        },
        None => {
            return Err(Box::new(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "No study population data could be loaded\n"
            )));
        }
    };

    if cli_args.propagate_counts{
        println!("Propagating counts up the Ontology graph\n");
        
        let ancestor_cache: GOAncestorCache = GOAncestorCache::new(
            &ontology_graph, 
            &ontology, 
            &go_id_to_node_index)?;

        study_population.propagate_counts(&taxon_ids, &ancestor_cache);
        background_population.propagate_counts(&taxon_ids, &ancestor_cache);
        
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
        &background_population.go_term_count,
        &study_population.go_term_count,
        &background_population.taxon_protein_count,
        &study_population.taxon_protein_count,
    );
    
    let significant_fishers_results = adjust_species_p_values(
        &enrichment_results, 
        &cli_args.correction_method, 
        Some(cli_args.significance_threshold));
        
    let taxid_species_map = taxid_to_species(DEFAULT_LINEAGE.to_string())?;
    write_single_taxon_results(&significant_fishers_results, &ontology, cli_args.min_odds_ratio, &taxid_species_map, &cli_args.output_dir)?;
    
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
    println!("Finished analysis\n");
    
    println!("Generating enrichment plots\n");
    
    let plots_result = Command::new("python") 
        .arg(&*ENRICHMENT_PLOTS_SCRIPT)
        .arg("-r")
        .arg(&cli_args.output_dir)
        .arg("-s")
        .arg(&cli_args.study_pop)
        .status();

    match plots_result {
        Ok(status) => {
            if status.success() {
                println!("Enrichment plots created successfully\n");
            } else {
                eprintln!("Plot generation failed with exit code: {:?}\n", status.code().unwrap());
            }
        },
        Err(e) => {
            eprintln!("Failed to execute Python script: {}\n", e);
        }
    }
    Ok(())
}