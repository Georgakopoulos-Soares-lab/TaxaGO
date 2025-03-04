use clap::Parser;
use std::collections::{HashMap, HashSet};
use std::env::var;
use lazy_static::lazy_static;
use std::fs::create_dir_all;
use std::path::PathBuf;
use dirs::home_dir;
use daggy::NodeIndex;
use std::error::Error;

use TaxaGO::parsers::obo_parser::*;
use TaxaGO::parsers::background_parser::*;
use TaxaGO::utils::semantic_similarity::*;
use TaxaGO::utils::common_ancestor::*;

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

#[derive(Parser, Debug)]
#[command(name = "semantic-similarity")]
struct CliArgs {
    #[arg(
        short = 'o',
        long = "obo",
        value_name = "OBO_FILE",
        help = "Path to the Gene Ontology file in OBO format.",
        default_value_t = DEFAULT_OBO_PATH.to_string(),
    
    )]
    obo_file: String,
    
    #[arg(
        short='t',
        long = "terms",
        value_name = "GO_TERMS_OR_FILE",
        help = "Either a comma-separated list of GO terms [e.g., GO:0016070,GO:0140187] or a path to a file containing GO terms (one per line)",
    )]
    go_terms_input: String,

    #[arg(
        short ='i',
        long = "ids",
        value_name = "TAXON_IDS",
        help = "Comma-separated list of Taxon IDs [e.g., 9606,10090]",
        default_value = "9606",

    )]
    taxon_ids: String,
    
    #[arg(
        short = 'b',
        long = "background",
        value_name = "BACKGROUND_DIR",
        help = "Directory containing background populations.",
        default_value_t = DEFAULT_BACKGROUND.to_string(),
    )]
    background_dir: String,

    #[arg(
        short = 'd',
        long = "dir",
        value_name = "RESULTS_DIR",
        help = "Directory to write results.",
        default_value = "./",

    )]
    output_dir: String,
    
    #[arg(
        short = 'm',
        long = "method",
        value_name = "METHOD",
        help = "Method to calculate semantic similarity between two GO terms. [available: resnik, lin, jiang-conrath, wang]",
        default_value = "resnik",
    )]
    method: String,

    #[arg(
        short = 'p',
        long = "propagate-counts",
        help = "Propagates GO term counts upwards the Ontology graph (from child to parent). [Must be specified to propagate the counts]",
        default_value_t = false
    )]
    propagate_counts: bool,
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli_args: CliArgs = CliArgs::parse();  
    create_dir_all(&cli_args.output_dir)?;
    
    println!("\nReading ontology information from: {}\n\nBuilding ontology graph\n", &cli_args.obo_file);

    let ontology = parse_obo_file(&cli_args.obo_file)?;
    let (ontology_graph, go_id_to_node_index) = build_ontology_graph(&ontology)?;
    
    let node_index_to_go_id: HashMap<NodeIndex, u32> = go_id_to_node_index
        .iter()
        .map(|(&go_id, &node_idx)| (node_idx, go_id))
        .collect();
    
    println!("Reading background populations from: {}\n", &cli_args.background_dir);
    
    let taxon_ids: HashSet<TaxonID> = cli_args.taxon_ids.split(',')
                                                        .map(|s| s.trim())
                                                        .filter_map(|s| s.parse::<u32>().ok())
                                                        .collect();

    let mut background_data: BackgroundPop = read_background_pop(
        &taxon_ids, 
        &cli_args.background_dir);

    if cli_args.propagate_counts {
        println!("Propagating counts up the Ontology graph\n");
        background_data.propagate_counts(&ontology_graph, &go_id_to_node_index);   
    }
    let go_term_count: HashMap<u32, HashMap<u32, usize>> = background_data.go_term_count;

    let go_terms = process_go_terms_input(&cli_args.go_terms_input)?;

    println!("Calculating Information Content (IC) for {} GO terms\n", go_terms.len());

    let mut expanded_terms = go_terms.clone();
    for &term in go_terms.iter() {
        if let Some(&node_idx) = go_id_to_node_index.get(&term) {
            let ancestry_path = collect_ancestry_path(&ontology_graph, node_idx);
            for (idx, _) in ancestry_path {
                let ancestor_id = node_index_to_go_id[&idx];
                expanded_terms.insert(ancestor_id);
            }
        }
    }
    let ic_results = calculate_information_content(
        &go_term_count,
        &expanded_terms,
        &go_id_to_node_index
    );

    println!("Finding Most Informative Common Ancestor (MICA)\n");
    
    println!("Calculating semantic similarity using {} method \n", cli_args.method);
    
    for &taxon_id in &taxon_ids {        
        let term_pairs = generate_term_pairs(
            &go_terms,
            taxon_id,
            &ic_results,
            &ontology_graph,
            &go_id_to_node_index,
            &node_index_to_go_id,
            &cli_args.method
        );
        write_similarity_to_tsv(&term_pairs, &go_terms, taxon_id, &cli_args.output_dir);
    }
    
    println!("All semantic similarity calculations completed successfully!\n");
    
    Ok(())
}