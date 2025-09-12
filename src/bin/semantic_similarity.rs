use clap::Parser;
use rustc_hash::{FxHashMap, FxHashSet};
use std::env::var;
use std::fs::create_dir_all;
use std::path::PathBuf;
use dirs::home_dir;
use daggy::NodeIndex;
use std::error::Error;
use petgraph::algo::toposort;

use TaxaGO::parsers::obo_parser::*;
use TaxaGO::parsers::background_parser::*;
use TaxaGO::utils::semantic_similarity::*;
use TaxaGO::utils::common_ancestor::*;
use TaxaGO::analysis::count_propagation::*;

fn get_default_asset_path(filename: &str) -> String {
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
        .join(filename)
        .to_string_lossy()
        .into_owned()
}

#[derive(Parser, Debug)]
#[command(name = "semantic-similarity")]
struct CliArgs {
    #[arg(
        short = 'o',
        long = "obo",
        value_name = "OBO_FILE",
        help = "Path to the Gene Ontology file in OBO format.",
    )]
    obo_file: Option<String>,
    
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
    )]
    background_dir: Option<String>,

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
        value_enum,
        help = "Method to calculate semantic similarity between two GO terms.",
        default_value_t = Method::Resnik,
    )]
    method: Method,

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
    
    let default_obo_path = get_default_asset_path("go.obo");
    let default_background_path = get_default_asset_path("background_pop");
    
    let obo_file = cli_args.obo_file.unwrap_or(default_obo_path);
    let background_dir = cli_args.background_dir.unwrap_or(default_background_path);
    
    create_dir_all(&cli_args.output_dir)?;
    
    println!("\nReading ontology information from: {}\n\nBuilding ontology graph\n", &obo_file);

    let obo_file_path = PathBuf::from(&obo_file);
    let ontology = parse_obo_file(&obo_file_path)?;
    let (ontology_graph, go_id_to_node_index) = build_ontology_graph(&ontology)?;
    
    let node_index_to_go_id: FxHashMap<NodeIndex, u32> = go_id_to_node_index
        .iter()
        .map(|(&go_id, &node_idx)| (node_idx, go_id))
        .collect();

    let topo_result = toposort(&ontology_graph, None).unwrap();

    let mut propagation_order: Vec<GOTermID> = topo_result
        .iter()
        .filter_map(|&node_idx| node_index_to_go_id.get(&node_idx).copied())
        .collect();
    
    propagation_order.reverse();
    
    println!("Reading background populations from: {}\n", &background_dir);
    
    let taxon_ids: FxHashSet<TaxonID> = cli_args.taxon_ids.split(',')
                                                        .map(|s| s.trim())
                                                        .filter_map(|s| s.parse::<u32>().ok())
                                                        .collect();
    let categories: Vec<EvidenceCategory> = vec![
        EvidenceCategory::Experimental,
        EvidenceCategory::Phylogenetic,
        EvidenceCategory::Computational,
        EvidenceCategory::Author,
        EvidenceCategory::Curator,
        EvidenceCategory::Electronic];

    let mut background_population = match BackgroundPop::read_background_pop(
        &taxon_ids, 
        &background_dir,
        &categories
    )? {
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

    if cli_args.propagate_counts {
        println!("Propagating counts up the Ontology graph\n");
        
        let ancestor_cache: GOAncestorCache = GOAncestorCache::new(
            &ontology_graph, 
            &ontology, 
            &go_id_to_node_index,
            &node_index_to_go_id
        )?;

        background_population.propagate_counts(&taxon_ids, &ancestor_cache);
    }

    let go_term_count: FxHashMap<u32, FxHashMap<u32, usize>> = background_population.go_term_count;

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
    
    for &taxon_id in &taxon_ids {        
        println!("Processing for Taxon ID: {}\n", taxon_id);
    
        let term_pairs = generate_term_pairs(
            &go_terms,
            taxon_id,
            &ic_results,
            &ontology_graph,
            &go_id_to_node_index,
            &node_index_to_go_id,
            &propagation_order,
            cli_args.method
        );
        
        if term_pairs.is_empty() && !go_terms.is_empty() {
            println!("No term pairs generated for taxon {}. This might happen if terms are not found or IC values are missing for IC-based methods.", taxon_id);
        }

        write_similarity_to_tsv(
            &term_pairs, 
            &go_terms, 
            taxon_id, 
            cli_args.method, 
            &cli_args.output_dir
        )
        .map_err(|e| format!("Failed to write similarity TSV for taxon {}: {}", taxon_id, e))?;
    }
    
    println!("All semantic similarity calculations completed successfully!\n");
    
    Ok(())
}
