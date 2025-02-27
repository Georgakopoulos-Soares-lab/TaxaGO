use clap::Parser;
use std::collections::{HashMap,HashSet};
use std::env::var;
use lazy_static::lazy_static;
use std::fs::{File, create_dir_all};
use std::path::{Path,PathBuf};
use dirs::home_dir;
use TaxaGO::parsers::{
    obo_parser::*,
    background_parser::*
};
use std::error::Error;
use std::io::{self, BufRead, BufReader};
use rayon::prelude::*;
use daggy::{NodeIndex, Walker};

pub type InformationContent = f64;

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

    )]
    output_dir: String,
    
    #[arg(
        short = 'm',
        long = "method",
        value_name = "METHOD",
        help = "Method to calculate semantic similarity between two GO terms. [available: resnik, wang, lin]",
        default_value = "wang",
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

#[derive(Debug, Default)]
pub struct BackgroundCounts {
    pub go_term_counts: HashMap<TaxonID, GOTermCount>,
}

impl BackgroundCounts {
    pub fn new() -> Self {
        Self {
            go_term_counts: HashMap::new(),
        }
    }

    pub fn propagate_counts(
        &mut self,
        graph: &OntologyGraph,
        go_id_to_node_index: &HashMap<u32, NodeIndex>
    ) {
        let ancestor_cache: HashMap<NodeIndex, HashSet<NodeIndex>> = go_id_to_node_index
            .values()
            .par_bridge()
            .map(|&node_idx| (node_idx, get_unique_ancestors(node_idx, graph)))
            .collect();

        let node_index_to_go_id: HashMap<NodeIndex, u32> = go_id_to_node_index
            .iter()
            .map(|(&go_id, &node_idx)| (node_idx, go_id))
            .collect();

        let updated_counts: HashMap<_, _> = self.go_term_counts
            .par_iter()
            .map(|(&taxon_id, go_term_counts)| {
                let mut propagated_counts = HashMap::with_capacity(go_term_counts.len() * 2);
                
                for (&go_id, &direct_count) in go_term_counts.iter() {
                    if let Some(&node_idx) = go_id_to_node_index.get(&go_id) {
                        *propagated_counts.entry(go_id).or_insert(0) += direct_count;
                        
                        if let Some(ancestors) = ancestor_cache.get(&node_idx) {
                            for &ancestor_idx in ancestors {
                                if let Some(&ancestor_id) = node_index_to_go_id.get(&ancestor_idx) {
                                    *propagated_counts.entry(ancestor_id).or_insert(0) += direct_count;
                                }
                            }
                        }
                    }
                }
                
                (taxon_id, propagated_counts)
            })
            .collect();
        
        self.go_term_counts = updated_counts;
    }
}

fn process_taxon_background(file_path: impl AsRef<Path>) -> io::Result<GOTermCount> {
    let file = File::open(file_path)?;
    let reader = BufReader::with_capacity(128 * 1024, file);
    let mut go_term_counts = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        
        if parts.len() < 2 {
            continue;
        }

        if let Some(go_str) = parts[1].strip_prefix("GO:") {
            if let Ok(go_id) = go_str.parse::<u32>() {
                *go_term_counts.entry(go_id).or_insert(0) += 1;
            }
        }
    }

    Ok(go_term_counts)
}

fn get_unique_ancestors(node_idx: NodeIndex, graph: &OntologyGraph) -> HashSet<NodeIndex> {
    let mut ancestors = HashSet::new();
    let mut to_visit = vec![node_idx];
    let mut visited = HashSet::new();
    
    while let Some(current_idx) = to_visit.pop() {
        if !visited.insert(current_idx) {
            continue;
        }
        
        let mut parents = graph.parents(current_idx);
        while let Some((edge_idx, parent_idx)) = parents.walk_next(graph) {
            match graph.edge_weight(edge_idx).unwrap() {
                Relationship::IsA | Relationship::PartOf => {
                    ancestors.insert(parent_idx);
                    to_visit.push(parent_idx);
                },
                _ => continue,
            }
        }
    }
    
    ancestors
}

pub fn load_background(
    background_path: String,
    taxon_ids: String,
) -> io::Result<BackgroundCounts> {
    let taxon_id_set: HashSet<u32> = taxon_ids
        .split(',')
        .map(|s| s.trim())
        .filter_map(|s| s.parse().ok())
        .collect();

    let mut background_counts = BackgroundCounts::new();

    let results: Vec<(TaxonID, io::Result<GOTermCount>)> = taxon_id_set
        .par_iter()
        .map(|&taxon_id| {
            let file_path = Path::new(&background_path)
                .join(format!("{}_background.txt", taxon_id));
            
            (taxon_id, process_taxon_background(&file_path))
        })
        .collect();

    for (taxon_id, result) in results {
        match result {
            Ok(counts) => {
                background_counts.go_term_counts.insert(taxon_id, counts);
            },
            Err(e) => {
                eprintln!("Error processing taxon {}: {}", taxon_id, e);
            }
        }
    }

    Ok(background_counts)
}

fn parse_single_go_term(term: &str) -> Result<u32, String> {
    let term = term.trim();
    
    let numeric_part = if term.to_lowercase().starts_with("go:") {
        &term[3..]
    } else {
        return Err(format!("Invalid GO term format: {}. Term must start with 'GO:'", term));
    };
    
    numeric_part.parse::<u32>().map_err(|_| {
        format!("Invalid GO term number: {}. Expected a valid number after 'GO:'", numeric_part)
    })
}

fn parse_go_terms(terms: &str) -> Result<Vec<u32>, String> {
    terms
        .split(',')
        .map(parse_single_go_term)
        .collect()
}

fn read_go_terms_from_file(file_path: &str) -> Result<Vec<u32>, String> {
    let file = File::open(file_path)
        .map_err(|e| format!("Failed to open terms file: {}", e))?;
    
    let reader = BufReader::new(file);
    let mut go_terms = Vec::new();
    
    for line_result in reader.lines() {
        let line = line_result
            .map_err(|e| format!("Error reading terms file: {}", e))?;
        
        let term = line.trim();
        if !term.is_empty() {
            match parse_single_go_term(term) {
                Ok(go_id) => go_terms.push(go_id),
                Err(e) => return Err(e),
            }
        }
    }
    
    if go_terms.is_empty() {
        return Err(format!("No valid GO terms found in file: {}", file_path));
    }
    
    Ok(go_terms)
}

fn process_go_terms_input(input: &str) -> Result<Vec<u32>, String> {
    if !input.contains(',') && Path::new(input).exists() {
        println!("Input appears to be a file path. Reading GO terms from file: {}\n", input);
        read_go_terms_from_file(input)
    } else {
        println!("Processing input as comma-separated GO terms\n");
        parse_go_terms(input)
    }
}

fn calculate_information_content(
    background_counts: &BackgroundCounts,
    go_terms: &[u32],
    go_id_to_node_index: &HashMap<u32, NodeIndex>,
) -> HashMap<TaxonID, HashMap<u32, InformationContent>> {
    background_counts.go_term_counts
        .par_iter()
        .map(|(&taxon_id, counts)| {

            let total_annotations: usize = counts.values().copied().sum();
            
            let ic_map: HashMap<u32, InformationContent> = go_terms
                .iter()
                .filter_map(|&go_id| {
                    if let Some(&count) = counts.get(&go_id) {
                        if count > 0 && total_annotations > 0 {
                            let probability = count as f64 / total_annotations as f64;
                            let ic = -probability.log2();
                            Some((go_id, ic))
                        } else {
                            if go_id_to_node_index.contains_key(&go_id) {
                                Some((go_id, f64::INFINITY))
                            } else {
                                None
                            }
                        }
                    } else {
                        None
                    }
                })
                .collect();
            
            (taxon_id, ic_map)
        })
        .collect()
}

fn main() -> Result<(), Box<dyn Error>> {
    let cli_args: CliArgs = CliArgs::parse();  
    create_dir_all(&cli_args.output_dir)?;
    
    println!("\nReading ontology information from: {}\n\nBuilding ontology graph\n", &cli_args.obo_file);

    let ontology = parse_obo_file(&cli_args.obo_file)?;
    let (ontology_graph, go_id_to_node_index) = build_ontology_graph(&ontology)?;
    
    println!("Reading background populations from: {}\n", &cli_args.background_dir);
    
    let mut background_counts = load_background(
        cli_args.background_dir,
        cli_args.taxon_ids
    )?;

    if cli_args.propagate_counts {
        println!("Propagating counts up the Ontology graph\n");
        background_counts.propagate_counts(&ontology_graph, &go_id_to_node_index);   
    }
    
    let go_terms = process_go_terms_input(&cli_args.go_terms_input)?;
    
    println!("Calculating information content for {} GO terms\n", go_terms.len());

    let ic_results = calculate_information_content(
        &background_counts,
        &go_terms,
        &go_id_to_node_index
    );
    println!("{:?}", &ic_results);
    Ok(())
}