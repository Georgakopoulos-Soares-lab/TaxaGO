use clap::Parser;
use rustc_hash::FxHashMap;
use std::error::Error;
use std::process::Command;
use lazy_static::lazy_static;
use std::path::PathBuf;
use dirs::home_dir;
use std::env::var;
use daggy::NodeIndex;

use TaxaGO::parsers::obo_parser::*;
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

#[derive(Parser, Debug)]
#[command(name = "common-ancestors")]
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
        short = 't',
        long = "terms",
        value_name = "GO_TERMS",
        help = "Comma-separated list of GO terms [e.g., GO:0016070,GO:0140187].",
        required = true,
    )]
    go_terms: String,

    #[arg(
        short = 'd',
        long = "dir",
        value_name = "RESULTS_DIR",
        help = "Directory to write results.",
        default_value="./"
    )]
    graph_path: String
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

fn main() -> Result<(), Box<dyn Error>> {
    let cli_args: CliArgs = CliArgs::parse();  
    
    let target_go_ids = match parse_go_terms(&cli_args.go_terms) {
        Ok(ids) => ids,
        Err(e) => {
            eprintln!("Error parsing GO terms: {}", e);
            std::process::exit(1);
        }
    };
    
    let obo_map = parse_obo_file(&cli_args.obo_file)?;
    let (ontology_graph, go_id_to_node_index) = build_ontology_graph(&obo_map)?;
    
    let node_index_to_go_id: FxHashMap<NodeIndex, u32> = go_id_to_node_index
        .iter()
        .map(|(&go_id, &node_idx)| (node_idx, go_id))
        .collect();

    for &go_id in &target_go_ids {
        if !go_id_to_node_index.contains_key(&go_id) {
            eprintln!("Error: GO:{:07} not found in the ontology", go_id);
            std::process::exit(1);
        }
    }
    
    let mut ancestry_paths = Vec::new();
    
    for &go_id in &target_go_ids {
        let target_node_idx = go_id_to_node_index[&go_id];        
        let path = collect_ancestry_path(&ontology_graph, target_node_idx);
        ancestry_paths.push(path);
    }
    
    let common_ancestors = find_common_ancestors(&ancestry_paths, &node_index_to_go_id);
    
    if common_ancestors.is_empty() {
        eprintln!("No common ancestors found between the provided GO terms");
        std::process::exit(1);
    }

    let first_common_ancestor = find_first_common_ancestor(
        &ancestry_paths, 
        &node_index_to_go_id,
        &ontology_graph
    );
    
    println!("\nAnalyzing GO terms: {}", target_go_ids.iter()
        .map(|id| format!("GO:{:07}", id))
        .collect::<Vec<_>>()
        .join(", "));
    
    println!("\nCommon ancestors:");
    for &go_id in &common_ancestors {
        let term = &obo_map[&go_id];
        println!("GO:{:07} - {}", go_id, term.name.replace("_", " "));
    }
    
    let mermaid_chart = generate_mermaid_chart(
        &ontology_graph,
        &target_go_ids,
        &go_id_to_node_index,
        &node_index_to_go_id,
        &obo_map,
        first_common_ancestor
    );

    let graph_file = cli_args.graph_path.to_string() + "ontology_graph.mmd";
    let graph_pdf = cli_args.graph_path.to_string() + "ontology_graph.pdf";

    std::fs::write(&graph_file, &mermaid_chart)?;
    println!("\nMermaid graph Markdown file has been written to: {}\n", &graph_file);

    let mermaid_cli = Command::new("mmdc")
    .arg("-i")
    .arg(&graph_file)
    .arg("-o")
    .arg(&graph_pdf)
    .arg("-f")
    .output()
    .expect("Failed to execute mmdc command\n");

    if mermaid_cli.status.success() {
        println!("Mermaid chart PDF file has been written to: {}\n", &graph_pdf);
    } else {
        eprintln!("Error: {}\n", String::from_utf8_lossy(&mermaid_cli.stderr));
    }
    Ok(())
}