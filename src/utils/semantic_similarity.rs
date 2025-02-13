use clap::Parser;
use std::env::var;
use lazy_static::lazy_static;
use std::path::PathBuf;
use dirs::home_dir;

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
        long = "terms",
        value_name = "GO_TERMS",
        help = "Comma-separated list of GO terms (e.g., GO:0016070,GO:0140187)",
        required = true
    )]
    go_terms: String,
    
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
        help = "Directory to write results.",
        required = true
    )]
    output_dir: String,
    
    #[arg(
        long = "method",
        value_name = "METHOD",
        help = "Method to calculate semantic similarity between two GO terms. [available: resnik, wang, lin]",
        default_value = "wang",
        required = true
    )]
    method: String,
}

fn main() {
    let _cli_args: CliArgs = CliArgs::parse();  
    todo!();
}