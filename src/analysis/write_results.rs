use rustc_hash::FxHashMap;
use std::fs::{self, File, create_dir_all};
use std::io::{self, BufWriter, Write};
use std::error::Error;
use std::path::PathBuf;
use std::fmt::Write as FmtWrite;
use lazy_static::lazy_static;

use crate::parsers::{
    background_parser::*,
    obo_parser::*
};
use crate::analysis::{
    enrichment_analysis::*,
    phylogenetic_meta_analysis::*
};

pub fn clean_directory(dir_path: &PathBuf) -> io::Result<()> {
    if dir_path.exists() {
        for entry in fs::read_dir(dir_path)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_file() {
                fs::remove_file(path)?;
            } else if path.is_dir() {
                fs::remove_dir_all(path)?;
            }
        }
    }
    Ok(())
}

const BUFFER_SIZE: usize = 8192 * 32;

lazy_static! {
    static ref NAMESPACE_MAPPING: FxHashMap<&'static str, &'static str> = {
        let mut go_term_class_map = FxHashMap::default();
        go_term_class_map.insert("biological_process", "Biological Process");
        go_term_class_map.insert("molecular_function", "Molecular Function");
        go_term_class_map.insert("cellular_component", "Cellular Component");
        
        go_term_class_map
    };
}

struct TermCache {
    go_terms: FxHashMap<u32, String>,
}

impl TermCache {
    fn new() -> Self {
        Self {
            go_terms: FxHashMap::with_capacity_and_hasher(
                10000,
                rustc_hash::FxBuildHasher::default()
            ),
        }
    }

    #[inline]
    fn get_go_term<'a>(&'a mut self, go_id: u32) -> &'a str {
        self.go_terms.entry(go_id).or_insert_with(|| {
            format!("GO:{:07}", go_id)
        })
    }
}

#[inline]
fn format_namespace(namespace: &str) -> &str {
    NAMESPACE_MAPPING.get(namespace).unwrap_or(&namespace)
}

pub fn write_single_taxon_results(
    data: &FxHashMap<u32, FxHashMap<GOTermID, GOTermResults>>,
    ontology: &FxHashMap<u32, OboTerm>,
    min_log_odds_ratio: f64,
    taxid_species_map: &FxHashMap<TaxonID, String>,
    output_dir: &PathBuf,
) -> Result<(), Box<dyn Error>> {
    let results_dir = PathBuf::from(output_dir).join("single_taxon_results");
    let plots_dir = PathBuf::from(&results_dir).join("plots"); 
    create_dir_all(&results_dir)?;
    create_dir_all(&plots_dir)?;

    println!("Writing single taxon results to: {}\n", results_dir.to_str().unwrap());
    
    let mut term_cache = TermCache::new();
    
    let mut line_buffer = String::with_capacity(256);
    
    for (taxon_id, go_terms) in data {
        let species_name = taxid_species_map.get(taxon_id)
            .unwrap_or(&taxon_id.to_string())
            .replace(" ", "_");

        let filename = results_dir.join(format!("{}_GOEA_results.txt", species_name));
        let file = File::create(&filename)?;
        let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
        
        writer.write_all(b"GO Term ID\tName\tNamespace\tlog(Odds Ratio)\tStatistical significance\n")?;
        // \tN Study with term\tN Study without term\tN Background with term\tN Background without term
        for (go_term, results) in go_terms {
            if results.log_odds_ratio >= min_log_odds_ratio {
                if let Some(term) = ontology.get(go_term) {
                    if !term.is_obsolete {
                        let formatted_go_term = term_cache.get_go_term(*go_term);
                        let namespace_str = match term.namespace {
                            NameSpace::BiologicalProcess => "biological_process",
                            NameSpace::MolecularFunction => "molecular_function",
                            NameSpace::CellularComponent => "cellular_component",
                        };
                        let formatted_namespace = format_namespace(namespace_str);
                        
                        line_buffer.clear();
                        
                        write!(
                            &mut line_buffer,
                            "{}\t{}\t{}\t{:.3}\t{:.5e}\n",
                            formatted_go_term,
                            term.name,
                            formatted_namespace,
                            results.log_odds_ratio,
                            results.p_value,
                        )?;
                        
                        writer.write_all(line_buffer.as_bytes())?;
                    }
                }
            }
        }
        
        writer.flush()?;
    }
    Ok(())
}

pub fn write_taxonomy_results(
    data: &FxHashMap<String, FxHashMap<u32, TaxonomyGOResult>>,
    ontology: &FxHashMap<u32, OboTerm>,
    min_log_odds_ratio: f64,
    output_dir: &PathBuf,
    level: &String,
) -> Result<(), Box<dyn Error>> {
    let results_dir = PathBuf::from(output_dir).join("combined_taxonomy_results");
    let plots_dir = PathBuf::from(&results_dir).join("plots");
    create_dir_all(&results_dir)?;
    create_dir_all(&plots_dir)?;

    println!("Writing {} results to: {}\n", level, results_dir.to_str().unwrap());
    
    let mut term_cache = TermCache::new();
    let mut line_buffer = String::with_capacity(256);
    
    for (taxonomy, go_terms) in data {
        let filename = results_dir.join(format!("{}_GOEA_results.txt", taxonomy));
        let file = File::create(&filename)?;
        let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
        
        writer.write_all(b"GO Term ID\tName\tNamespace\tlog(Odds Ratio)\tStatistical significance\n")?;
        // \tHeterogeneity\tSpecies Percentage\tN with GO term\tN in taxonomy
        for (go_term, result) in go_terms {
            if result.log_odds_ratio >= min_log_odds_ratio {
                if let Some(term) = ontology.get(go_term) {
                    if !term.is_obsolete {
                        let formatted_go_term = term_cache.get_go_term(*go_term);
                        let namespace_str = match term.namespace {
                            NameSpace::BiologicalProcess => "biological_process",
                            NameSpace::MolecularFunction => "molecular_function",
                            NameSpace::CellularComponent => "cellular_component",
                        };
                        let formatted_namespace = format_namespace(namespace_str);
                        
                        line_buffer.clear();
                        write!(
                            &mut line_buffer,
                            "{}\t{}\t{}\t{:.3}\t{:.5e}\n",
                            formatted_go_term,
                            term.name,
                            formatted_namespace,
                            result.log_odds_ratio,
                            result.p_value
                        )?;
                        
                        writer.write_all(line_buffer.as_bytes())?;
                    }
                }
            }
        }
        
        writer.flush()?;
    }
    Ok(())
}