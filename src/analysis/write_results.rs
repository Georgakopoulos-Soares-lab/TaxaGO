use std::collections::HashMap;
use std::fs::{self, File, create_dir_all};
use std::io::{self, BufWriter, Write};
use std::error::Error;
use std::path::{PathBuf, Path};
use std::fmt::Write as FmtWrite;
use lazy_static::lazy_static;
use crate::parsers::obo_parser::{OboTerm, NameSpace};
use crate::analysis::enrichment_analysis::GOTermResults;
use crate::analysis::multiple_testing_correction::TaxonomyGOResult;

fn clean_directory(dir_path: &str) -> io::Result<()> {
    let path = Path::new(dir_path);
    if path.exists() {
        for entry in fs::read_dir(path)? {
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
    static ref NAMESPACE_MAPPING: HashMap<&'static str, &'static str> = {
        let mut go_term_class_map = HashMap::new();
        go_term_class_map.insert("biological_process", "Biological Process");
        go_term_class_map.insert("molecular_function", "Molecular Function");
        go_term_class_map.insert("cellular_component", "Cellular Component");
        
        go_term_class_map
    };
}

struct TermCache {
    go_terms: HashMap<u32, String>,
}

impl TermCache {
    fn new() -> Self {
        Self {
            go_terms: HashMap::new(),
        }
    }

    #[inline]
    fn get_go_term(&mut self, go_id: u32) -> String {
        if let Some(term) = self.go_terms.get(&go_id) {
            term.clone()
        } else {
            let term = format!("GO:{:07}", go_id);
            self.go_terms.insert(go_id, term.clone());
            term
        }
    }
}

#[inline]
fn format_namespace(namespace: &str) -> &str {
    NAMESPACE_MAPPING.get(namespace).unwrap_or(&namespace)
}

pub fn write_single_taxon_results(
    data: &HashMap<u32, HashMap<u32, GOTermResults>>,
    ontology: &HashMap<u32, OboTerm>,
    min_log_odds_ratio: f64,
    taxid_species_map: &HashMap<u32, String>,
    output_dir: &str
) -> Result<(), Box<dyn Error>> {
    let results_dir = PathBuf::from(output_dir).join("single_taxon_results");
    let plots_dir = PathBuf::from(&results_dir).join("plots");
    clean_directory(&output_dir)?; 
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
        
        writer.write_all(b"GO Term ID\tName\tNamespace\tlog(Odds Ratio)\tStatistical significance\tN Study with term\tN Study without term\tN Background with term\tN Background without term\n")?;
        
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
                            "{}\t{}\t{}\t{:.3}\t{:.5e}\t{}\t{}\t{}\t{}\n",
                            formatted_go_term,
                            term.name,
                            formatted_namespace,
                            results.log_odds_ratio,
                            results.p_value,
                            results.contingency_table[0],
                            results.contingency_table[1],
                            results.contingency_table[2],
                            results.contingency_table[3]
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
    data: &HashMap<String, HashMap<u32, TaxonomyGOResult>>,
    ontology: &HashMap<u32, OboTerm>,
    min_log_odds_ratio: f64,
    output_dir: &str,
    level: &String,
) -> Result<(), Box<dyn Error>> {
    let results_dir = PathBuf::from(output_dir).join("combined_taxonomy_results");
    let plots_dir = PathBuf::from(&results_dir).join("plots");
    clean_directory(&output_dir)?; 
    create_dir_all(&results_dir)?;
    create_dir_all(&plots_dir)?;

    println!("Writing {} results to: {}\n", level, results_dir.to_str().unwrap());
    
    let mut term_cache = TermCache::new();
    let mut line_buffer = String::with_capacity(256);
    
    for (taxonomy, go_terms) in data {
        let filename = results_dir.join(format!("{}_GOEA_results.txt", taxonomy));
        let file = File::create(&filename)?;
        let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
        
        writer.write_all(b"GO Term ID\tName\tNamespace\tlog(Odds Ratio)\tStatistical significance\tHeterogeneity\tSpecies Percentage\tN with GO term\tN in taxonomy\n")?;
        
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
                            "{}\t{}\t{}\t{:.3}\t{:.5e}\t{:.5e}\t{:.3}\t{}\t{}\n",
                            formatted_go_term,
                            term.name,
                            formatted_namespace,
                            result.log_odds_ratio,
                            result.p_value,
                            result.tau_squared,
                            result.species_percentage,
                            result.species_count,
                            result.total_species
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