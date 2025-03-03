use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::Path;
use std::io::{BufRead, BufReader};
use daggy::NodeIndex;

use crate::parsers::background_parser::*;

pub type InformationContent = f64;

pub fn parse_single_go_term(term: &str) -> Result<u32, String> {
    let term = term.trim();
    
    if !term.to_lowercase().starts_with("go:") {
        return Err(format!("Invalid GO term format: {}. Term must start with 'GO:'", term));
    }
    
    term[3..].parse::<u32>().map_err(|_| {
        format!("Invalid GO term number: {}. Expected a valid number after 'GO:'", &term[3..])
    })
}

pub fn parse_go_terms_from_iter<'a, I>(terms: I) -> Result<HashSet<u32>, String> 
where 
    I: Iterator<Item = &'a str>
{
    let mut results = HashSet::new();
    
    for term in terms {
        let term = term.trim();
        if !term.is_empty() {
            results.insert(parse_single_go_term(term)?);
        }
    }
    
    if results.is_empty() {
        return Err("No valid GO terms found".to_string());
    }
    
    Ok(results)
}

pub fn parse_go_terms(terms: &str) -> Result<HashSet<u32>, String> {
    parse_go_terms_from_iter(terms.split(','))
}

pub fn read_go_terms_from_file(file_path: &str) -> Result<HashSet<u32>, String> {
    let file = File::open(file_path)
        .map_err(|e| format!("Failed to open terms file: {}", e))?;
    
    let reader = BufReader::new(file);
    let lines = reader.lines()
        .map(|line_result| line_result.map_err(|e| format!("Error reading terms file: {}", e)))
        .collect::<Result<Vec<String>, String>>()?;
    
    parse_go_terms_from_iter(lines.iter().map(|s| s.as_str()))
        .map_err(|e| format!("{} in file: {}", e, file_path))
}

pub fn process_go_terms_input(input: &str) -> Result<HashSet<u32>, String> {
    if !input.contains(',') && Path::new(input).exists() {
        println!("Input appears to be a file path. Reading GO terms from file: {}\n", input);
        read_go_terms_from_file(input)
    } else {
        println!("Processing input as comma-separated GO terms\n");
        parse_go_terms(input)
    }
}

pub fn calculate_information_content(
    background_go_term_counts: &HashMap<u32, HashMap<u32, usize>>,
    go_terms: &HashSet<u32>,
    go_id_to_node_index: &HashMap<u32, NodeIndex>,
) -> HashMap<TaxonID, HashMap<u32, InformationContent>> {
    background_go_term_counts
        .iter()
        .map(|(&taxon_id, counts)| {
            let total_annotations: usize = counts.values().copied().sum();
            
            let ic_map: HashMap<u32, InformationContent> = go_terms
                .iter()
                .filter_map(|&go_id| {
                    if let Some(&count) = counts.get(&go_id) {
                        if count > 0 && total_annotations > 0 {
                            let probability = count as f64 / total_annotations as f64;
                            let ic = -probability.ln();
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