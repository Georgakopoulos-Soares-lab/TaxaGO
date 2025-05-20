use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::{self,File};
use std::path::Path;
use std::io::{BufRead, BufReader, Write};
use daggy::NodeIndex;
use crate::parsers::{
    background_parser::*,
    obo_parser::*,
};
use crate::utils::common_ancestor::*;

pub type InformationContent = f64;
#[derive(Debug)]
pub struct TermPair {
    pub term1: u32,
    pub term2: u32,
    pub ic_term1: f64,
    pub ic_term2: f64,
    pub mica: (u32, f64),
    pub similarity: f64,
}

pub fn parse_single_go_term(term: &str) -> Result<u32, String> {
    let term = term.trim();
    
    if !term.to_lowercase().starts_with("go:") {
        return Err(format!("Invalid GO term format: {}. Term must start with 'GO:'", term));
    }
    
    term[3..].parse::<u32>().map_err(|_| {
        format!("Invalid GO term number: {}. Expected a valid number after 'GO:'", &term[3..])
    })
}

pub fn parse_go_terms_from_iter<'a, I>(terms: I) -> Result<FxHashSet<u32>, String> 
where 
    I: Iterator<Item = &'a str>
{
    let mut results = FxHashSet::default();
    
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

pub fn parse_go_terms(terms: &str) -> Result<FxHashSet<u32>, String> {
    parse_go_terms_from_iter(terms.split(','))
}

pub fn read_go_terms_from_file(file_path: &str) -> Result<FxHashSet<u32>, String> {
    let file = File::open(file_path)
        .map_err(|e| format!("Failed to open terms file: {}", e))?;
    
    let reader = BufReader::new(file);
    let lines = reader.lines()
        .map(|line_result| line_result.map_err(|e| format!("Error reading terms file: {}", e)))
        .collect::<Result<Vec<String>, String>>()?;
    
    parse_go_terms_from_iter(lines.iter().map(|s| s.as_str()))
        .map_err(|e| format!("{} in file: {}", e, file_path))
}

pub fn process_go_terms_input(input: &str) -> Result<FxHashSet<u32>, String> {
    if !input.contains(',') && Path::new(input).exists() {
        println!("Input appears to be a file path. Reading GO terms from file: {}\n", input);
        read_go_terms_from_file(input)
    } else {
        println!("Processing input as comma-separated GO terms\n");
        parse_go_terms(input)
    }
}

pub fn calculate_information_content(
    background_go_term_counts: &FxHashMap<u32, FxHashMap<u32, usize>>,
    go_terms: &FxHashSet<u32>,
    go_id_to_node_index: &FxHashMap<u32, NodeIndex>,
) -> FxHashMap<TaxonID, FxHashMap<u32, InformationContent>> {
    background_go_term_counts
        .iter()
        .map(|(&taxon_id, counts)| {
            let total_annotations: usize = counts.values().copied().sum();
            
            let ic_map: FxHashMap<u32, InformationContent> = go_terms
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

impl TermPair {
    pub fn resnik_similarity(&self) -> f64 {
        self.mica.1
    }
    
    pub fn lin_similarity(&self) -> f64 {
        if self.ic_term1 == 0.0 || self.ic_term2 == 0.0 {
            return 0.0;
        }
        
        (2.0 * self.mica.1) / (self.ic_term1 + self.ic_term2)
    }
    
    pub fn jiang_conrath_similarity(&self) -> f64 {
        let distance = self.ic_term1 + self.ic_term2 - 2.0 * self.mica.1;
        1.0 / (1.0 + distance.max(0.0))  
    }
    
    pub fn calculate_similarity(&mut self, method: &str) -> f64 {
        let similarity = match method.to_lowercase().as_str() {
            "resnik" => self.resnik_similarity(),
            "lin" => self.lin_similarity(),
            "jiang-conrath" => self.jiang_conrath_similarity(),
            _ => {
                println!("Warning: Unsupported method '{}', defaulting to Resnik", method);
                self.resnik_similarity()
            }
        };
        
        self.similarity = similarity;
        similarity
    }
    
    
    pub fn new(
        term1: u32, 
        term2: u32, 
        ic_term1: f64, 
        ic_term2: f64, 
        mica: (u32, f64), 
        method: &str) -> Self {
            let mut pair = Self {
                term1,
                term2,
                ic_term1,
                ic_term2,
                mica,
                similarity: 0.0,
            };
        
        pair.calculate_similarity(method);
        pair
    }
}

pub fn find_mica_for_pair(
    term1: u32,
    term2: u32,
    ontology_graph: &OntologyGraph,
    go_id_to_node_index: &FxHashMap<u32, NodeIndex>,
    node_index_to_go_id: &FxHashMap<NodeIndex, u32>,
    ic_values: &FxHashMap<u32, f64>,
) -> Option<(u32, f64)> {
    if term1 == term2 {
        if let Some(&ic) = ic_values.get(&term1) {
            return Some((term1, ic));
        }
    }
    
    let node_idx1 = match go_id_to_node_index.get(&term1) {
        Some(&idx) => idx,
        None => return None,
    };
    
    let node_idx2 = match go_id_to_node_index.get(&term2) {
        Some(&idx) => idx,
        None => return None,
    };
    
    let path1 = collect_ancestry_path(ontology_graph, node_idx1);
    let path2 = collect_ancestry_path(ontology_graph, node_idx2);
    
    let ancestry_paths = vec![path1, path2];
    let common_ancestors = find_common_ancestors(&ancestry_paths, node_index_to_go_id);

    let mut max_ic = f64::NEG_INFINITY;
    let mut mica_id = 0;
    
    for ancestor_id in common_ancestors {  
        if let Some(&ic) = ic_values.get(&ancestor_id) {
            if ic > max_ic {
                max_ic = ic;
                mica_id = ancestor_id;
            }
        }
    }
    
    if mica_id != 0 {
        Some((mica_id, max_ic))
    } else {
        None
    }
}

pub fn generate_term_pairs(
    go_terms: &FxHashSet<u32>,
    taxon_id: TaxonID,
    ic_results: &FxHashMap<TaxonID, FxHashMap<u32, f64>>,
    ontology_graph: &OntologyGraph,
    go_id_to_node_index: &FxHashMap<u32, NodeIndex>,
    node_index_to_go_id: &FxHashMap<NodeIndex, u32>,
    method: &str,
) -> Vec<TermPair> {
    let terms: Vec<u32> = go_terms.iter().cloned().collect();
    let mut pairs = Vec::new();
    
    let ic_values = match ic_results.get(&taxon_id) {
        Some(values) => values,
        None => {
            println!("Warning: No IC values found for taxon ID: {}\n", taxon_id);
            return pairs;
        }
    };
    
    println!("Generating all pairwise term combinations for taxon {}\n", taxon_id);
    
    for i in 0..terms.len() {
        for j in i..terms.len() {
            let term1 = terms[i];
            let term2 = terms[j];
            
            let ic_term1 = ic_values.get(&term1).copied().unwrap_or(f64::INFINITY);
            let ic_term2 = ic_values.get(&term2).copied().unwrap_or(f64::INFINITY);
            
            let mica = find_mica_for_pair(
                term1, term2,
                ontology_graph,
                go_id_to_node_index,
                node_index_to_go_id,
                ic_values
            ).unwrap_or((0, 0.0));
            
            pairs.push(TermPair::new(
                term1,
                term2,
                ic_term1,
                ic_term2,
                mica,
                method
            ));
            
        }
    }
    
    pairs
}

pub fn write_similarity_to_tsv(
    term_pairs: &[TermPair],
    go_terms: &FxHashSet<u32>,
    taxon_id: TaxonID,
    output_dir: &str,
) {
    let output_path = Path::new(output_dir);
    fs::create_dir_all(output_path).unwrap();

    let filename = format!("{}/similarity_taxon_{}.tsv", output_dir, taxon_id);
    println!("Writing similarity matrix to {}\n", filename);

    let mut file = File::create(&filename).unwrap();
   
    let terms: Vec<String> = go_terms.iter()
        .map(|&id| format!("GO:{:07}", id))
        .collect();

    let mut similarity_map: FxHashMap<(u32, u32), f64> = FxHashMap::default();
    for pair in term_pairs {
        similarity_map.insert((pair.term1, pair.term2), pair.similarity);
        similarity_map.insert((pair.term2, pair.term1), pair.similarity);
    }

    write!(file, "\t").unwrap();
    for term in &terms {
        write!(file, "{}\t", term).unwrap();
    }
    writeln!(file).unwrap();

    for (_, row_term_str) in terms.iter().enumerate() {
        let row_term = parse_single_go_term(&row_term_str).unwrap();
        
        write!(file, "{}\t", row_term_str).unwrap();
        
        for (_, col_term_str) in terms.iter().enumerate() {
            let col_term = parse_single_go_term(&col_term_str).unwrap();
            
            let similarity = similarity_map.get(&(row_term, col_term)).copied().unwrap_or(0.0);
            
            write!(file, "{:.6}\t", similarity).unwrap();
        }
        writeln!(file).unwrap();
    }

}