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
use std::collections::VecDeque;
use daggy::Walker;

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

impl TermPair {
    pub fn new_for_wang(term1: u32, term2: u32, wang_similarity_score: f64) -> Self {
        Self {
            term1,
            term2,
            ic_term1: 0.0,
            ic_term2: 0.0,
            mica: (0, 0.0),
            similarity: wang_similarity_score,
        }
    }

    pub fn new_for_ic(
        term1: u32,
        term2: u32,
        ic_term1: f64,
        ic_term2: f64,
        mica: (u32, f64),
        method: &str, 
    ) -> Self {
        let mut pair = Self {
            term1,
            term2,
            ic_term1,
            ic_term2,
            mica,
            similarity: 0.0,
        };
        pair.similarity = match method.to_lowercase().as_str() {
            "resnik" => pair.mica.1, 
            "lin" => {
                if pair.ic_term1 == 0.0 || pair.ic_term2 == 0.0 || (pair.ic_term1 + pair.ic_term2) == 0.0 { 0.0 }
                else { (2.0 * pair.mica.1) / (pair.ic_term1 + pair.ic_term2) }
            }
            "jiang-conrath" => {
                let distance = pair.ic_term1 + pair.ic_term2 - 2.0 * pair.mica.1;
                1.0 / (1.0 + distance.max(0.0))
            }
            _ => {
                eprintln!("Warning: Unsupported IC-based method '{}' in TermPair::new_for_ic, defaulting to 0.0", method);
                0.0
            }
        };
        pair
    }
}


pub const IS_A_WEIGHT: f64 = 0.8;
pub const PART_OF_WEIGHT: f64 = 0.6;

pub fn get_edge_weight(relationship: &Relationship) -> f64 {
    match relationship {
        Relationship::IsA => IS_A_WEIGHT,
        Relationship::PartOf => PART_OF_WEIGHT,
        Relationship::Regulates => 0.0,
        Relationship::PositivelyRegulates => 0.0,
        Relationship::NegativelyRegulates => 0.0,
        Relationship::OccursIn => 0.0,
    }
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
    global_rev_topo_order: &[GOTermID], 
    method: &str,
) -> Vec<TermPair> {
    let terms_vec: Vec<u32> = go_terms.iter().cloned().collect();
    let mut pairs = Vec::new();

    let method_lower = method.to_lowercase();

    if method_lower == "wang" {
        println!("Calculating Wang's similarity for term pairs...");
        for i in 0..terms_vec.len() {
            for j in i..terms_vec.len() {
                let term1 = terms_vec[i];
                let term2 = terms_vec[j];

                match wang_similarity(
                    term1,
                    term2,
                    ontology_graph,
                    go_id_to_node_index,
                    node_index_to_go_id,
                    global_rev_topo_order,
                ) {
                    Ok(sim_score) => {
                        pairs.push(TermPair::new_for_wang(term1, term2, sim_score));
                    }
                    Err(e) => {
                        eprintln!("Error calculating Wang's similarity for GO:{:07} and GO:{:07}: {}", term1, term2, e);
                        pairs.push(TermPair::new_for_wang(term1, term2, 0.0));
                    }
                }
            }
        }
    } else {
        let ic_values_for_taxon = match ic_results.get(&taxon_id) {
            Some(values) => values,
            None => {
                eprintln!("Warning: No IC values found for taxon ID: {} for method {}. Returning empty pairs.", taxon_id, method);
                return pairs;
            }
        };

        println!("Calculating {} similarity for term pairs using IC...", method);
        for i in 0..terms_vec.len() {
            for j in i..terms_vec.len() {
                let term1 = terms_vec[i];
                let term2 = terms_vec[j];

                let ic_term1 = ic_values_for_taxon.get(&term1).copied().unwrap_or(0.0);
                let ic_term2 = ic_values_for_taxon.get(&term2).copied().unwrap_or(0.0);

                let mica = find_mica_for_pair(
                    term1, term2,
                    ontology_graph,
                    go_id_to_node_index,
                    node_index_to_go_id,
                    ic_values_for_taxon,
                ).unwrap_or((0, 0.0));

                pairs.push(TermPair::new_for_ic(
                    term1, term2, ic_term1, ic_term2, mica, &method_lower,
                ));
            }
        }
    }
    pairs
}

pub fn calculate_s_values(
    term_id: GOTermID,
    ontology_graph: &OntologyGraph,
    go_id_to_node_index: &FxHashMap<GOTermID, NodeIndex>,
    node_index_to_go_id: &FxHashMap<NodeIndex, GOTermID>,
    topo_order: &[GOTermID],
) -> Result<FxHashMap<GOTermID, f64>, String> {
    let term_a_node_idx = match go_id_to_node_index.get(&term_id) {
        Some(&idx) => idx,
        None => return Err(format!("Term GO:{:07} not found in ontology", term_id)),
    };

    let mut ancestors: FxHashSet<GOTermID> = FxHashSet::default();
    let mut ancestor_traversal_queue = VecDeque::new();
    ancestor_traversal_queue.push_back(term_a_node_idx);
    ancestors.insert(term_id);
    let mut head = 0;
    while head < ancestor_traversal_queue.len() {
        let current_node_idx = ancestor_traversal_queue[head];
        head += 1;
        let mut parents = ontology_graph.parents(current_node_idx);
        while let Some((edge_idx, parent_node_idx)) = parents.walk_next(ontology_graph) {
            if let Some(relationship) = ontology_graph.edge_weight(edge_idx) {
                match relationship {
                    Relationship::IsA | Relationship::PartOf => {
                        if let Some(parent_go_id) = node_index_to_go_id.get(&parent_node_idx) {
                            if ancestors.insert(*parent_go_id) {
                                ancestor_traversal_queue.push_back(parent_node_idx);
                            }
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    let mut s_values: FxHashMap<GOTermID, f64> = FxHashMap::default();
    for &ancestor_id_in_dag_a in &ancestors {
        s_values.insert(ancestor_id_in_dag_a, 0.0);
    }
    s_values.insert(term_id, 1.0);

    for &current_go_id in topo_order.iter() {
        if !ancestors.contains(&current_go_id) {
            continue;
        }

        if current_go_id == term_id {
            continue;
        }

        let current_node_idx = match go_id_to_node_index.get(&current_go_id) {
            Some(idx) => *idx,
            None => continue, 
        };
        
        let mut max_s_contrib_for_current_go_id: f64 = 0.0;

        let mut children_walker = ontology_graph.children(current_node_idx);
        while let Some((edge_to_child_idx, child_node_idx)) = children_walker.walk_next(ontology_graph) {
            if let Some(child_go_id_in_ontology) = node_index_to_go_id.get(&child_node_idx) {
                if ancestors.contains(child_go_id_in_ontology) {
                    if let Some(s_value_of_child) = s_values.get(child_go_id_in_ontology) {
                        if let Some(relationship) = ontology_graph.edge_weight(edge_to_child_idx) {
                            match relationship {
                                Relationship::IsA | Relationship::PartOf => {
                                    let weight: f64 = get_edge_weight(relationship);
                                    max_s_contrib_for_current_go_id = max_s_contrib_for_current_go_id.max(*s_value_of_child * weight);
                                }
                                _ => {} 
                            }
                        }
                    }
                }
            }
        }
        s_values.insert(current_go_id, max_s_contrib_for_current_go_id);
    }

    Ok(s_values)
}

pub fn calculate_semantic_value(
    s_values: &FxHashMap<GOTermID, f64>,
) -> f64 {
    s_values.values().sum()
}


pub fn wang_similarity(
    term1_id: GOTermID,
    term2_id: GOTermID,
    ontology_graph: &OntologyGraph,
    go_id_to_node_index: &FxHashMap<GOTermID, NodeIndex>,
    node_index_to_go_id: &FxHashMap<NodeIndex, GOTermID>,
    global_rev_topo_order: &[GOTermID], 
) -> Result<f64, String> {
    if term1_id == term2_id {
        return Ok(1.0); 
    }

    let s_values_term1 = calculate_s_values(
        term1_id,
        ontology_graph,
        go_id_to_node_index,
        node_index_to_go_id,
        global_rev_topo_order,
    )?;

    let s_values_term2 = calculate_s_values(
        term2_id,
        ontology_graph,
        go_id_to_node_index,
        node_index_to_go_id,
        global_rev_topo_order,
    )?;

    let sv_term1 = calculate_semantic_value(&s_values_term1);

    let sv_term2 = calculate_semantic_value(&s_values_term2);

    if sv_term1 == 0.0 || sv_term2 == 0.0 {
        return Ok(0.0);
    }

    let mut sum_common_s_values = 0.0;

    for (ancestor_t_id, s_value_for_t_in_term1) in &s_values_term1 {
        if let Some(s_value_for_t_in_term2) = s_values_term2.get(ancestor_t_id) {
            sum_common_s_values += s_value_for_t_in_term1 + s_value_for_t_in_term2;
        }
    }

    let denominator = sv_term1 + sv_term2;
    if denominator == 0.0 {
        return Ok(0.0);
    }

    let similarity = sum_common_s_values / denominator;

    Ok(similarity)
}

pub fn write_similarity_to_tsv(
    term_pairs: &[TermPair],
    go_terms: &FxHashSet<GOTermID>,
    taxon_id: TaxonID,
    method: &str,
    output_dir: &str,
) -> Result<(), String> {
    let output_path = Path::new(output_dir);
    fs::create_dir_all(output_path)
        .map_err(|e| format!("Failed to create output directory {}: {}", output_dir, e))?;

    let method_filename_part = method
        .to_lowercase()
        .replace(|c: char| !c.is_alphanumeric() && c != '_', "_")
        .replace('-', "_");

    let filename = format!(
        "{}/similarity_{}_taxon_{}.tsv",
        output_dir, method_filename_part, taxon_id
    );
    println!("Writing {} similarity matrix to {}", method, filename);

    let mut file = File::create(&filename)
        .map_err(|e| format!("Failed to create output file {}: {}", filename, e))?;

    let mut sorted_go_ids: Vec<GOTermID> = go_terms.iter().cloned().collect();
    sorted_go_ids.sort_unstable();


    let mut similarity_map: FxHashMap<(GOTermID, GOTermID), f64> = FxHashMap::default();
    for pair in term_pairs {
        similarity_map.insert((pair.term1, pair.term2), pair.similarity);
        if pair.term1 != pair.term2 {
            similarity_map.insert((pair.term2, pair.term1), pair.similarity);
        }
    }

    write!(file, "\t")
        .map_err(|e| format!("Failed to write to file {}: {}", filename, e))?;

    for (i, go_id) in sorted_go_ids.iter().enumerate() {
        write!(file, "GO:{:07}", go_id)
            .map_err(|e| format!("Failed to write to file {}: {}", filename, e))?;
        if i < sorted_go_ids.len() - 1 {
            write!(file, "\t")
                .map_err(|e| format!("Failed to write to file {}: {}", filename, e))?;
        }
    }
    writeln!(file)
        .map_err(|e| format!("Failed to write to file {}: {}", filename, e))?;

    for row_go_id in &sorted_go_ids {
        write!(file, "GO:{:07}\t", row_go_id)
            .map_err(|e| format!("Failed to write to file {}: {}", filename, e))?;

        for (j, col_go_id) in sorted_go_ids.iter().enumerate() {
            let similarity = similarity_map
                .get(&(*row_go_id, *col_go_id))
                .copied()
                .unwrap_or(0.0); 

            write!(file, "{:.6}", similarity)
                .map_err(|e| format!("Failed to write to file {}: {}", filename, e))?;
            if j < sorted_go_ids.len() - 1 {
                write!(file, "\t")
                    .map_err(|e| format!("Failed to write to file {}: {}", filename, e))?;
            }
        }
        writeln!(file)
            .map_err(|e| format!("Failed to write to file {}: {}", filename, e))?;
    }

    println!("Successfully wrote similarity matrix to {}", filename);
    Ok(())
}