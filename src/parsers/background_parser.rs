use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use rayon::prelude::*;
use daggy::{NodeIndex, Walker};
use crate::parsers::obo_parser::*;

pub type TaxonID = u32;
pub type ProteinCount = HashMap<TaxonID, usize>;
pub type ProteinToGO = HashMap<String, HashSet<u32>>;
pub type GOTermCount = HashMap<u32, usize>;

#[derive(Debug, Default, Clone)]
pub struct BackgroundPop {
    pub taxon_protein_count: ProteinCount,
    pub protein_to_go: HashMap<TaxonID, ProteinToGO>,
    pub go_term_count: HashMap<TaxonID, GOTermCount>,
}

impl BackgroundPop {
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

        let updated_counts: HashMap<_, _> = self.go_term_count
            .par_iter()
            .map(|(&taxon_id, go_term_counts)| {
                let mut updated_counts = HashMap::with_capacity(go_term_counts.len() * 2);
                
                for (&go_id, &direct_count) in go_term_counts.iter() {
                    if let Some(&node_idx) = go_id_to_node_index.get(&go_id) {
                        *updated_counts.entry(go_id).or_insert(0) += direct_count;
                        
                        if let Some(ancestors) = ancestor_cache.get(&node_idx) {
                            for &ancestor_idx in ancestors {
                                if let Some(&ancestor_id) = node_index_to_go_id.get(&ancestor_idx) {
                                    *updated_counts.entry(ancestor_id).or_insert(0) += direct_count;
                                }
                            }
                        }
                    }
                }
                
                (taxon_id, updated_counts)
            })
            .collect();
        
        self.go_term_count = updated_counts;
    }
}

fn get_unique_ancestors(
    node_idx: NodeIndex,
    graph: &OntologyGraph,
) -> HashSet<NodeIndex> {
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

fn process_single_taxon(taxon_background_path: impl AsRef<Path>) -> Result<Option<(usize, ProteinToGO, GOTermCount)>, String> {
    let path = taxon_background_path.as_ref();
    
    if !path.is_file() {
        return Ok(None);
    }

    let file = File::open(path).map_err(|e| format!("Failed to open file: {}", e))?;
    let reader = BufReader::with_capacity(128 * 1024, file);

    let mut protein_to_go_map = HashMap::new();
    let mut go_term_counts = HashMap::new();

    for line_result in reader.lines() {
        let line = line_result.map_err(|e| format!("Error reading line: {}", e))?;
        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() < 2 {
            continue;
        }

        let protein_id = parts[0].to_string();
        if let Some(go_str) = parts[1].strip_prefix("GO:") {
            if let Ok(go_id) = go_str.parse::<u32>() {
                protein_to_go_map
                    .entry(protein_id)  
                    .or_insert_with(HashSet::new)
                    .insert(go_id);

                *go_term_counts.entry(go_id).or_insert(0) += 1;
            }
        }
    }

    Ok(Some((protein_to_go_map.len(), protein_to_go_map, go_term_counts)))
}

pub fn read_background_pop(taxon_ids: &HashSet<TaxonID>, dir: &str) -> BackgroundPop {
    let results: Vec<(TaxonID, Option<(usize, ProteinToGO, GOTermCount)>)> = taxon_ids
        .par_iter()
        .map(|&taxon_id| {
            let taxon_background_file_path = format!("{}/{}_background.txt", dir, taxon_id);
            let result = process_single_taxon(&taxon_background_file_path);
            (taxon_id, result.unwrap_or(None))
        })
        .collect();

    let mut background_pop = BackgroundPop::default();
    
    for (taxon_id, result) in results {
        if let Some((protein_count, protein_to_go_map, go_term_counts)) = result {
            background_pop.taxon_protein_count.insert(taxon_id, protein_count);
            background_pop.protein_to_go.insert(taxon_id, protein_to_go_map);
            background_pop.go_term_count.insert(taxon_id, go_term_counts);
        } else {
            eprintln!("No valid data found for taxon {}", taxon_id);
        }
    }

    background_pop
}
