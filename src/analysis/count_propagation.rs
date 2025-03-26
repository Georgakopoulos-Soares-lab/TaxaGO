use std::collections::{HashMap, HashSet};
use std::io::Result;
use crate::parsers::{background_parser::*, study_parser::*};
use crate::parsers::obo_parser::*;
use daggy::{NodeIndex, Walker};
use::rayon::prelude::*;
use petgraph::algo::toposort;

#[derive(Debug, Default, Clone)]
pub struct GOAncestorCache {
    pub parent_map: HashMap<GOTermID, HashSet<GOTermID>>,
    pub propagation_order: Vec<GOTermID>, 
}

// NEEDS REFACTORING TO MINIMIZE MEMORY USAGE
impl GOAncestorCache {
    pub fn new(
        ontology_graph: &OntologyGraph,
        go_term_map: &OboMap,
        go_id_to_node_index: &HashMap<GOTermID, NodeIndex>,
        node_index_to_go_id: &HashMap<NodeIndex, GOTermID> 
    ) -> Result<Self> {
       
        let go_terms: Vec<&GOTermID> = go_term_map.keys().collect();

        let parent_map: HashMap<GOTermID, HashSet<GOTermID>> = go_terms
            .par_iter()
            .filter_map(|&go_term| {
                go_id_to_node_index.get(go_term).map(|&node_idx| {
                    let ancestors = get_unique_ancestors(
                        node_idx, ontology_graph, &node_index_to_go_id);
                    
                    (*go_term, ancestors)
                })
            })
            .collect();
        
        let topo_result = toposort(ontology_graph, None).unwrap();

        let mut propagation_order: Vec<GOTermID> = topo_result.iter()
                    .filter_map(|&node_idx| node_index_to_go_id.get(&node_idx).copied())
                    .collect();
                

        propagation_order.reverse();

        Ok(Self {
            parent_map,
            propagation_order
        })
    }
    
}

impl StudyPop {
    pub fn propagate_counts(
        &mut self,
        taxon_ids: &HashSet<TaxonID>,
        ancestor_cache: &GOAncestorCache
    ) -> () {

        let results: Vec<(TaxonID, GOTermCount, HashMap<GOTermID, HashSet<String>>)> = taxon_ids
            .par_iter()  
            .filter_map(|taxon_id| {

                let go_term_count = self.go_term_count.get(taxon_id)?;
                let go_term_protein_sets = self.go_term_to_protein_set.get(taxon_id)?;
                
                let mut local_counts = go_term_count.clone();
                let mut local_protein_sets = go_term_protein_sets.clone();
                
                for go_term in &ancestor_cache.propagation_order {

                    let current_proteins = match local_protein_sets.get(go_term) {
                        Some(proteins) => proteins.clone(), 
                        None => continue,
                    };
                    
                    let parent_terms = match ancestor_cache.parent_map.get(go_term) {
                        Some(parents) => parents.clone(),
                        None => continue,
                    };
                    
                    for parent in parent_terms {
                        let parent_proteins = local_protein_sets
                            .entry(parent)
                            .or_insert_with(HashSet::new);
                                                
                        for protein in &current_proteins {
                            parent_proteins.insert(protein.clone());
                        }
                        
                        let after_count = parent_proteins.len();
                        *local_counts.entry(parent).or_insert(0) = after_count;
                    }
                }
                
                Some((*taxon_id, local_counts, local_protein_sets))
            })
            .collect();
        
        for (taxon_id, counts, protein_sets) in results {
            if let Some(taxon_counts) = self.go_term_count.get_mut(&taxon_id) {
                *taxon_counts = counts;
            }
            
            if let Some(taxon_protein_sets) = self.go_term_to_protein_set.get_mut(&taxon_id) {
                *taxon_protein_sets = protein_sets;
            }
        }
    }
}

impl BackgroundPop {
    pub fn propagate_counts(
        &mut self,
        taxon_ids: &HashSet<TaxonID>,
        ancestor_cache: &GOAncestorCache
    ) -> () {
        let results: Vec<(TaxonID, GOTermCount, HashMap<GOTermID, HashSet<String>>)> = taxon_ids
            .par_iter()  
            .filter_map(|taxon_id| {
                let go_term_count = self.go_term_count.get(taxon_id)?;
                let go_term_protein_sets = self.go_term_to_protein_set.get(taxon_id)?;
                
                let mut local_counts = go_term_count.clone();
                let mut local_protein_sets = go_term_protein_sets.clone();
                
                for go_term in &ancestor_cache.propagation_order {
                    let current_proteins = match local_protein_sets.get(go_term) {
                        Some(proteins) => proteins.clone(),
                        None => continue,
                    };
                    
                    let parent_terms = match ancestor_cache.parent_map.get(go_term) {
                        Some(parents) => parents.clone(), 
                        None => continue,
                    };
                    
                    for parent in parent_terms {
                        let parent_proteins = local_protein_sets
                            .entry(parent)
                            .or_insert_with(HashSet::new);
                                                
                        for protein in &current_proteins {
                            parent_proteins.insert(protein.clone());
                        }
                        
                        let after_count = parent_proteins.len();
                        *local_counts.entry(parent).or_insert(0) = after_count;
                    }
                }
                
                Some((*taxon_id, local_counts, local_protein_sets))
            })
            .collect();
        
        for (taxon_id, counts, protein_sets) in results {
            if let Some(taxon_counts) = self.go_term_count.get_mut(&taxon_id) {
                *taxon_counts = counts;
            }
            
            if let Some(taxon_protein_sets) = self.go_term_to_protein_set.get_mut(&taxon_id) {
                *taxon_protein_sets = protein_sets;
            }
        }
    }
}

fn get_unique_ancestors(
    node_idx: NodeIndex,
    ontology_graph: &OntologyGraph,
    node_index_to_go_id: &HashMap<NodeIndex, GOTermID>
) -> HashSet<GOTermID> {

    let mut ancestors = HashSet::new();
    let mut to_visit = vec![node_idx];
    let mut visited = HashSet::new();
    
    while let Some(current_idx) = to_visit.pop() {
        if !visited.insert(current_idx) {
            continue;
        }
        
        let mut parents = ontology_graph.parents(current_idx);
        while let Some((edge_idx, parent_idx)) = parents.walk_next(ontology_graph) {
            match ontology_graph.edge_weight(edge_idx).unwrap() {
                Relationship::IsA | Relationship::PartOf => {
                    let go_term_id = node_index_to_go_id.get(&parent_idx).unwrap();
                    ancestors.insert(*go_term_id);
                    to_visit.push(parent_idx);
                },
                _ => continue,
            }
        }
    }
    
    ancestors
}