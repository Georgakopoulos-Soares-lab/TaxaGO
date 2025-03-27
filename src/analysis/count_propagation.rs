use rustc_hash::{FxHashMap, FxHashSet};
use std::io::Result;
use crate::parsers::{background_parser::*, study_parser::*};
use crate::parsers::obo_parser::*;
use daggy::{NodeIndex, Walker};
use std::sync::Arc;
use::rayon::prelude::*;
use petgraph::algo::toposort;

#[derive(Debug, Default, Clone)]
pub struct GOAncestorCache {
    pub parent_map: FxHashMap<GOTermID, FxHashSet<GOTermID>>,
    pub propagation_order: Vec<GOTermID>, 
}

impl GOAncestorCache {
    pub fn new(
        ontology_graph: &OntologyGraph,
        go_term_map: &OboMap,
        go_id_to_node_index: &FxHashMap<GOTermID, NodeIndex>,
        node_index_to_go_id: &FxHashMap<NodeIndex, GOTermID> 
    ) -> Result<Self> {
        
        let go_terms: Vec<&GOTermID> = go_term_map.keys().collect();

        let parent_map: FxHashMap<GOTermID, FxHashSet<GOTermID>> = go_terms
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

        let mut propagation_order: Vec<GOTermID> = topo_result
            .iter()
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
        taxon_ids: &FxHashSet<TaxonID>,
        ancestor_cache: &GOAncestorCache
    ) {
        taxon_ids.iter().for_each(|taxon_id| {
            let (Some(go_term_count), Some(go_term_protein_sets)) = (
                self.go_term_count.get_mut(taxon_id),
                self.go_term_to_protein_set.get_mut(taxon_id)
            ) else { return };
            
            ancestor_cache.propagation_order.iter().for_each(|&go_term| {
                let Some(parent_terms) = ancestor_cache.parent_map.get(&go_term) else { return };
                if parent_terms.is_empty() {
                    return;
                }
                
                let Some(proteins) = go_term_protein_sets.get(&go_term) else { return };
                if proteins.is_empty() {
                    return;
                }
                
                let term_proteins: Vec<Protein> = proteins.iter().cloned().collect();
                
                parent_terms.iter().for_each(|&parent| {
                    let parent_proteins = go_term_protein_sets
                        .entry(parent)
                        .or_insert_with(FxHashSet::default);
                                        
                    term_proteins.iter().for_each(|protein| {
                        parent_proteins.insert(Arc::clone(protein));
                    });
                    
                    *go_term_count.entry(parent).or_insert(0) = parent_proteins.len();
                });
            });
        });
    }
}

impl BackgroundPop {
    pub fn propagate_counts(
        &mut self,
        taxon_ids: &FxHashSet<TaxonID>,
        ancestor_cache: &GOAncestorCache
    ) {
        taxon_ids.iter().for_each(|taxon_id| {
            let (Some(go_term_count), Some(go_term_protein_sets)) = (
                self.go_term_count.get_mut(taxon_id),
                self.go_term_to_protein_set.get_mut(taxon_id)
            ) else { return };
            
            ancestor_cache.propagation_order.iter().for_each(|&go_term| {
                let Some(parent_terms) = ancestor_cache.parent_map.get(&go_term) else { return };
                if parent_terms.is_empty() {
                    return;
                }
                
                let Some(proteins) = go_term_protein_sets.get(&go_term) else { return };
                if proteins.is_empty() {
                    return;
                }
                
                let term_proteins: Vec<Protein> = proteins.iter().cloned().collect();
                
                parent_terms.iter().for_each(|&parent| {
                    let parent_proteins = go_term_protein_sets
                        .entry(parent)
                        .or_insert_with(FxHashSet::default);
                                        
                    term_proteins.iter().for_each(|protein| {
                        parent_proteins.insert(Arc::clone(protein));
                    });
                    
                    *go_term_count.entry(parent).or_insert(0) = parent_proteins.len();
                });
            });
        });
    }
}

fn get_unique_ancestors(
    node_idx: NodeIndex,
    ontology_graph: &OntologyGraph,
    node_index_to_go_id: &FxHashMap<NodeIndex, GOTermID>
) -> FxHashSet<GOTermID> {

    let mut ancestors = FxHashSet::default();
    let mut to_visit = vec![node_idx];
    let mut visited = FxHashSet::default();
    
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