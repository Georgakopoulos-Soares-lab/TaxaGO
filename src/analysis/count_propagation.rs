use std::collections::{HashMap, HashSet};
use crate::parsers::study_parser::StudyPop;
use crate::parsers::obo_parser::*;
use daggy::{NodeIndex, Walker};
use::rayon::prelude::*;

impl StudyPop {
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