use rustc_hash::{FxHashMap, FxHashSet};
use crate::parsers::{
    background_parser::*,
    study_parser::*,
    obo_parser::*
};
use crate::analysis::enrichment_analysis::*;
use rayon::prelude::*;
use daggy::{NodeIndex, Walker};

type Weights = FxHashMap<GOTermID, FxHashMap<Protein, f64>>;

fn sig_ratio_log_p_values(
    p_value1: f64,
    p_value2: f64
) -> f64 {

    if p_value1 == 0.0 || p_value2 == 0.0{
        if p_value1 == 0.0 && p_value2 == 0.0 { return 1.0; }
        if p_value1 == 0.0 && p_value2 != 0.0 { return 1.0; }
        if p_value2 == 0.0 && p_value1 != 0.0{ return 0.0; }
    }

    let log1 = p_value1.ln();
    let log2 = p_value2.ln();
    
    if log2 == 0.0 { 
        if log1 == 0.0 { return 1.0; } 
        return 1.0;
    }

    if log1 == 0.0 { 
        return 0.0; 
    }
    
    let factor = log1 / log2;
    factor.max(0.0).min(1.0)
}

impl EnrichmentAnalysis {
    pub fn weight(
        &self,
        taxon_ids: &FxHashSet<TaxonID>,
        study_pop: &StudyPop,
        background_pop: &BackgroundPop,
        level_to_go_term: &LevelToTerms,
        ontology_graph: &OntologyGraph,
        go_id_to_node_index: &FxHashMap<GOTermID, NodeIndex>,
        node_index_to_go_id: &FxHashMap<NodeIndex, GOTermID>,
    ) -> FxHashMap<TaxonID, FxHashMap<GOTermID, GOTermResults>> {

        let max_level = match level_to_go_term.keys().max() {
            Some(&level) => level,
            None => return FxHashMap::default(),
        };

        taxon_ids
            .par_iter()
            .map(|&taxon_id| {

                let mut node_weights: Weights = FxHashMap::default();
                
                let taxon_study_go_term_to_proteins = match study_pop.go_term_to_protein_set.get(&taxon_id) {
                    Some(data) => data,
                    None => return (taxon_id, FxHashMap::default()),
                };

                for (go_term_id, proteins_in_term) in taxon_study_go_term_to_proteins.iter() {
                    let mut term_weights: FxHashMap<Protein, f64> = FxHashMap::default();
                    for protein in proteins_in_term.iter() {
                        term_weights.insert(protein.clone(), 1.0);
                    }
                    node_weights.insert(*go_term_id, term_weights);
                }

                let mut go_term_results: FxHashMap<GOTermID, GOTermResults> = FxHashMap::default();

                for level in (1..=max_level).rev() {
                    if let Some(go_terms_at_level) = level_to_go_term.get(&level) {
                        for &term_id in go_terms_at_level {
                            if taxon_study_go_term_to_proteins.contains_key(&term_id) && !go_term_results.contains_key(&term_id) {
                                let initial_score_result = self.calculate_weighted_score(
                                    term_id,
                                    taxon_id,
                                    study_pop,
                                    background_pop,
                                    node_weights.get(&term_id),
                                );
                                go_term_results.insert(term_id, initial_score_result);
                            }
                        }
                    }
                }
                
                for level in (1..=max_level).rev() {
                    if let Some(go_terms_at_level) = level_to_go_term.get(&level) {
                        for &current_go_term_id in go_terms_at_level {
                            if !taxon_study_go_term_to_proteins.contains_key(&current_go_term_id) {
                                continue;
                            }
                            let parent_idx= go_id_to_node_index.get(&current_go_term_id).unwrap();
                            let children_ids = self.get_direct_children(
                                ontology_graph,
                                *parent_idx,
                                node_index_to_go_id
                            );
                            
                            self.compute_term_sig(
                                current_go_term_id,
                                &children_ids,
                                &mut node_weights,
                                &mut go_term_results,
                                taxon_id,
                                study_pop,
                                background_pop,
                                ontology_graph,
                                go_id_to_node_index,
                                node_index_to_go_id
                            );
                        }
                    }
                }
                (taxon_id, go_term_results)
            })
            .collect()
    }

    #[allow(clippy::too_many_arguments)]
    fn compute_term_sig(
        &self,
        term_id: GOTermID,
        current_children_ids: &FxHashSet<GOTermID>,
        node_weights: &mut Weights,
        go_term_results: &mut FxHashMap<GOTermID, GOTermResults>,
        taxon_id: TaxonID,
        study_pop: &StudyPop,
        background_pop: &BackgroundPop,
        ontology_graph: &OntologyGraph,
        go_id_to_node_index: &FxHashMap<GOTermID, NodeIndex>,
        node_index_to_go_id: &FxHashMap<NodeIndex, GOTermID>,
    ) {
        let results = self.calculate_weighted_score(
            term_id,
            taxon_id,
            study_pop,
            background_pop,
            node_weights.get(&term_id),
        );

        let current_p_value = results.p_value;
        go_term_results.insert(term_id, results);

        if current_children_ids.is_empty() {
            return;
        }

        let mut more_significant_children_ids: FxHashSet<GOTermID> = FxHashSet::default();
        for &child_id in current_children_ids {
            if let Some(child_result) = go_term_results.get(&child_id) {
                if child_result.p_value < current_p_value {
                    more_significant_children_ids.insert(child_id);
                }
            }
        }

        if more_significant_children_ids.is_empty() {
            for &child_id in current_children_ids {
                if let Some(child_result) = go_term_results.get(&child_id) {
                    let p_value_child = child_result.p_value;
                    if p_value_child > current_p_value {
                        let factor = sig_ratio_log_p_values(p_value_child, current_p_value);
                        
                        Self::apply_weight_to_term_proteins(
                            child_id,
                            factor,
                            taxon_id,
                            study_pop,
                            node_weights,
                        );

                        let updated_child_results = self.calculate_weighted_score(
                            child_id, taxon_id, study_pop, background_pop, node_weights.get(&child_id),
                        );
                        go_term_results.insert(child_id, updated_child_results);
                    }
                }
            }
            return;
        }
        else {
            for &better_child_id in &more_significant_children_ids {
                if let Some(better_child_result) = go_term_results.get(&better_child_id) {
                    let p_value_better_child = better_child_result.p_value;
                    let factor = sig_ratio_log_p_values(current_p_value, p_value_better_child);
                    let node_idx= go_id_to_node_index.get(&better_child_id).unwrap();
                    let ancestors_of_u_inclusive = self.get_ancestors_including_term(
                        *node_idx,
                        ontology_graph, 
                        node_index_to_go_id, 
                    );
                    
                    if let Some(proteins_in_better_child_set) = study_pop.go_term_to_protein_set.get(&taxon_id).and_then(|tg| tg.get(&better_child_id)) {
                        for ancestor_id in ancestors_of_u_inclusive {
                            if let Some(proteins_in_ancestor_set) = study_pop.go_term_to_protein_set.get(&taxon_id).and_then(|tg| tg.get(&ancestor_id)) {
                                if let Some(ancestor_specific_weights) = node_weights.get_mut(&ancestor_id) {
                                    for protein_in_ancestor in proteins_in_ancestor_set.iter() {
                                        if proteins_in_better_child_set.contains(protein_in_ancestor) {
                                            if let Some(weight_val) = ancestor_specific_weights.get_mut(protein_in_ancestor) {
                                                *weight_val *= factor;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            let remaining_children_ids: FxHashSet<GOTermID> = current_children_ids
                .difference(&more_significant_children_ids)
                .copied()
                .collect();
            
            self.compute_term_sig(
                term_id,
                &remaining_children_ids,
                node_weights,
                go_term_results,
                taxon_id,
                study_pop,
                background_pop,
                ontology_graph,
                go_id_to_node_index,
                node_index_to_go_id,
            );
        }
    }

    fn apply_weight_to_term_proteins(
        term_id: GOTermID,
        factor: f64,
        taxon_id: TaxonID,
        study_pop: &StudyPop,
        node_weights: &mut Weights,
    ) {

        if let Some(proteins_in_term_map) = study_pop.go_term_to_protein_set.get(&taxon_id).and_then(|tg| tg.get(&term_id)) {
            if let Some(term_weights) = node_weights.get_mut(&term_id) {
                for protein in proteins_in_term_map.iter() {
                    if let Some(weight) = term_weights.get_mut(protein) {
                        *weight *= factor;
                    }
                }
            }
        }
    }

    fn calculate_weighted_score(
        &self,
        term_id: GOTermID,
        taxon_id: TaxonID,
        study_pop: &StudyPop,
        background_pop: &BackgroundPop,
        term_specific_protein_weights: Option<&FxHashMap<Protein, f64>>,
    ) -> GOTermResults {

        let study_proteins_in_term = study_pop
            .go_term_to_protein_set
            .get(&taxon_id)
            .and_then(|tg| tg.get(&term_id))
            .cloned()
            .unwrap();

        let all_study_proteins = study_pop
            .taxon_map
            .get(&taxon_id)
            .cloned()
            .unwrap();

        let mut weighted_study_proteins_in_term_float = 0.0;

        if let Some(weights) = term_specific_protein_weights {
            for protein_in_term in &study_proteins_in_term {
                let weight = weights.get(protein_in_term).copied().unwrap_or(1.0);
                weighted_study_proteins_in_term_float += weight;
            }
        } else {
            weighted_study_proteins_in_term_float += study_proteins_in_term.len() as f64;
        }
        let weighted_study_proteins_in_term = weighted_study_proteins_in_term_float.round() as usize;

        let background_total_in_term_count = background_pop
            .go_term_count
            .get(&taxon_id)
            .and_then(|tgc| tgc.get(&term_id))
            .copied()
            .unwrap();

        let total_study_proteins_count = all_study_proteins.len();
        
        let total_background_proteins_count = background_pop
            .taxon_protein_count
            .get(&taxon_id)
            .copied()
            .unwrap();

        let contingency_table = create_contingency_table(
            weighted_study_proteins_in_term,
            background_total_in_term_count,
            total_study_proteins_count,
            total_background_proteins_count,
        );

        GOTermResults {
            log_odds_ratio: calculate_log_odds_ratio(&contingency_table),
            p_value: calculate_p_value(&contingency_table, self.test_type),
            contingency_table,
            variance: calculate_variance(&contingency_table),
        }
    }

    fn get_ancestors_including_term(
        &self,
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
        ancestors.insert(*node_index_to_go_id.get(&node_idx).unwrap()); 
        ancestors
    }

    fn get_direct_children(
        &self,
        ontology_graph: &OntologyGraph,
        parent_idx: NodeIndex,
        node_index_to_go_id: &FxHashMap<NodeIndex, GOTermID>
    ) -> FxHashSet<GOTermID> {

        let mut children_set: FxHashSet<GOTermID>  = FxHashSet::default();
        let mut children = ontology_graph.children(parent_idx);
       
        while let Some((edge_idx, child_idx)) = children.walk_next(ontology_graph) {
            match ontology_graph.edge_weight(edge_idx).unwrap() {
                Relationship::IsA | Relationship::PartOf => {
                    let go_term_id = node_index_to_go_id.get(&child_idx).unwrap();
                    children_set.insert(*go_term_id);
                },
                _ => continue,
            }
        }
        
        children_set
    }
}