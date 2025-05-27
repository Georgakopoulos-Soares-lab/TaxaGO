use rustc_hash::{FxHashMap, FxHashSet};
use crate::parsers::{
    background_parser::*,
    study_parser::*,
    obo_parser::*
};
use crate::analysis::enrichment_analysis::*;
use rayon::prelude::*;
use daggy::NodeIndex;

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
    pub fn weight_analysis(
        &self,
        taxon_ids: &FxHashSet<TaxonID>,
        study_pop: &StudyPop,
        background_pop: &BackgroundPop,
        level_to_go_term: &LevelToTerms,
        obo_map: &OboMap,
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
                        term_weights.insert(*protein, 1.0);
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
                                    *study_pop,
                                    background_pop,
                                    node_weights.get(&term_id),
                                );
                                go_term_results.insert(term_id, initial_score_result);
                            }
                        }
                    }
                }
                
                // Second pass: Apply the recursive weight adjustment logic.
                for level in (1..=max_level).rev() {
                    if let Some(go_terms_at_level) = level_to_go_term.get(&level) {
                        for &current_go_term_id in go_terms_at_level {
                             // Ensure we only process terms relevant to the study for this taxon
                            if !taxon_study_go_term_to_proteins.contains_key(&current_go_term_id) {
                                continue;
                            }
                            let children_ids = get_direct_children_from_graph( // Assumes this utility exists
                                ontology_graph,
                                current_go_term_id,
                                go_id_to_node_index,
                                node_index_to_go_id,
                                obo_map
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
                                node_index_to_go_id,
                                obo_map,
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
        term_id_u: GOTermID,
        current_children_ids: &FxHashSet<GOTermID>, // Children of u to consider in this iteration
        node_weights: &mut Weights,
        go_term_results: &mut FxHashMap<GOTermID, GOTermResults>,
        taxon_id: TaxonID,
        study_pop: &StudyPop,
        background_pop: &BackgroundPop,
        ontology_graph: &OntologyGraph,
        go_id_to_node_index: &FxHashMap<GOTermID, NodeIndex>,
        node_index_to_go_id: &FxHashMap<NodeIndex, GOTermID>,
        obo_map: &OboMap,
    ) {
        // 1. (Re)calculate current p-value for term_id_u based on its current weights.
        let results_u = self.calculate_weighted_score(
            term_id_u,
            taxon_id,
            *study_pop,
            background_pop,
            node_weights.get(&term_id_u),
        );
        let current_p_value_u = results_u.p_value;
        go_term_results.insert(term_id_u, results_u);

        // Base Case for recursion: if no children to consider in this iteration.
        if current_children_ids.is_empty() {
            return;
        }

        // 3. Identify children (from current_children_ids) that are more significant than u.
        let mut more_significant_children_ids: FxHashSet<GOTermID> = FxHashSet::default();
        for &child_id in current_children_ids {
            if let Some(child_result) = go_term_results.get(&child_id) {
                if child_result.p_value < current_p_value_u {
                    more_significant_children_ids.insert(child_id);
                }
            }
        }

        // 4. Scenario A (Prose: "If u is more significant, genes contained in the children are down-weighted")
        // This means u is better than or equal to all children considered in this iteration.
        if more_significant_children_ids.is_empty() {
            for &child_id in current_children_ids { // All these children are less significant or equal to u
                if let Some(child_result) = go_term_results.get(&child_id) {
                    let p_value_child = child_result.p_value;
                    // Only down-weight child if it's strictly less significant (worse) than u
                    if p_value_child > current_p_value_u {
                        let factor = sig_ratio_log_p_values(p_value_child, current_p_value_u); // log(worse)/log(better) -> <1
                        
                        Self::apply_weight_to_term_genes(
                            child_id,
                            factor,
                            taxon_id,
                            study_pop,
                            node_weights,
                        );

                        let updated_child_results = self.calculate_weighted_score(
                            child_id, taxon_id, *study_pop, background_pop, node_weights.get(&child_id),
                        );
                        go_term_results.insert(child_id, updated_child_results);
                    }
                }
            }
            return; // All children processed relative to this u's score.
        }
        // 5. Scenario B (Prose: "If at least one child is more significant than u...")
        else {
            for &better_child_id in &more_significant_children_ids {
                if let Some(better_child_result) = go_term_results.get(&better_child_id) {
                    let p_value_better_child = better_child_result.p_value;
                    // Factor to down-weight u because better_child explains genes better.
                    // current_p_value_u is "worse" than p_value_better_child.
                    let factor = sig_ratio_log_p_values(current_p_value_u, p_value_better_child); // log(worse)/log(better) -> <1

                    let ancestors_of_u_inclusive = get_ancestors_including_self_from_graph(
                        ontology_graph, term_id_u, go_id_to_node_index, node_index_to_go_id, obo_map
                    );
                    
                    if let Some(proteins_in_better_child_set) = study_pop.go_term_to_protein_set.get(&taxon_id).and_then(|tg| tg.get(&better_child_id)) {
                        for ancestor_id in ancestors_of_u_inclusive {
                            if let Some(proteins_in_ancestor_set) = study_pop.go_term_to_protein_set.get(&taxon_id).and_then(|tg| tg.get(&ancestor_id)) {
                                if let Some(ancestor_specific_weights) = node_weights.get_mut(&ancestor_id) {
                                    for protein_in_ancestor in proteins_in_ancestor_set.iter() {
                                        if proteins_in_better_child_set.contains(protein_in_ancestor) { // Common gene
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

            // Recursive call for u with children that were NOT more significant than the original current_p_value_u
            let remaining_children_ids: FxHashSet<GOTermID> = current_children_ids
                .difference(&more_significant_children_ids)
                .copied()
                .collect();
            
            self.compute_term_sig(
                term_id_u,
                &remaining_children_ids, // Process u against the less significant children
                node_weights,
                go_term_results,
                taxon_id,
                study_pop,
                background_pop,
                ontology_graph,
                go_id_to_node_index,
                node_index_to_go_id,
                obo_map,
            );
        }
    }

    fn apply_weight_to_term_genes(
        term_id: GOTermID,
        factor: f64,
        taxon_id: TaxonID,
        study_pop: &StudyPop,
        node_weights: &mut Weights,
    ) {
        if factor >= 1.0 { return; } // Only apply actual down-weighting factors

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
        study_pop: StudyPop,
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

        let mut weighted_study_sig_in_term_float = 0.0;

        if let Some(weights) = term_specific_protein_weights {
            for protein_in_term in &study_proteins_in_term {
                let weight = weights.get(protein_in_term).copied().unwrap_or(1.0);
                weighted_study_sig_in_term_float += weight;
            }
        } else {
            weighted_study_sig_in_term_float += study_proteins_in_term.len() as f64;
        }
        let study_sig_in_term_count = weighted_study_sig_in_term_float.ceil() as usize;

        let background_total_in_term_count = background_pop
            .go_term_count
            .get(&taxon_id)
            .and_then(|tgc| tgc.get(&term_id))
            .copied()
            .unwrap_or(0);

        let total_significant_study_proteins_count = all_study_proteins.len();
        
        let total_background_proteins_count = background_pop
            .taxon_protein_count
            .get(&taxon_id)
            .copied()
            .unwrap();

        let contingency_table = create_contingency_table(
            study_sig_in_term_count,
            background_total_in_term_count,
            total_significant_study_proteins_count,
            total_background_proteins_count,
        );

        GOTermResults {
            log_odds_ratio: calculate_log_odds_ratio(&contingency_table),
            p_value: calculate_p_value(&contingency_table, self.test_type),
            contingency_table,
            variance: calculate_variance(&contingency_table),
        }
    }
}