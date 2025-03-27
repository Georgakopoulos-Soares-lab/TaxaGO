use rustc_hash::{FxHashMap, FxHashSet};
use crate::parsers::{
    background_parser::*,
    study_parser::*,
    obo_parser::*
};
use crate::analysis::
    enrichment_analysis::*;
use rayon::prelude::*;

impl EnrichmentAnalysis{
    pub fn elim_analysis(
        &self,
        taxon_ids: &FxHashSet<TaxonID>, 
        significance_threshold: f64,         
        study_pop: &StudyPop,
        background_pop: &BackgroundPop,
        level_to_go_term: &LevelToTerms
    ) -> FxHashMap<TaxonID, FxHashMap<GOTermID, GOTermResults>> {
        
        let max_level = match level_to_go_term.keys().max() {
            Some(&level) => level,
            None => return FxHashMap::default()
        };
    
        taxon_ids.par_iter()
            .map(|&taxon_id| {
        
                let taxon_study_go_term_proteins = match study_pop.go_term_to_protein_set.get(&taxon_id) {
                    Some(proteins) => proteins,
                    None => return (taxon_id, FxHashMap::default())
                };
                
                let taxon_study_total_count = match study_pop.taxon_protein_count.get(&taxon_id) {
                    Some(&count) => count,
                    None => return (taxon_id, FxHashMap::default())
                };
        
                let taxon_background_go_term_count = match background_pop.go_term_count.get(&taxon_id) {
                    Some(counts) => counts,
                    None => return (taxon_id, FxHashMap::default())
                };
                
                let taxon_background_total_count = match background_pop.taxon_protein_count.get(&taxon_id) {
                    Some(&count) => count,
                    None => return (taxon_id, FxHashMap::default())
                };
        
                let mut marked_proteins: FxHashSet<&Protein> = FxHashSet::default();
                let mut go_term_results = FxHashMap::with_capacity_and_hasher(
                    level_to_go_term.values().map(|v| v.len()).sum(),
                    rustc_hash::FxBuildHasher::default()
                );
                
                for level in (1..=max_level).rev() {
                    if let Some(go_terms) = level_to_go_term.get(&level) {
                        for &go_term in go_terms {
                            let original_study_proteins = match taxon_study_go_term_proteins.get(&go_term) {
                                Some(proteins) => proteins,
                                None => continue
                            };
                            
                            let background_counts = match taxon_background_go_term_count.get(&go_term) {
                                Some(&count) => count,
                                None => continue
                            };
                            
                            let unmarked_proteins: Vec<&Protein> = original_study_proteins
                                .iter()
                                .filter(|protein| !marked_proteins.contains(protein))
                                .collect();
                            
                            let study_counts = unmarked_proteins.len();
                
                            let results = analyze_single_go_term(
                                study_counts,
                                background_counts,
                                taxon_study_total_count,
                                taxon_background_total_count,
                                self.test_type
                            );
                
                            if results.p_value <= significance_threshold {
                                marked_proteins.extend(unmarked_proteins);
                            }
                
                            go_term_results.insert(go_term, results);
                        } 
                    }
                }
                
                (taxon_id, go_term_results)
            })
            .collect()
    }
}