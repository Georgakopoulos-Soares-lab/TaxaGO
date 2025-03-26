use std::collections::{HashMap, HashSet};
use crate::parsers::{
    background_parser::*,
    study_parser::*,
    obo_parser::*
};
use crate::analysis::
    enrichment_analysis::*;
// use rayon::prelude::*;


impl EnrichmentAnalysis{
    pub fn elim_analysis(
        &self,
        taxon_ids: &HashSet<TaxonID>, 
        significance_threshold: f64,         
        study_pop: &StudyPop,
        background_pop: &BackgroundPop,
        level_to_go_term: &LevelToTerms
    ) -> HashMap<TaxonID, HashMap<GOTermID, GOTermResults>> {
        
        let mut enrichment_results: HashMap<TaxonID, HashMap<GOTermID, GOTermResults>> = HashMap::with_capacity(taxon_ids.len());
       
        let max_level = *level_to_go_term.keys().max().unwrap();
    
        for taxon_id in taxon_ids {
            let taxon_study_go_term_proteins = study_pop.go_term_to_protein_set.get(taxon_id).unwrap();
            let taxon_study_total_count = study_pop.taxon_protein_count.get(taxon_id).unwrap();
    
            let taxon_background_go_term_count = background_pop.go_term_count.get(taxon_id).unwrap();
            let taxon_background_total_count = background_pop.taxon_protein_count.get(taxon_id).unwrap();
    
            let mut marked_proteins: HashSet<&Protein> = HashSet::new();
            let mut go_term_results: HashMap<GOTermID, GOTermResults> = HashMap::new();
            
            for level in (1..=max_level).rev() {
                let go_terms = level_to_go_term.get(&level).unwrap();
    
                for go_term in go_terms {
                    let original_study_proteins = taxon_study_go_term_proteins.get(go_term).unwrap();
                    let background_counts = *taxon_background_go_term_count.get(go_term).unwrap();
                    
                    let study_counts = original_study_proteins
                        .iter()
                        .filter(|protein| !marked_proteins.contains(protein))
                        .count();
    
                    let results = analyze_single_go_term(
                        study_counts,
                        background_counts,
                        *taxon_study_total_count,
                        *taxon_background_total_count,
                        self.test_type
                    );
    
                    if results.p_value <= significance_threshold {
                        let significant_proteins: Vec<&Protein> = original_study_proteins
                            .iter()
                            .filter(|protein| !marked_proteins.contains(protein))
                            .collect();
                            
                        marked_proteins.extend(significant_proteins);
                    }
    
                    go_term_results.insert(*go_term, results);
                } 
            }
            enrichment_results.insert(*taxon_id, go_term_results);
        }
        
        enrichment_results
    }
}
