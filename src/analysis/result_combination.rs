use rustc_hash::FxHashMap;
use crate::{analysis::enrichment_analysis::*, parsers::background_parser::{GOTermID, TaxonID}};

pub fn group_results_by_taxonomy(
    family_taxa: &FxHashMap<String, Vec<TaxonID>>,
    fisher_results: &FxHashMap<TaxonID, FxHashMap<GOTermID, GOTermResults>>,
    threshold: f64,
) -> FxHashMap<String, FxHashMap<TaxonID, FxHashMap<GOTermID, GOTermResults>>> {
    let mut result = FxHashMap::default();
    let mut go_term_counts: FxHashMap<String, FxHashMap<GOTermID, usize>> = FxHashMap::default();
    let mut actual_species_counts: FxHashMap<String, usize> = FxHashMap::default();

    count_species_and_go_terms(family_taxa, fisher_results, &mut actual_species_counts, &mut go_term_counts);

    process_and_filter_results(family_taxa, fisher_results, threshold, &actual_species_counts, &go_term_counts, &mut result);

    result
}

fn count_species_and_go_terms(
    family_taxa: &FxHashMap<String, Vec<u32>>,
    fisher_results: &FxHashMap<u32, FxHashMap<u32, GOTermResults>>,
    actual_species_counts: &mut FxHashMap<String, usize>,
    go_term_counts: &mut FxHashMap<String, FxHashMap<u32, usize>>,
) {
    family_taxa.iter().for_each(|(family, taxa)| {
        actual_species_counts.insert(family.to_string(), taxa.len());
        
        let family_go_terms = go_term_counts.entry(family.clone()).or_default();
        
        taxa.iter()
            .filter_map(|taxon_id| fisher_results.get(taxon_id))
            .flat_map(|go_terms| go_terms.keys())
            .for_each(|go_term_id| {
                *family_go_terms.entry(*go_term_id).or_default() += 1;
            });
    });
}

fn process_and_filter_results(
    family_taxa: &FxHashMap<String, Vec<u32>>,
    fisher_results: &FxHashMap<u32, FxHashMap<u32, GOTermResults>>,
    threshold: f64,
    actual_species_counts: &FxHashMap<String, usize>,
    go_term_counts: &FxHashMap<String, FxHashMap<u32, usize>>,
    result: &mut FxHashMap<String, FxHashMap<TaxonID, FxHashMap<GOTermID, GOTermResults>>>,
) {
    family_taxa.iter().for_each(|(family, taxa)| {
        let species_count = *actual_species_counts.get(family).unwrap_or(&0);

        if species_count == 0 {
            return;
        }

        let empty_map = FxHashMap::default();
        let family_go_terms = go_term_counts.get(family).unwrap_or(&empty_map);
        
        let taxa_map: FxHashMap<TaxonID, FxHashMap<GOTermID, GOTermResults>> = taxa.iter()
            .filter_map(|taxon_id| {
                fisher_results.get(taxon_id).map(|go_terms| {
                    let go_terms_with_variance: FxHashMap<GOTermID, GOTermResults> = go_terms.iter()
                        .filter_map(|(go_term_id, go_term_result)| {
                            let go_term_count = family_go_terms.get(go_term_id).unwrap_or(&0);
                            let prevalence = *go_term_count as f64 / species_count as f64;
                            
                            if prevalence >= threshold {
                                Some((*go_term_id, go_term_result.clone()))
                            } else {
                                None
                            }
                        })
                        .collect();
                    
                    if go_terms_with_variance.is_empty() {
                        None
                    } else {
                        Some((*taxon_id, go_terms_with_variance))
                    }
                })
            })
            .flatten()
            .collect();
        
        if !taxa_map.is_empty() {
            result.insert(family.clone(), taxa_map);
        }
    });
}