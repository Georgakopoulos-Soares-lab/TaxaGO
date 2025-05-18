use rustc_hash::FxHashMap;
use adjustp::{adjust, Procedure};
use crate::{
    parsers::background_parser::*,
    analysis::enrichment_analysis::*,
    analysis::phylogenetic_meta_analysis::*
};
use clap::ValueEnum;
use std::hash::Hash;

type SpeciesResults = FxHashMap<TaxonID, FxHashMap<GOTermID, GOTermResults>>;
type TaxonomyResults = FxHashMap<String, FxHashMap<GOTermID, TaxonomyGOResult>>;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum AdjustmentMethod {
    None,
    Bonferroni,
    BenjaminiHochberg,
    BenjaminiYekutieli,
}

impl AdjustmentMethod {
    fn to_procedure(&self) -> Option<Procedure> {
        match self {
            AdjustmentMethod::None => None,
            AdjustmentMethod::Bonferroni => Some(Procedure::Bonferroni),
            AdjustmentMethod::BenjaminiHochberg => Some(Procedure::BenjaminiHochberg),
            AdjustmentMethod::BenjaminiYekutieli => Some(Procedure::BenjaminiYekutieli),
        }
    }
}

trait PValueAdjustable {
    type Key;
    
    fn extract_p_value(&self) -> f64;
    fn extract_log_odds_ratio(&self) -> f64;
    fn with_adjusted_p_value(&self, new_p_value: f64) -> Self;
}

impl PValueAdjustable for GOTermResults {
    type Key = u32;
    
    fn extract_p_value(&self) -> f64 {
        self.p_value
    }

    fn extract_log_odds_ratio(&self) -> f64 {
        self.log_odds_ratio
    }
    
    fn with_adjusted_p_value(&self, new_p_value: f64) -> Self {
        Self {
            log_odds_ratio: self.log_odds_ratio,
            p_value: new_p_value,
            contingency_table: self.contingency_table,
            variance: self.variance
        }
    }
}

impl PValueAdjustable for TaxonomyGOResult {
    type Key = String;
    
    fn extract_p_value(&self) -> f64 {
        self.p_value
    }

    fn extract_log_odds_ratio(&self) -> f64 {
        self.log_odds_ratio
    }

    
    fn with_adjusted_p_value(&self, new_p_value: f64) -> Self {
        Self {
            log_odds_ratio: self.log_odds_ratio,
            p_value: new_p_value,
            species_number: self.species_number
        }
    }
}

fn adjust_p_values<T: PValueAdjustable + Clone>(
    results: &FxHashMap<T::Key, FxHashMap<u32, T>>,
    method: AdjustmentMethod,
    significance_threshold: Option<f64>,
    log_odds_ratio_threshold: f64
) -> FxHashMap<T::Key, FxHashMap<u32, T>>
where
    T::Key: Clone + Hash + Eq,
{
    if matches!(method, AdjustmentMethod::None) {
        println!("No p-value adjustment performed. Applying p-value and log odds ratio filter to original values.\n");
        let mut filtered_results = FxHashMap::default();
        for (key, go_terms) in results {
            for (&go_id, result_item) in go_terms {
                let p_value_passes = significance_threshold.map_or(true, |thresh| result_item.extract_p_value() <= thresh);
                let log_odds_passes = result_item.extract_log_odds_ratio() >= log_odds_ratio_threshold;

                if p_value_passes && log_odds_passes {
                    filtered_results
                        .entry(key.clone())
                        .or_insert_with(FxHashMap::default)
                        .insert(go_id, result_item.clone());
                }
            }
        }
        return filtered_results;
    }
    
    let mut all_pvalues = Vec::new();
    let mut pvalue_mapping = Vec::new();
    
    for (key, go_terms) in results {
        for (&go_id, result) in go_terms {
            all_pvalues.push(result.extract_p_value());
            pvalue_mapping.push((key.clone(), go_id, result.clone()));
        }
    }
    
    let adjusted_pvalues = if let Some(proc) = method.to_procedure() {
        adjust(&all_pvalues, proc)
    } else {
        all_pvalues
    };
    
    let mut adjusted_results = FxHashMap::default();
    for (index, (key, go_id, result)) in pvalue_mapping.into_iter().enumerate() {
        let adjusted_p = adjusted_pvalues[index];
        let p_value_passes = significance_threshold.map_or(true, |threshold| adjusted_p <= threshold);
        let log_odds_passes = result.extract_log_odds_ratio() >= log_odds_ratio_threshold;

        if p_value_passes && log_odds_passes {
            adjusted_results
                .entry(key)
                .or_insert_with(FxHashMap::default)
                .insert(go_id, result.with_adjusted_p_value(adjusted_p));
        }
    }
    
    adjusted_results
}

pub fn adjust_species_p_values(
    results: &SpeciesResults,
    adjustment_method: AdjustmentMethod,
    significance_threshold: Option<f64>,
    log_odds_ratio_threshold: f64,
) -> SpeciesResults {
    println!("Adjusting single taxon p-values using method: {:?}\n", adjustment_method);
    adjust_p_values(
        results, 
        adjustment_method, 
        significance_threshold,
        log_odds_ratio_threshold
    )
}

pub fn adjust_taxonomy_p_values(
    results: &TaxonomyResults,
    adjustment_method: AdjustmentMethod,
    significance_threshold: Option<f64>,
    log_odds_ratio_threshold: f64,
    level: &String,
) -> TaxonomyResults {
    println!("Adjusting p-values at {} level using method: {:?}\n", level, adjustment_method);
    adjust_p_values(
        results, 
        adjustment_method, 
        significance_threshold,
        log_odds_ratio_threshold
    )
}