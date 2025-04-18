use rustc_hash::FxHashMap;
use adjustp::{adjust, Procedure};
use crate::analysis::enrichment_analysis::GOTermResults;

#[derive(Debug, Clone)]
pub struct TaxonomyGOResult {
    pub log_odds_ratio: f64,
    pub p_value: f64,
    pub tau_squared: f64,
    pub species_percentage: f64,
    pub species_count: usize,
    pub total_species: usize,
}

type SpeciesResults = FxHashMap<u32, FxHashMap<u32, GOTermResults>>;
type TaxonomyResults = FxHashMap<String, FxHashMap<u32, TaxonomyGOResult>>;

#[derive(Debug)]
pub enum AdjustmentMethod {
    None,
    Bonferroni,
    BenjaminiHochberg,
    BenjaminiYekutieli,
}

impl From<&str> for AdjustmentMethod {
    fn from(s: &str) -> Self {
        match s.to_lowercase().as_str() {
            "none" => AdjustmentMethod::None,
            "bonferroni" => AdjustmentMethod::Bonferroni,
            "bh" => AdjustmentMethod::BenjaminiHochberg,
            "by" => AdjustmentMethod::BenjaminiYekutieli,
            _ => AdjustmentMethod::None,
        }
    }
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
    fn with_adjusted_p_value(&self, new_p_value: f64) -> Self;
}

impl PValueAdjustable for GOTermResults {
    type Key = u32;
    
    fn extract_p_value(&self) -> f64 {
        self.p_value
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
    
    fn with_adjusted_p_value(&self, new_p_value: f64) -> Self {
        Self {
            log_odds_ratio: self.log_odds_ratio,
            p_value: new_p_value,
            tau_squared: self.tau_squared,
            species_percentage: self.species_percentage,
            species_count: self.species_count,
            total_species: self.total_species,
        }
    }
}

fn adjust_p_values<T: PValueAdjustable + Clone>(
    results: &FxHashMap<T::Key, FxHashMap<u32, T>>,
    method: AdjustmentMethod,
    significance_threshold: Option<f64>,
) -> FxHashMap<T::Key, FxHashMap<u32, T>>
where
    T::Key: Clone + std::hash::Hash + Eq,
{
    if matches!(method, AdjustmentMethod::None) {
        println!("No p-value adjustment performed\n");
        return results.clone();
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
        
        if significance_threshold.map_or(true, |threshold| adjusted_p <= threshold) {
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
    adjustment_method: &str,
    significance_threshold: Option<f64>,
) -> SpeciesResults {
    println!("Adjusting single taxon p-values using method: {:?}\n", adjustment_method);
    adjust_p_values(results, adjustment_method.into(), significance_threshold)
}

pub fn adjust_taxonomy_p_values(
    results: &TaxonomyResults,
    adjustment_method: &str,
    significance_threshold: Option<f64>,
    level: &String,
) -> TaxonomyResults {
    println!("Adjusting p-values at {} level using method: {:?}\n", level, adjustment_method);
    adjust_p_values(results, adjustment_method.into(), significance_threshold)
}