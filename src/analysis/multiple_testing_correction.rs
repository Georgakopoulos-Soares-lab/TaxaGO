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

    fn extract_p_value(&self) -> f64 { self.p_value }
    fn extract_log_odds_ratio(&self) -> f64 { self.log_odds_ratio }

    fn with_adjusted_p_value(&self, new_p_value: f64) -> Self {
        Self {
            log_odds_ratio: self.log_odds_ratio,
            p_value: new_p_value,
            contingency_table: self.contingency_table,
            variance: self.variance,
        }
    }
}

impl PValueAdjustable for TaxonomyGOResult {
    type Key = String;

    fn extract_p_value(&self) -> f64 { self.p_value }
    fn extract_log_odds_ratio(&self) -> f64 { self.log_odds_ratio }

    fn with_adjusted_p_value(&self, new_p_value: f64) -> Self {
        Self {
            log_odds_ratio: self.log_odds_ratio,
            p_value: new_p_value,
            species_number: self.species_number,
        }
    }
}

fn adjust_p_values_grouped<T>(
    results: &FxHashMap<T::Key, FxHashMap<u32, T>>,
    method: AdjustmentMethod,
    significance_threshold: Option<f64>,
    log_odds_ratio_threshold: f64,
) -> FxHashMap<T::Key, FxHashMap<u32, T>>
where
    T: PValueAdjustable + Clone,
    T::Key: Clone + Hash + Eq,
{
    if matches!(method, AdjustmentMethod::None) {
        let mut filtered = FxHashMap::default();
        for (key, go_terms) in results {
            for (&go_id, res) in go_terms {
                let pass_p = significance_threshold.map_or(true, |thr| res.extract_p_value() <= thr);
                let pass_es = res.extract_log_odds_ratio() >= log_odds_ratio_threshold;
                if pass_p && pass_es {
                    filtered.entry(key.clone())
                        .or_insert_with(FxHashMap::default)
                        .insert(go_id, res.clone());
                }
            }
        }
        return filtered;
    }

    // Group by species/taxonomy key only (not by namespace)
    let mut groups: FxHashMap<T::Key, Vec<(u32, T)>> = FxHashMap::default();
    for (key, go_terms) in results.iter() {
        for (&go_id, res) in go_terms.iter() {
            groups.entry(key.clone())
                  .or_default()
                  .push((go_id, res.clone()));
        }
    }

    let mut out: FxHashMap<T::Key, FxHashMap<u32, T>> = FxHashMap::default();
    for (key, rows) in groups {
        // Collect all p-values for this species/taxonomy
        let pvals: Vec<f64> = rows.iter().map(|(_, r)| r.extract_p_value()).collect();
        
        // Apply multiple testing correction across all GO terms for this species
        let adj = if let Some(proc_) = method.to_procedure() {
            adjust(&pvals, proc_)
        } else {
            pvals
        };

        // Filter results based on adjusted p-values and effect size
        for ((go_id, res), q) in rows.into_iter().zip(adj.into_iter()) {
            let pass_p = significance_threshold.map_or(true, |thr| q <= thr);
            let pass_es = res.extract_log_odds_ratio() >= log_odds_ratio_threshold;
            if pass_p && pass_es {
                out.entry(key.clone())
                    .or_insert_with(FxHashMap::default)
                    .insert(go_id, res.with_adjusted_p_value(q));
            }
        }
    }
    out
}

pub fn adjust_species_p_values(
    results: &SpeciesResults,
    adjustment_method: AdjustmentMethod,
    significance_threshold: Option<f64>,
    log_odds_ratio_threshold: f64,
) -> SpeciesResults {
    println!("Adjusting single taxon p-values using method: {:?} (grouped by species)\n", adjustment_method);
    adjust_p_values_grouped(
        results,
        adjustment_method,
        significance_threshold,
        log_odds_ratio_threshold,
    )
}

pub fn adjust_taxonomy_p_values(
    results: &TaxonomyResults,
    adjustment_method: AdjustmentMethod,
    significance_threshold: Option<f64>,
    log_odds_ratio_threshold: f64,
    level: &String,
) -> TaxonomyResults {
    println!("Adjusting p-values at {} level using method: {:?} (grouped by taxonomy)\n", level, adjustment_method);
    adjust_p_values_grouped(
        results,
        adjustment_method,
        significance_threshold,
        log_odds_ratio_threshold,
    )
}