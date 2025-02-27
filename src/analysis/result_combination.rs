use std::collections::HashMap;
use std::f64::consts::PI;
use crate::analysis::enrichment_analysis::GOTermResults;
use crate::analysis::multiple_testing_correction::TaxonomyGOResult;
use statrs::distribution::{Normal, ContinuousCDF};

#[derive(Copy, Clone)]
pub enum PValueCombinationMethod {
    Fisher,
    Stouffer,
    Cauchy,
}
pub fn group_results_by_taxonomy(
    family_taxa: &HashMap<String, Vec<u32>>,
    fisher_results: &HashMap<u32, HashMap<u32, GOTermResults>>,
    threshold: f64,
) -> HashMap<String, HashMap<u32, HashMap<u32, (f64, f64, [usize; 4], f64)>>> {
    let mut result = HashMap::new();
    let mut go_term_counts: HashMap<String, HashMap<u32, usize>> = HashMap::new();
    let mut actual_species_counts: HashMap<String, usize> = HashMap::new();

    count_species_and_go_terms(family_taxa, fisher_results, &mut actual_species_counts, &mut go_term_counts);

    process_and_filter_results(family_taxa, fisher_results, threshold, &actual_species_counts, &go_term_counts, &mut result);

    result
}

fn count_species_and_go_terms(
    family_taxa: &HashMap<String, Vec<u32>>,
    fisher_results: &HashMap<u32, HashMap<u32, GOTermResults>>,
    actual_species_counts: &mut HashMap<String, usize>,
    go_term_counts: &mut HashMap<String, HashMap<u32, usize>>,
) {
    for (family, taxa) in family_taxa {
        let species_count = taxa.len();
        actual_species_counts.insert(family.clone(), species_count);

        let family_go_terms = go_term_counts.entry(family.clone()).or_default();

        for taxon_id in taxa {
            if let Some(go_terms) = fisher_results.get(taxon_id) {
                for go_term_id in go_terms.keys() {
                    *family_go_terms.entry(*go_term_id).or_default() += 1;
                }
            }
        }
    }
}

fn process_and_filter_results(
    family_taxa: &HashMap<String, Vec<u32>>,
    fisher_results: &HashMap<u32, HashMap<u32, GOTermResults>>,
    threshold: f64,
    actual_species_counts: &HashMap<String, usize>,
    go_term_counts: &HashMap<String, HashMap<u32, usize>>,
    result: &mut HashMap<String, HashMap<u32, HashMap<u32, (f64, f64, [usize; 4], f64)>>>,
) {
    for (family, taxa) in family_taxa {
        let species_count = *actual_species_counts.get(family).unwrap_or(&0);

        if species_count == 0 {
            continue;
        }

        let empty_map = HashMap::new();
        let family_go_terms = go_term_counts.get(family).unwrap_or(&empty_map);
        let mut taxa_map = HashMap::new();

        for taxon_id in taxa {
            if let Some(go_terms) = fisher_results.get(taxon_id) {
                let mut go_terms_with_variance = HashMap::new();

                for (go_term_id, go_term_result) in go_terms {
                    let go_term_count = family_go_terms.get(go_term_id).unwrap_or(&0);
                    let prevalence = *go_term_count as f64 / species_count as f64;

                    if prevalence >= threshold {
                        let variance = calculate_variance(&go_term_result.contingency_table);
                        go_terms_with_variance.insert(
                            *go_term_id, 
                            (
                                go_term_result.log_odds_ratio,
                                go_term_result.p_value,
                                go_term_result.contingency_table,
                                variance
                            )
                        );
                    }
                }

                if !go_terms_with_variance.is_empty() {
                    taxa_map.insert(*taxon_id, go_terms_with_variance);
                }
            }
        }

        if !taxa_map.is_empty() {
            result.insert(family.clone(), taxa_map);
        }
    }
}

fn calculate_variance(contingency: &[usize; 4]) -> f64 {
    if contingency.iter().all(|&x| x > 0) {
        1.0 / contingency[0] as f64 +
        1.0 / contingency[1] as f64 +
        1.0 / contingency[2] as f64 +
        1.0 / contingency[3] as f64
    } else {
        f64::INFINITY
    }
}

pub fn combine_taxonomic_results(
    organized_results: &HashMap<String, HashMap<u32, HashMap<u32, (f64, f64, [usize; 4], f64)>>>,
    tolerance: f64,
    max_iterations: usize,
    combination_method: &PValueCombinationMethod,
) -> HashMap<String, HashMap<u32, TaxonomyGOResult>> {
    let mut final_estimates = HashMap::new();

    for (family, taxa_results) in organized_results {
        let mut go_term_estimates = HashMap::new();
        let total_species = taxa_results.len();

        let all_go_terms: Vec<u32> = taxa_results.values().flat_map(|go_maps| go_maps.keys().cloned()).collect();
        let unique_go_terms: Vec<u32> = {
            let mut temp = all_go_terms.clone();
            temp.sort_unstable();
            temp.dedup();
            temp
        };

        for go_term in unique_go_terms {
            let mut effect_sizes = Vec::new();
            let mut variances = Vec::new();
            let mut p_values = Vec::new();

            let species_with_go = taxa_results.values()
                .filter(|taxa_map| taxa_map.contains_key(&go_term))
                .count();
            
            let species_percentage = if total_species > 0 {
                (species_with_go as f64 / total_species as f64) * 100.0
            } else {
                0.0
            };

            for taxa_map in taxa_results.values() {
                if let Some(&(lor, p, _, var)) = taxa_map.get(&go_term) {
                    effect_sizes.push(lor);
                    variances.push(var);
                    p_values.push(p);
                }
            }

            match effect_sizes.len() {
                0 => continue,
                1 => {
                    let tau_squared = 0.0;
                    let mean_effect = effect_sizes[0];
                    let p_value = p_values[0];
                    go_term_estimates.insert(go_term, TaxonomyGOResult {
                        log_odds_ratio: mean_effect,
                        p_value,
                        tau_squared,
                        species_percentage,
                        species_count: species_with_go,
                        total_species,
                    });
                }
                _ => {
                    let (tau_squared, mean_effect) = paule_mandel_estimator(&effect_sizes, &variances, tolerance, max_iterations);
                    let combined_p = combine_pvalues(&p_values, &combination_method);
                    go_term_estimates.insert(go_term, TaxonomyGOResult {
                        log_odds_ratio: mean_effect,
                        p_value: combined_p,
                        tau_squared,
                        species_percentage,
                        species_count: species_with_go,
                        total_species,
                    });
                }
            }
        }

        final_estimates.insert(family.clone(), go_term_estimates);
    }

    final_estimates
}

pub fn paule_mandel_estimator(
    effect_sizes: &[f64],
    variances: &[f64],
    tolerance: f64,
    max_iterations: usize,
) -> (f64, f64) {
    assert_eq!(effect_sizes.len(), variances.len(), "Effect sizes and variances must have the same length");

    let mut tau_squared = 0.0;
    let mut iteration = 0;

    loop {
        let weights: Vec<f64> = variances.iter().map(|&v| 1.0 / (v + tau_squared)).collect();
        let weighted_mean: f64 = effect_sizes.iter().zip(&weights).map(|(&y, &w)| y * w).sum::<f64>() / weights.iter().sum::<f64>();

        let numerator: f64 = effect_sizes.iter().zip(&weights).map(|(&y, &w)| w * w * (y - weighted_mean).powi(2)).sum::<f64>() - weights.iter().sum::<f64>();
        let denominator: f64 = weights.iter().sum::<f64>();

        let new_tau_squared = (numerator / denominator).max(0.0);

        if (new_tau_squared - tau_squared).abs() < tolerance || iteration >= max_iterations {
            let final_weights: Vec<f64> = variances.iter().map(|&v| 1.0 / (v + new_tau_squared)).collect();
            let final_mean: f64 = effect_sizes.iter().zip(&final_weights).map(|(&y, &w)| y * w).sum::<f64>() / final_weights.iter().sum::<f64>();

            return (new_tau_squared, final_mean);
        }

        tau_squared = new_tau_squared;
        iteration += 1;
    }
}

pub fn combine_pvalues(
    p_values: &[f64], 
    method: &PValueCombinationMethod
) -> f64 {
    const EPSILON: f64 = 1e-300;

    if p_values.is_empty() {
        return 1.0;
    }

    match method {
        PValueCombinationMethod::Cauchy => {
            let t = p_values
                .iter()
                .map(|&p| (0.5 - p.clamp(EPSILON, 1.0 - EPSILON)) * PI)
                .map(|term| term.tan())
                .sum::<f64>();

            let abs_t = t.abs();
            let survival = 0.5 - (abs_t.atan() / PI);
            let combined_p = 2.0 * survival;

            combined_p.clamp(0.0, 1.0)
        },
        PValueCombinationMethod::Stouffer => {
            let normal = Normal::new(0.0, 1.0).unwrap();
            let n = p_values.len() as f64;
            
            let z_sum = p_values
                .iter()
                .map(|&p| p.clamp(EPSILON, 1.0 - EPSILON))
                .map(|p| normal.inverse_cdf(1.0 - p))
                .sum::<f64>();
            
            let z_combined = z_sum / n.sqrt();
            
            1.0 - normal.cdf(z_combined)
        },
        PValueCombinationMethod::Fisher => {0.5} // to update 
    }
}