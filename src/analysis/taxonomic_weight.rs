use rustc_hash::FxHashMap;
use crate::parsers::background_parser::*;
use crate::analysis::multiple_testing_correction::*;
use rayon::prelude::*;

/// Calculate taxonomic distance between two species
/// 
/// Returns a value between 0.0 (identical) and 1.0 (maximally distant)
/// This calculation is based on shared taxonomic ranks below the specified level
fn calculate_taxonomic_distance(
    lineage_map: &FxHashMap<u32, Vec<String>>,
    taxon_a: TaxonID,
    taxon_b: TaxonID,
    level_index: usize 
) -> f64 {
    // Get lineages for both taxa
    let lineage_a = match lineage_map.get(&taxon_a) {
        Some(l) => l,
        None => return 1.0, // Maximum distance if lineage unknown
    };
    
    let lineage_b = match lineage_map.get(&taxon_b) {
        Some(l) => l,
        None => return 1.0,
    };
    
    // Count shared ranks below the specified level
    let mut shared_ranks: usize = 0;
    let max_ranks: usize = level_index; // Only consider ranks below the specified level
    
    for i in 0..max_ranks {
        if i < lineage_a.len() && i < lineage_b.len() && lineage_a[i] == lineage_b[i] {
            shared_ranks += 1;
        } else {
            break; // Stop at first difference
        }
    }
    
    // If looking at the most specific level (e.g., genus), distance is 0
    if max_ranks == 0 {
        return 0.0;
    }
    
    // Convert shared ranks to a distance
    1.0 - (shared_ranks as f64 / max_ranks as f64)
}

/// Determine the index corresponding to a taxonomic level
fn get_taxonomic_level_index(level: &str) -> Option<usize> {
    match level.to_lowercase().as_str() {
        "genus" => Some(0),
        "family" => Some(1),
        "order" => Some(2),
        "class" => Some(3),
        "phylum" => Some(4),
        "kingdom" => Some(5),
        "superkingdom" => Some(6),
        _ => None,
    }
}

/// Build distance matrices for each taxonomic group
/// 
/// Returns a map of taxonomic group names to their distance matrices
fn build_distance_matrices_by_group(
    lineage_map: &FxHashMap<u32, Vec<String>>,
    family_taxa: &FxHashMap<String, Vec<u32>>,
    taxonomic_level: &str
) -> FxHashMap<String, FxHashMap<(TaxonID, TaxonID), f64>> {
    // Get the index for the taxonomic level
    let level_index = match get_taxonomic_level_index(taxonomic_level) {
        Some(idx) => idx,
        None => {
            eprintln!("Invalid taxonomic level: {} \n", taxonomic_level);
            return FxHashMap::default();
        }
    };
    
    // Create distance matrices for each group
    let mut group_matrices = FxHashMap::default();
    
    for (group_name, taxa) in family_taxa {
        let n = taxa.len();
        
        // Generate all unique taxon pairs
        let mut pairs = Vec::with_capacity(n * (n + 1) / 2);
        for i in 0..n {
            for j in i..n {
                pairs.push((taxa[i], taxa[j]));
            }
        }
        
        // Calculate distances in parallel
        let distances: Vec<((TaxonID, TaxonID), f64)> = pairs
            .par_iter()
            .map(|&(taxon_a, taxon_b)| {
                let distance = calculate_taxonomic_distance(
                    lineage_map,
                    taxon_a,
                    taxon_b,
                    level_index
                );
                
                ((taxon_a, taxon_b), distance)
            })
            .collect();
        
        // Build the distance matrix from results
        let mut distance_matrix = FxHashMap::default();
        for ((taxon_a, taxon_b), distance) in distances {
            distance_matrix.insert((taxon_a, taxon_b), distance);
            
            // Add symmetric pair if not the same taxon
            if taxon_a != taxon_b {
                distance_matrix.insert((taxon_b, taxon_a), distance);
            }
        }
        
        group_matrices.insert(group_name.clone(), distance_matrix);
    }
    
    group_matrices
}

/// Calculate taxonomic weights for species within a group
/// 
/// Returns a map of taxon IDs to their weights
fn calculate_taxonomic_weights(
    distance_matrix: &FxHashMap<(TaxonID, TaxonID), f64>,
    taxa: &[TaxonID]
) -> FxHashMap<TaxonID, f64> {
    let weights: Vec<(TaxonID, f64)> = taxa.par_iter()
    .map(|&taxon_id| {
        let mut total_distance = 0.0;
        let mut count = 0;
        
        for &other_id in taxa {
            if taxon_id != other_id {
                if let Some(&distance) = distance_matrix.get(&(taxon_id, other_id)) {
                    total_distance += distance;
                    count += 1;
                }
            }
        }
        
        let avg_distance = if count > 0 { total_distance / count as f64 } else { 0.0 };
        let weight = 0.5 + avg_distance;
        
        (taxon_id, weight)
    })
    .collect();

    // Convert to hashmap
    let mut taxon_weights = FxHashMap::default();
    for (id, weight) in weights {
        taxon_weights.insert(id, weight);
    }

    // Normalize weights (must be done sequentially after all weights are calculated)
    let total_weight: f64 = taxon_weights.values().sum();
    let norm_factor = taxa.len() as f64 / total_weight;

    for weight in taxon_weights.values_mut() {
        *weight *= norm_factor;
    }

    taxon_weights
    }

/// Taxonomically weighted Paule-Mandel estimator
/// 
/// Incorporates both statistical weights (based on variance) and 
/// biological weights (based on taxonomic uniqueness)
pub fn taxonomic_paule_mandel_estimator(
    effect_sizes: &[f64],
    variances: &[f64],
    bio_weights: &[f64],
    tolerance: f64,
    max_iterations: usize,
) -> (f64, f64) {
    assert_eq!(effect_sizes.len(), variances.len(), "Effect sizes and variances must have the same length");
    assert_eq!(effect_sizes.len(), bio_weights.len(), "Effect sizes and biological weights must have the same length");

    let mut tau_squared = 0.0;
    let mut iteration = 0;

    loop {
        // Calculate combined weights (statistical × biological)
        let combined_weights: Vec<f64> = variances.iter()
            .zip(bio_weights)
            .map(|(&v, &b)| b / (v + tau_squared))
            .collect();
            
        // Calculate weighted mean effect size
        let sum_weights = combined_weights.iter().sum::<f64>();
        let weighted_mean = effect_sizes.iter()
            .zip(&combined_weights)
            .map(|(&y, &w)| y * w)
            .sum::<f64>() / sum_weights;
        
        // Calculate numerator and denominator for tau² estimate
        let numerator = effect_sizes.iter()
            .zip(&combined_weights)
            .map(|(&y, &w)| w * w * (y - weighted_mean).powi(2))
            .sum::<f64>() - sum_weights;
            
        let denominator = combined_weights.iter()
            .map(|&w| w * w)
            .sum::<f64>();
        
        // Calculate new tau² estimate
        let new_tau_squared = (numerator / denominator).max(0.0);
        
        // Check convergence
        if (new_tau_squared - tau_squared).abs() < tolerance || iteration >= max_iterations {
            // Calculate final weights with converged tau²
            let final_weights: Vec<f64> = variances.iter()
                .zip(bio_weights)
                .map(|(&v, &b)| b / (v + new_tau_squared))
                .collect();
                
            // Calculate final weighted mean
            let final_mean = effect_sizes.iter()
                .zip(&final_weights)
                .map(|(&y, &w)| y * w)
                .sum::<f64>() / final_weights.iter().sum::<f64>();
                
            return (new_tau_squared, final_mean);
        }
        
        tau_squared = new_tau_squared;
        iteration += 1;
    }
}

/// Combine p-values with taxonomic weights
/// 
/// Uses a weighted variant of the Simes procedure for multiple testing correction
pub fn combine_weighted_pvalues(
    p_values: &[f64],
    weights: &[f64]
) -> f64 {
    if p_values.is_empty() {
        return 1.0;
    }
    
    // Create (p-value, weight) pairs and sort by p-value
    let mut pairs: Vec<(f64, f64)> = p_values.iter()
        .zip(weights.iter())
        .map(|(&p, &w)| (p, w))
        .collect();
    
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    
    // Apply weighted Simes procedure
    let total_weight: f64 = weights.iter().sum();
    let mut min_adjusted_p: f64 = 1.0;
    let mut cumulative_weight: f64 = 0.0;
    
    for (p, weight) in pairs {
        cumulative_weight += weight;
        let weighted_p = p * total_weight / cumulative_weight;
        min_adjusted_p = min_adjusted_p.min(weighted_p);
    }
    
    // Ensure result is in valid p-value range
    f64::min(1.0, f64::max(0.0, min_adjusted_p))
}

/// Combine GO term enrichment results with taxonomic weighting
/// 
/// This function enhances the standard result combination with biologically
/// meaningful weights based on taxonomic relationships
pub fn combine_taxonomic_results_with_weighting(
    organized_results: &FxHashMap<String, FxHashMap<u32, FxHashMap<u32, (f64, f64, [usize; 4], f64)>>>,
    lineage_map: &FxHashMap<u32, Vec<String>>,
    family_taxa: &FxHashMap<String, Vec<u32>>,
    taxonomic_level: &str,
    tolerance: f64,
    max_iterations: usize,
) -> FxHashMap<String, FxHashMap<u32, TaxonomyGOResult>> {
    // Build distance matrices for each taxonomic group
    let group_distance_matrices = build_distance_matrices_by_group(
        lineage_map, 
        family_taxa, 
        taxonomic_level
    );
    
    // Final results container
    let mut final_estimates = FxHashMap::default();

    // Process each taxonomic group
    for (group, taxa_results) in organized_results {
        let mut go_term_estimates = FxHashMap::default();
        
        // Get taxa for this group
        let taxa = match family_taxa.get(group) {
            Some(t) => t,
            None => continue, // Skip if no taxa found
        };
        
        // Get distance matrix for this group
        let distance_matrix = match group_distance_matrices.get(group) {
            Some(m) => m,
            None => continue, // Skip if no distance matrix found
        };
        
        // Calculate taxonomic weights for this group
        let taxonomic_weights = calculate_taxonomic_weights(distance_matrix, taxa);
        
        // Find all GO terms for this group
        let total_species = taxa_results.len();
        let all_go_terms: Vec<u32> = taxa_results.values()
            .flat_map(|go_maps| go_maps.keys().cloned())
            .collect();
            
        let unique_go_terms: Vec<u32> = {
            let mut temp = all_go_terms.clone();
            temp.sort_unstable();
            temp.dedup();
            temp
        };

        // Process each GO term
        for go_term in unique_go_terms {
            let mut effect_sizes = Vec::new();
            let mut variances = Vec::new();
            let mut bio_weights = Vec::new();
            let mut p_values = Vec::new();
            let mut taxon_ids = Vec::new();

            // Count species with this GO term
            let species_with_go = taxa_results.iter()
                .filter(|(taxon_id, go_map)| go_map.contains_key(&go_term))
                .count();
                
            let species_percentage = if total_species > 0 {
                (species_with_go as f64 / total_species as f64) * 100.0
            } else {
                0.0
            };

            // Collect data for this GO term across all species
            for (taxon_id, go_map) in taxa_results {
                if let Some(&(lor, p, _, var)) = go_map.get(&go_term) {
                    effect_sizes.push(lor);
                    variances.push(var);
                    p_values.push(p);
                    taxon_ids.push(*taxon_id);
                    
                    // Get taxonomic weight
                    let weight = taxonomic_weights.get(taxon_id).copied().unwrap_or(1.0);
                    bio_weights.push(weight);
                }
            }

            // Apply weighted combination based on number of species
            if !effect_sizes.is_empty() {
                let (tau_squared, mean_effect) = if effect_sizes.len() == 1 {
                    // If only one species, use its values directly
                    (0.0, effect_sizes[0])
                } else {
                    // Apply taxonomically weighted Paule-Mandel estimator
                    taxonomic_paule_mandel_estimator(
                        &effect_sizes, 
                        &variances, 
                        &bio_weights, 
                        tolerance, 
                        max_iterations
                    )
                };
                
                // Combine p-values with taxonomic weights
                let combined_p = combine_weighted_pvalues(&p_values, &bio_weights);
                
                // Store final result
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

        // Add group results to final estimates
        if !go_term_estimates.is_empty() {
            final_estimates.insert(group.clone(), go_term_estimates);
        }
    }

    final_estimates
}