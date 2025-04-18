use polars::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use crate::parsers::background_parser::*;
use crate::analysis::enrichment_analysis::*;
use ndarray::{Array1, Array2, ArrayBase, OwnedRepr, Dim, Axis};
use ndarray_linalg::{Inverse, SVD};


pub fn read_vcv_matrix(
    matrix_path: &String,
) -> Result<DataFrame, PolarsError> {

    let mut vcv_matrix = CsvReadOptions::default()
        .try_into_reader_with_file_path(Some(matrix_path.into()))?
        .finish()?;

    vcv_matrix = vcv_matrix
        .lazy()
        .with_columns([
            col("taxa")
                .cast(DataType::String)
        ])
        .collect()?;
    
    Ok(vcv_matrix)
}

pub fn filter_vcv_matrix(
    vcv_matrix: DataFrame,
    taxon_ids: &FxHashSet<TaxonID>
) -> Result<DataFrame, PolarsError> {
    let mut taxon_ids_str: FxHashSet<String> = taxon_ids
        .iter()
        .map(|id| id.to_string())
        .collect();
    
    taxon_ids_str.insert("taxa".to_string());
    
    let taxon_ids_vec: Vec<String> = taxon_ids_str.iter().cloned().collect();
    let taxon_ids_series = Series::new(" ".into(), &taxon_ids_vec);

    let expression = col("taxa").is_in(lit(taxon_ids_series));
    let filtered_matrix = vcv_matrix.lazy()
        .filter(expression)
        .collect()?;
    
    let all_columns = filtered_matrix
        .get_column_names()
        .to_vec()
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<String>>();

    let all_columns_set = FxHashSet::from_iter(all_columns);
    let columns_to_drop = all_columns_set
        .difference(&taxon_ids_str)
        .map(|s| s.to_string())
        .collect::<Vec<String>>();

    let result_matrix = filtered_matrix.drop_many(&columns_to_drop);
    
    Ok(result_matrix)
}

pub fn svd_transform(
    vcv_matrix: &DataFrame
) -> ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>>{
    
    let array= vcv_matrix.to_ndarray::<Float64Type>(IndexOrder::C).unwrap();
    
    let svdcomp = array.svd(true, true).unwrap();
    let u = &svdcomp.0;
    let d = &svdcomp.1;

    let n = d.len();
    let mut sqrt_d_diag: ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>> = Array2::zeros((n, n));
    for i in 0..n {
        sqrt_d_diag[[i, i]] = d[i].sqrt();
    }
    
    let binding = u.clone().unwrap();
    let ut = &binding.t();
    let tnew = u.clone().unwrap().dot(&sqrt_d_diag).dot(ut);
    
    let mDnew: ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>> = tnew.inv().unwrap();

    mDnew

}

pub fn weighted_phylogenetic_regression(
    log_odds_ratios: &Array1<f64>,
    variances:&Array1<f64>,
    vcv_matrix: &DataFrame
) -> f64{
    let design_matrix: Array1<f64> = Array1::ones(log_odds_ratios.len());
    let weights: Array1<f64> = variances.mapv(|x| 1.0 / x);

    let w_size = weights.len();
    let mut w_matrix = Array2::zeros((w_size, w_size));
    for i in 0..w_size {
        w_matrix[[i, i]] = weights[i];
    }

    let mDnew = svd_transform(vcv_matrix);

    let log_odds_2d = log_odds_ratios.view().insert_axis(Axis(1));
    let design_matrix_2d = design_matrix.view().insert_axis(Axis(1));
    
    let y_new = mDnew.dot(&log_odds_2d);
    let x_new = mDnew.dot(&design_matrix_2d);

    let xt_w_x = x_new.t().dot(&w_matrix).dot(&x_new);
    let xt_w_e = x_new.t().dot(&w_matrix).dot(&y_new);

    let xt_w_x_inv = 1.0 / xt_w_x[[0, 0]];
    
    let b_pma = xt_w_x_inv * xt_w_e[[0, 0]];

    b_pma

}

pub fn phylogenetic_meta_analysis(
    taxon_ids: &FxHashSet<TaxonID>,
    lineage_results: FxHashMap<String, FxHashMap<TaxonID, FxHashMap<GOTermID, GOTermResults>>>,
    matrix_path: &String
) -> FxHashMap<String, FxHashMap<GOTermID, f64>> {
    
    // Read and filter the main VCV matrix
    let mut vcv_matrix = read_vcv_matrix(matrix_path).unwrap();
    vcv_matrix = filter_vcv_matrix(vcv_matrix, taxon_ids).unwrap();

    // Get the taxon order from the filtered matrix
    let taxon_order: Vec<TaxonID> = vcv_matrix
        .column("taxa")
        .unwrap()
        .str()
        .unwrap()
        .into_iter()
        .filter_map(|opt_s| opt_s.and_then(|s| s.parse::<u32>().ok()))
        .collect();

    let mut results = FxHashMap::default();

    // Process each taxonomic level
    for (level, taxon_map) in lineage_results {
        let mut level_results = FxHashMap::default();
        
        // Get all GO terms for this level
        let all_go_terms: FxHashSet<GOTermID> = taxon_map
            .values()
            .flat_map(|go_term_map| go_term_map.keys().cloned())
            .collect();
        
        // Process each GO term
        for go_term in all_go_terms {
            let mut relevant_taxon_ids = FxHashSet::default();
            let mut log_odds_ratios = Vec::new();
            let mut variances = Vec::new();
            
            // Collect data for taxa that have this GO term
            for &taxon_id in &taxon_order {
                if let Some(go_term_map) = taxon_map.get(&taxon_id) {
                    if let Some(result) = go_term_map.get(&go_term) {
                        log_odds_ratios.push(result.log_odds_ratio);
                        variances.push(result.variance);
                        relevant_taxon_ids.insert(taxon_id);
                    }
                }
            }
            
            // Skip if we don't have any data
            if relevant_taxon_ids.is_empty() {
                continue;
            }
            
            let b_pma: f64;
            
            if log_odds_ratios.len() == 1 {
                // If we have just one data point, use it directly
                b_pma = log_odds_ratios[0];
            } else {
                // For multiple data points, perform phylogenetic regression
                let go_term_vcv_matrix = filter_vcv_matrix(vcv_matrix.clone(), &relevant_taxon_ids).unwrap();
                
                let log_odds_array = Array1::from(log_odds_ratios);
                let variance_array = Array1::from(variances);

                b_pma = weighted_phylogenetic_regression(
                    &log_odds_array,
                    &variance_array,
                    &go_term_vcv_matrix
                );
            }
            
            // Store the b_pma value
            level_results.insert(go_term, b_pma);
        }
        
        // Add results for this level if we have any
        if !level_results.is_empty() {
            results.insert(level, level_results);
        }
    }
    
    println!("Phylogenetic meta-analysis results: {:?}", results); 
    results
}

// FxHashMap<Taxonomic level, FxHashMap<Taxon ID, FxHashMap<GO Term ID, (log(Odds Ratio), P-value, Contingency Table, Variance)>>>

// Process taxonomic levels in parallel 
/*
For each taxonomic level, we will:
    1. We will process each GO term
    For each term we need to: 
        1.1 Find the taxon ids that have this term 
        1.2 Filter the VCV matrix to only include these taxon ids
        1.3 We need create an Array with log(odds ratio)        DONE
        1.4 We need create an Array with weights (1/variance)       DONE
        1.5 We need create an Array with ones           DONE
        IMPORTANT: All these arrays must have the values in the same order as the taxa appear in the VCV matrix
        1.6 Calculate the SVD of the VCV matrix             DONE
        1.7 Perform the calculations with the SVD       DONE
        1.8 Fit the model       DONE
        1.9 Combine p-values
        1.10 Store the results in the FxHashMap
*/