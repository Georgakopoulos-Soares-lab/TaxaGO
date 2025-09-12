use polars::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use crate::parsers::background_parser::*;
use crate::analysis::enrichment_analysis::*;
use ndarray::{Array1, Array2, ArrayBase, OwnedRepr, Dim, Axis};
use nalgebra::DMatrix;
use rand::{
    seq::SliceRandom,
    SeedableRng,
    rngs::StdRng
};
use rayon::prelude::*;
use std::path::PathBuf;
#[derive(Debug, Clone)]
pub struct TaxonomyGOResult {
    pub log_odds_ratio: f64,
    pub p_value: f64,
    pub species_number: usize
}
fn ndarray2_to_nalgebra(arr: &Array2<f64>) -> DMatrix<f64> {
    let (nrows, ncols) = arr.dim();
    DMatrix::from_row_slice(nrows, ncols, arr.as_slice().expect("Input ndarray was not contiguous"))
}
fn nalgebra_to_ndarray2(mat: &DMatrix<f64>) -> Array2<f64> {
    let (nrows, ncols) = mat.shape();
    Array2::from_shape_fn((nrows, ncols), |(r, c)| mat[(r, c)])
}

pub fn read_vcv_matrix(
    matrix_path: PathBuf,
) -> Result<DataFrame, PolarsError> {
    let mut vcv_matrix = CsvReadOptions::default()
        .with_has_header(true)
        .with_infer_schema_length(Some(15000))
        .try_into_reader_with_file_path(Some(matrix_path))?
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
    vcv_matrix_df: &DataFrame
) -> ArrayBase<OwnedRepr<f64>, Dim<[usize; 2]>> { 
    
    let nd_array_original = vcv_matrix_df.to_ndarray::<Float64Type>(IndexOrder::C).unwrap();
    let na_matrix = ndarray2_to_nalgebra(&nd_array_original);
    
    let svd_result = na_matrix.svd(true, true);
    let u_nalgebra = svd_result.u.expect("SVD U matrix not computed by nalgebra");
    let s_nalgebra_vec = svd_result.singular_values;

    let dim = s_nalgebra_vec.len();
    
    let eps = 1e-6;
    let inv_sqrt = DMatrix::from_fn(dim, dim, |r,c| if r==c { 
        if s_nalgebra_vec[r] > eps {
            1.0 / s_nalgebra_vec[r].sqrt()
        } else {
            0.0
        }
    } else { 
        0.0 
    });
    let mDnew_nalgebra = &u_nalgebra * inv_sqrt * u_nalgebra.transpose();

    nalgebra_to_ndarray2(&mDnew_nalgebra)
}

pub fn phylogenetic_meta_analysis_calculation(
    log_odds_array: &Array1<f64>,
    variance_array: &Array1<f64>,
    vcv_matrix: &DataFrame,
    permutations: u32,
    seed: u64
) -> (f64, f64) {
    let weights = variance_array.mapv(|x| 1.0 / x);
    
    let n = log_odds_array.len();
    
    let design_matrix: Array1<f64> = Array1::ones(n);
    let design_matrix_2d = design_matrix.view().insert_axis(Axis(1));
    
    let mDnew = svd_transform(vcv_matrix);
    
    let x_new = mDnew.dot(&design_matrix_2d);
    
    let log_odds_2d = log_odds_array.view().insert_axis(Axis(1));
        
    let mut w_matrix = Array2::zeros((n, n));
    for i in 0..n {
        w_matrix[[i, i]] = weights[i];
    }
    
    let y_new = mDnew.dot(&log_odds_2d);
    
    let xt_w_x = x_new.t().dot(&w_matrix).dot(&x_new);
    let xt_w_e = x_new.t().dot(&w_matrix).dot(&y_new);
    let xt_w_x_inv = 1.0 / xt_w_x[[0, 0]];
    
    let b_pma = xt_w_x_inv * xt_w_e[[0, 0]];
    
    let indices: Vec<usize> = (0..n).collect();
    
    let exceeds_count: u32 = (0..permutations)
        .into_par_iter()
        .map(|perm_idx| {
            let thread_seed = seed.wrapping_add(perm_idx as u64);
            let mut thread_rng = StdRng::seed_from_u64(thread_seed);
            
            let mut indices_rand = indices.clone();
            indices_rand.shuffle(&mut thread_rng);
            
            let mut rand_data_temp: Vec<(usize, f64, f64)> = Vec::with_capacity(n);
            for i in 0..n {
                rand_data_temp.push((indices_rand[i], log_odds_array[i], weights[i]));
            }
            
            rand_data_temp.sort_by_key(|&(id, _, _)| id);
            
            let log_odds_rand: Vec<f64> = rand_data_temp.iter().map(|&(_, lor, _)| lor).collect();
            let weights_rand: Vec<f64> = rand_data_temp.iter().map(|&(_, _, w)| w).collect();
            
            let mut w_matrix_rand = Array2::zeros((n, n));
            for i in 0..n {
                w_matrix_rand[[i, i]] = weights_rand[i];
            }
            
            let log_odds_rand_array = Array1::from(log_odds_rand);
            let log_odds_rand_2d = log_odds_rand_array.view().insert_axis(Axis(1));
            
            let y_new_rand = mDnew.dot(&log_odds_rand_2d);
            
            let xt_w_x_rand = x_new.t().dot(&w_matrix_rand).dot(&x_new);
            let xt_w_e_rand = x_new.t().dot(&w_matrix_rand).dot(&y_new_rand);
            let xt_w_x_inv_rand = 1.0 / xt_w_x_rand[[0, 0]];
            let bpma_rand = xt_w_x_inv_rand * xt_w_e_rand[[0, 0]];
            
            if bpma_rand >= b_pma { 1 } else { 0 }
        })
        .sum();

    let p_bpma = exceeds_count + 1;
    let p_value = p_bpma as f64 / (permutations as f64 + 1.0);
    
    (b_pma, p_value)
}


pub fn phylogenetic_meta_analysis(
    taxon_ids: &FxHashSet<TaxonID>,
    lineage_results: FxHashMap<String, FxHashMap<TaxonID, FxHashMap<GOTermID, GOTermResults>>>,
    superkingdom_vcv_matrix: DataFrame,
    permutations: u32

) -> FxHashMap<String, FxHashMap<GOTermID, TaxonomyGOResult>> {
    
    let vcv_matrix = filter_vcv_matrix(superkingdom_vcv_matrix, taxon_ids).unwrap();

    let taxon_order: Vec<TaxonID> = vcv_matrix
        .column("taxa")
        .unwrap()
        .str()
        .unwrap()
        .into_iter()
        .filter_map(|opt_s| opt_s.and_then(|s| s.parse::<u32>().ok()))
        .collect();

    let mut results = FxHashMap::default();

    for (level, taxon_map) in lineage_results {
        let mut level_results = FxHashMap::default();
        
        let all_go_terms: FxHashSet<GOTermID> = taxon_map
            .values()
            .flat_map(|go_term_map| go_term_map.keys().cloned())
            .collect();
        
        for go_term in all_go_terms {
            let mut relevant_taxon_ids = FxHashSet::default();
            let mut log_odds_ratios = Vec::new();
            let mut p_values = Vec::new();
            let mut variances = Vec::new();
            
            for &taxon_id in &taxon_order {
                if let Some(go_term_map) = taxon_map.get(&taxon_id) {
                    if let Some(result) = go_term_map.get(&go_term) {
                        log_odds_ratios.push(result.log_odds_ratio);
                        p_values.push(result.p_value);
                        variances.push(result.variance);
                        relevant_taxon_ids.insert(taxon_id);
                    }
                }
            }
            
            if relevant_taxon_ids.is_empty() {
                continue;
            }
            
            let b_pma: f64;
            let p_value: f64;
            
            if log_odds_ratios.len() == 1 {
                b_pma = log_odds_ratios[0];
                p_value = p_values[0];
            } else {
                let go_term_vcv_matrix = filter_vcv_matrix(vcv_matrix.clone(), &relevant_taxon_ids).unwrap();
                
                let log_odds_array = Array1::from(log_odds_ratios);
                let variance_array = Array1::from(variances);
                
                let (b_pma_result, p_value_result) = phylogenetic_meta_analysis_calculation(
                    &log_odds_array,
                    &variance_array,
                    &go_term_vcv_matrix,
                    permutations,
                    42
                );
                
                b_pma = b_pma_result;
                p_value = p_value_result;
            }
            let num_species_with_go_term = relevant_taxon_ids.len() as u32;
            let go_result = TaxonomyGOResult {
                log_odds_ratio: b_pma,
                p_value: p_value,
                species_number: (num_species_with_go_term as usize),
            };

            level_results.insert(go_term, go_result);

        }
        
        if !level_results.is_empty() {
            results.insert(level, level_results);
        }
    }
    
    results
}