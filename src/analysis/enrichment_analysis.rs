use rustc_hash::{FxHashMap, FxHashSet};
use fishers_exact::fishers_exact;
use crate::parsers::background_parser::{GOTermCount, GOTermID, TaxonID};
use statrs::distribution::{Hypergeometric, DiscreteCDF};
use rayon::prelude::*;

pub type ContingencyTable = [usize; 4];

#[derive(Debug, Clone)]
pub struct GOTermResults {
    pub log_odds_ratio: f64,
    pub p_value: f64,
    pub contingency_table: ContingencyTable,
}

#[derive(Debug, Clone, Copy)]
pub enum StatisticalTest {
    Fishers,
    Hypergeometric,
}

pub fn create_contingency_table(
    study_with_go: usize,
    background_with_go: usize,
    total_study: usize,
    total_background: usize,
) -> ContingencyTable {
    let a = study_with_go;
    let b = total_study.saturating_sub(study_with_go);
    let c = background_with_go.saturating_sub(study_with_go);
    let d = total_background
        .saturating_sub(background_with_go)
        .saturating_sub(total_study.saturating_sub(study_with_go));
    
    [a + 1, b + 1, c + 1, d + 1]
}

pub fn calculate_log_odds_ratio(counts: &ContingencyTable) -> f64 {
    let [a, b, c, d] = *counts;
    ((f64::from(a as u32) * f64::from(d as u32)) / 
     (f64::from(b as u32) * f64::from(c as u32))).ln()
}

pub fn calculate_p_value(counts: &ContingencyTable, test_type: StatisticalTest) -> f64 {
    match test_type {
        StatisticalTest::Fishers => fishers_test(counts),
        StatisticalTest::Hypergeometric => hypergeometric_test(counts),
    }
}

pub fn fishers_test(counts: &ContingencyTable) -> f64 {
    let counts_u32 = [
        counts[0] as u32,
        counts[1] as u32,
        counts[2] as u32,
        counts[3] as u32,
    ];
    
    match fishers_exact(&counts_u32) {
        Ok(result) => result.greater_pvalue,
        Err(_) => 1.0,
    }
}

pub fn hypergeometric_test(counts: &ContingencyTable) -> f64 {
    let k = counts[0] as u32;
    let n = (counts[0] + counts[1]) as u32;
    let K = (counts[0] + counts[2]) as u32;
    let N = (counts[0] + counts[1] + counts[2] + counts[3]) as u32;

    match Hypergeometric::new(N.into(), K.into(), n.into()) {
        Ok(dist) => dist.sf((k - 1).into()),
        Err(_) => 1.0,
    }
}

pub fn analyze_single_go_term(
    study_with_go: usize,
    background_with_go: usize,
    total_study: usize,
    total_background: usize,
    test_type: StatisticalTest,
) -> GOTermResults {
    let contingency_table = create_contingency_table(
        study_with_go,
        background_with_go,
        total_study,
        total_background,
    );
    
    GOTermResults {
        log_odds_ratio: calculate_log_odds_ratio(&contingency_table),
        p_value: calculate_p_value(&contingency_table, test_type),
        contingency_table,
    }
}

pub struct EnrichmentAnalysis {
    pub test_type: StatisticalTest,
}

impl EnrichmentAnalysis {
    pub fn new(
        test_type: StatisticalTest) -> Self {
        Self { test_type }
    }

    pub fn classic(
        &self,
        taxon_ids: &FxHashSet<TaxonID>,
        background_go_counts: &FxHashMap<TaxonID, GOTermCount>,
        study_go_counts: &FxHashMap<TaxonID, GOTermCount>,
        background_totals: &FxHashMap<TaxonID, usize>,
        study_totals: &FxHashMap<TaxonID, usize>,
    ) -> FxHashMap<TaxonID, FxHashMap<GOTermID, GOTermResults>> {
        let contingency_tables = self.create_contingency_tables(
            taxon_ids,
            background_go_counts,
            study_go_counts,
            background_totals,
            study_totals,
        );
        
        self.calculate_statistics(&contingency_tables)
    }

    fn create_contingency_tables(
        &self,
        taxon_ids: &FxHashSet<u32>,
        background_go_counts: &FxHashMap<TaxonID, GOTermCount>,
        study_go_counts: &FxHashMap<TaxonID, GOTermCount>,
        background_totals: &FxHashMap<TaxonID, usize>,
        study_totals: &FxHashMap<TaxonID, usize>,
    ) -> FxHashMap<TaxonID, FxHashMap<u32, ContingencyTable>> {
        taxon_ids.par_iter()
        .filter_map(|&taxon_id| {
            let tables = self.create_taxon_contingency_tables(
                taxon_id,
                background_go_counts,
                study_go_counts,
                background_totals,
                study_totals,
            );
            
            tables.map(|t| (taxon_id, t))
        })
        .collect() 
    }

    fn create_taxon_contingency_tables(
        &self,
        taxon_id: TaxonID,
        background_go_counts: &FxHashMap<TaxonID, GOTermCount>,
        study_go_counts: &FxHashMap<TaxonID, GOTermCount>,
        background_totals: &FxHashMap<TaxonID, usize>,
        study_totals: &FxHashMap<TaxonID, usize>,
    ) -> Option<FxHashMap<u32, ContingencyTable>> {
        let (background_counts, study_counts) = match (
            background_go_counts.get(&taxon_id),
            study_go_counts.get(&taxon_id)
        ) {
            (Some(bg), Some(study)) => (bg, study),
            _ => return None,
        };

        let total_background = *background_totals.get(&taxon_id).unwrap_or(&0);
        let total_study = *study_totals.get(&taxon_id).unwrap_or(&0);

        let tables: FxHashMap<u32, ContingencyTable> = study_counts
            .iter()
            .map(|(&go_id, &study_with_go)| {
                let background_with_go = *background_counts.get(&go_id).unwrap_or(&0);
                
                let table = create_contingency_table(
                    study_with_go,
                    background_with_go,
                    total_study,
                    total_background,
                );
                (go_id, table)
            })
            .collect();

        if tables.is_empty() {
            None
        } else {
            Some(tables)
        }
    }

    fn calculate_statistics(
        &self,
        go_counts: &FxHashMap<TaxonID, FxHashMap<u32, ContingencyTable>>,
    ) -> FxHashMap<TaxonID, FxHashMap<u32, GOTermResults>> {
        go_counts
            .par_iter()
            .map(|(&taxon_id, go_terms)| {
                let term_results = go_terms
                    .iter()
                    .map(|(&go_id, counts)| {
                        let stats = GOTermResults {
                            log_odds_ratio: calculate_log_odds_ratio(counts),
                            p_value: calculate_p_value(counts, self.test_type),
                            contingency_table: *counts,
                        };
                        (go_id, stats)
                    })
                    .collect();
                (taxon_id, term_results)
            })
            .collect()
    }
}