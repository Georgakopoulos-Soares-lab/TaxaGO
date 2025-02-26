use std::collections::{HashMap, HashSet};
use fishers_exact::fishers_exact;
use crate::parsers::background_parser::{GOTermCount, TaxonID};
use statrs::distribution::{Hypergeometric, DiscreteCDF};

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

pub struct EnrichmentAnalysis {
    min_count: usize,
    test_type: StatisticalTest,
}

impl EnrichmentAnalysis {
    pub fn new(min_count: usize, test_type: StatisticalTest) -> Self {
        Self { min_count, test_type }
    }

    pub fn analyze(
        &self,
        taxon_ids: &HashSet<u32>,
        background_go_counts: &HashMap<TaxonID, GOTermCount>,
        study_go_counts: &HashMap<TaxonID, GOTermCount>,
        background_totals: &HashMap<TaxonID, usize>,
        study_totals: &HashMap<TaxonID, usize>,
    ) -> HashMap<TaxonID, HashMap<u32, GOTermResults>> {
        let contingency_tables = self.create_contingency_tables(
            taxon_ids,
            background_go_counts,
            study_go_counts,
            background_totals,
            study_totals,
        );
        
        let filtered_tables = self.filter_go_terms(&contingency_tables);
        self.calculate_statistics(&filtered_tables)
    }

    fn create_contingency_tables(
        &self,
        taxon_ids: &HashSet<u32>,
        background_go_counts: &HashMap<TaxonID, GOTermCount>,
        study_go_counts: &HashMap<TaxonID, GOTermCount>,
        background_totals: &HashMap<TaxonID, usize>,
        study_totals: &HashMap<TaxonID, usize>,
    ) -> HashMap<TaxonID, HashMap<u32, ContingencyTable>> {
        let mut result = HashMap::with_capacity(taxon_ids.len());

        for &taxon_id in taxon_ids {
            if let Some(tables) = self.create_taxon_contingency_tables(
                taxon_id,
                background_go_counts,
                study_go_counts,
                background_totals,
                study_totals,
            ) {
                result.insert(taxon_id, tables);
            }
        }

        result
    }

    fn create_taxon_contingency_tables(
        &self,
        taxon_id: TaxonID,
        background_go_counts: &HashMap<TaxonID, GOTermCount>,
        study_go_counts: &HashMap<TaxonID, GOTermCount>,
        background_totals: &HashMap<TaxonID, usize>,
        study_totals: &HashMap<TaxonID, usize>,
    ) -> Option<HashMap<u32, ContingencyTable>> {
        let (background_counts, study_counts) = match (
            background_go_counts.get(&taxon_id),
            study_go_counts.get(&taxon_id)
        ) {
            (Some(bg), Some(study)) => (bg, study),
            _ => return None,
        };

        let total_background = *background_totals.get(&taxon_id).unwrap_or(&0);
        let total_study = *study_totals.get(&taxon_id).unwrap_or(&0);

        let tables: HashMap<u32, ContingencyTable> = study_counts
            .iter()
            .map(|(&go_id, &study_with_go)| {
                let background_with_go = *background_counts.get(&go_id).unwrap_or(&0);
                let study_with_go = study_with_go;
                let table = self.calculate_contingency_table(
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

    fn calculate_contingency_table(
        &self,
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

    fn filter_go_terms(
        &self,
        go_counts: &HashMap<TaxonID, HashMap<u32, ContingencyTable>>,
    ) -> HashMap<TaxonID, HashMap<u32, ContingencyTable>> {
        go_counts
            .iter()
            .filter_map(|(&taxon_id, go_terms)| {
                let filtered_terms: HashMap<u32, ContingencyTable> = go_terms
                    .iter()
                    .filter(|(_, counts)| counts[0] >= self.min_count)
                    .map(|(&go_id, counts)| (go_id, *counts))
                    .collect();

                if filtered_terms.is_empty() {
                    None
                } else {
                    Some((taxon_id, filtered_terms))
                }
            })
            .collect()
    }

    fn calculate_statistics(
        &self,
        go_counts: &HashMap<TaxonID, HashMap<u32, ContingencyTable>>,
    ) -> HashMap<TaxonID, HashMap<u32, GOTermResults>> {
        go_counts
            .iter()
            .map(|(&taxon_id, go_terms)| {
                let term_results = go_terms
                    .iter()
                    .map(|(&go_id, counts)| {
                        let stats = GOTermResults {
                            log_odds_ratio: self.calculate_log_odds_ratio(counts),
                            p_value: self.calculate_p_value(counts),
                            contingency_table: *counts,
                        };
                        (go_id, stats)
                    })
                    .collect();
                (taxon_id, term_results)
            })
            .collect()
    }

    fn calculate_log_odds_ratio(&self, counts: &ContingencyTable) -> f64 {
        let [a, b, c, d] = *counts;
        ((f64::from(a as u32) * f64::from(d as u32)) / 
         (f64::from(b as u32) * f64::from(c as u32))).ln()
    }

    fn calculate_p_value(&self, counts: &[usize; 4]) -> f64 {
        match self.test_type {
            StatisticalTest::Fishers => self.fishers_test(counts),
            StatisticalTest::Hypergeometric => self.hypergeometric_test(counts),
        }
    }

    fn fishers_test(&self, counts: &[usize; 4]) -> f64 {
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

    fn hypergeometric_test(&self, counts: &[usize; 4]) -> f64 {
        let k = counts[0] as u32; // Genes in study list with GO term
        let n = (counts[0] + counts[1]) as u32; // Total genes in study
        let K = (counts[0] + counts[2]) as u32; // Total genes with GO term
        let N = (counts[0] + counts[1] + counts[2] + counts[3]) as u32; // Total genes

        match Hypergeometric::new(N.into(), K.into(), n.into()) {
            Ok(dist) => dist.sf((k - 1).into()),
            Err(_) => 1.0,
        }
    }
}
