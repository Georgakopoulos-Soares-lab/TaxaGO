use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Result, Error, ErrorKind};
use std::path::Path;
use rayon::prelude::*;
use crate::parsers::study_parser::*;

pub type TaxonID = u32;
pub type GOTermID = u32;
pub type Protein = String;

pub type ProteinCount = HashMap<TaxonID, usize>;
pub type ProteinToGO = HashMap<Protein, HashSet<GOTermID>>;
pub type GOTermCount = HashMap<GOTermID, usize>;
pub type GOTermToProteinSet = HashMap<GOTermID, HashSet<Protein>>;

#[derive(Debug, Default, Clone)]
pub struct BackgroundPop {
    pub taxon_protein_count: ProteinCount,
    pub protein_to_go: HashMap<TaxonID, ProteinToGO>,
    pub go_term_count: HashMap<TaxonID, GOTermCount>,
    pub go_term_to_protein_set: HashMap<TaxonID, GOTermToProteinSet>
}

impl BackgroundPop {

    pub fn read_background_pop(
        taxon_ids: &HashSet<TaxonID>, 
        dir: &str) -> Result<Option<Self>> {
        
            let results: Vec<(TaxonID, Option<(usize, ProteinToGO, GOTermCount, GOTermToProteinSet)>)> = taxon_ids
            .par_iter()
            .map(|&taxon_id| {
                let taxon_background_file_path = format!("{}/{}_background.txt", dir, taxon_id);
                
                let result = match process_single_taxon(&taxon_background_file_path) {
                    Ok(Some(data)) => Some(data),
                    Ok(None) => None,
                    Err(err) => {
                        eprintln!("Error processing taxon {}: {}", taxon_id, err);
                        None
                    }
                };
                
                (taxon_id, result)
            })
            .collect();

        let mut taxon_protein_count: ProteinCount = HashMap::new();
        let mut protein_to_go: HashMap<TaxonID, ProteinToGO> = HashMap::new();
        let mut go_term_count: HashMap<TaxonID, GOTermCount> = HashMap::new();
        let mut go_term_to_protein_set: HashMap<TaxonID, GOTermToProteinSet> = HashMap::new();

        for (taxon_id, result) in results {
            if let Some((
                protein_count, 
                protein_to_go_map, 
                go_term_counts,
                go_term_to_protein
            )) = result {
                taxon_protein_count.insert(taxon_id, protein_count);
                protein_to_go.insert(taxon_id, protein_to_go_map);
                go_term_count.insert(taxon_id, go_term_counts);
                go_term_to_protein_set.insert(taxon_id, go_term_to_protein);
            } else {
                eprintln!("No valid data found for taxon {}", taxon_id);
            }
        }

        Ok(Some(Self {
            taxon_protein_count,
            protein_to_go,
            go_term_count,
            go_term_to_protein_set
        }))
    }

    pub fn filter_by_study_population(
        &mut self,
        taxon_ids: &HashSet<TaxonID>,
        study_pop: &StudyPop) {
        
        for taxon_id in taxon_ids {

            let study_go_terms = match study_pop.go_term_count.get(&taxon_id) {
                Some(terms) => terms.keys().collect::<HashSet<&GOTermID>>(),
                None => continue,
            };
            
            let terms_to_remove = if let Some(bg_terms) = self.go_term_count.get(&taxon_id) {
                bg_terms
                    .keys()
                    .filter(|term_id| !study_go_terms.contains(term_id))
                    .cloned()
                    .collect::<Vec<GOTermID>>()
            } else {
                continue;
            };
            
            if terms_to_remove.is_empty() {
                continue;
            }
            
            if let Some(count_map) = self.go_term_count.get_mut(&taxon_id) {
                for term_id in &terms_to_remove {
                    count_map.remove(term_id);
                }
            }
            
            if let Some(term_map) = self.go_term_to_protein_set.get_mut(&taxon_id) {
                for term_id in &terms_to_remove {
                    term_map.remove(term_id);
                }
            }
            
            if let Some(protein_go_map) = self.protein_to_go.get_mut(&taxon_id) {
                for protein_go_terms in protein_go_map.values_mut() {
                    protein_go_terms.retain(|go_term| !terms_to_remove.contains(go_term));
                }
            }
        }
    }   
}

fn process_single_taxon(
    taxon_background_path: impl AsRef<Path>) -> Result<Option<(usize, ProteinToGO, GOTermCount, GOTermToProteinSet)>> {
    
    if !taxon_background_path.as_ref().is_file() {
        return Ok(None);
    }

    let file = match File::open(taxon_background_path.as_ref()) {
        Ok(f) => f,
        Err(e) => return Err(Error::new(
            ErrorKind::Other,
            format!("Failed to open file: {}", e)
        ))
    };

    let reader = BufReader::with_capacity(128 * 1024, file);

    let mut protein_to_go_map: HashMap<Protein, HashSet<GOTermID>> = HashMap::new();
    let mut go_term_counts: HashMap<GOTermID, usize> = HashMap::new();
    let mut go_term_to_protein_set: HashMap<GOTermID, HashSet<Protein>> = HashMap::new();

    for line_result in reader.lines() {
        let line = match line_result {
            Ok(l) => l,
            Err(e) => return Err(Error::new(
                ErrorKind::Other,
                format!("Error reading line: {}", e)
            ))
        };
        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() < 2 {
            continue;
        }

        let protein_id = parts[0].to_string();
        if let Some(go_str) = parts[1].strip_prefix("GO:") {
            if let Ok(go_id) = go_str.parse::<u32>() {
                protein_to_go_map
                    .entry(protein_id.clone())  
                    .or_insert_with(HashSet::new)
                    .insert(go_id);

                let is_new_association = go_term_to_protein_set
                    .entry(go_id)
                    .or_insert_with(HashSet::new)
                    .insert(protein_id);
                
                if is_new_association {
                    *go_term_counts.entry(go_id).or_insert(0) += 1;
                }

            }
        }
    }

    Ok(Some((
        protein_to_go_map.len(), 
        protein_to_go_map, 
        go_term_counts,
        go_term_to_protein_set
    )))
}

