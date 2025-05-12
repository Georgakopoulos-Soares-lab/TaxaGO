use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Result, Error, ErrorKind};
use std::path::Path;
use std::sync::Arc;
use rayon::prelude::*;
use crate::parsers::study_parser::*;
use compact_str::CompactString;

pub type TaxonID = u32;
pub type GOTermID = u32;
pub type Protein = Arc<CompactString>;

pub type ProteinCount = FxHashMap<TaxonID, usize>;
pub type ProteinToGO = FxHashMap<CompactString, FxHashSet<GOTermID>>;
pub type GOTermCount = FxHashMap<GOTermID, usize>;
pub type GOTermToProteinSet = FxHashMap<GOTermID, FxHashSet<Protein>>;

#[derive(Debug, Default, Clone)]
pub struct BackgroundPop {
    pub taxon_protein_count: ProteinCount,
    pub protein_to_go: FxHashMap<TaxonID, ProteinToGO>,
    pub go_term_count: FxHashMap<TaxonID, GOTermCount>,
    pub go_term_to_protein_set: FxHashMap<TaxonID, GOTermToProteinSet>
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EvidenceCategory {
    Experimental,
    Phylogenetic,
    Computational,
    Author,
    Curator,
    Electronic
}

pub fn map_code_to_category(
    code: &CompactString
) -> Result<EvidenceCategory> {
    match code.as_str() {
        "EXP" | "IDA" | "IPI" | "IMP" | "IGI" | "IEP" | "HTP" | "HDA" | "HMP" | "HGI" | "HEP" => Ok(EvidenceCategory::Experimental),
        "IBA" | "IBD" | "IKR" | "IRD" => Ok(EvidenceCategory::Phylogenetic),
        "ISS" | "ISO" | "ISA" | "ISM" | "IGC" | "RCA" => Ok(EvidenceCategory::Computational),
        "TAS" | "NAS" => Ok(EvidenceCategory::Author),
        "IC" | "ND" => Ok(EvidenceCategory::Curator),
        "IEA" => Ok(EvidenceCategory::Electronic),
        _ => Err(Error::new(
            ErrorKind::InvalidData,
            format!("Unrecognized evidence code: {}", code)
        ))
    }
}

pub fn map_input_to_category(
    cli_input: String
) -> Result<Vec<EvidenceCategory>> {
    let parts: Vec<String> = cli_input
        .to_lowercase()
        .split(',')
        .map(|part| part.trim().to_string())
        .filter(|part| !part.is_empty())
        .collect();

    let contains_all = parts.iter().any(|part| part == "all");

    if contains_all {
        Ok(vec![
            EvidenceCategory::Experimental,
            EvidenceCategory::Phylogenetic,
            EvidenceCategory::Computational,
            EvidenceCategory::Author,
            EvidenceCategory::Curator,
            EvidenceCategory::Electronic
        ])
    } else {
        let mut categories = Vec::new();
        
        for trimmed_part in parts {
            let category = match trimmed_part.as_str() {
                "experimental" => EvidenceCategory::Experimental,
                "phylogenetic" => EvidenceCategory::Phylogenetic,
                "computational" => EvidenceCategory::Computational,
                "author" => EvidenceCategory::Author,
                "curator" => EvidenceCategory::Curator,
                "electronic" => EvidenceCategory::Electronic,
                _ => return Err(Error::new(
                    ErrorKind::InvalidInput,
                    format!("Unrecognized evidence category: {}. Valid categories are: experimental, phylogenetic, computational, author, curator, automatic, or all.", trimmed_part)
                ))
            };
            
            categories.push(category);
        }
        
        if categories.is_empty() {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "No valid evidence categories provided. Valid categories are: experimental, phylogenetic, computational, author, curator, automatic, or all."
            ));
        }
        
        Ok(categories)
    }
}

impl BackgroundPop {
    pub fn read_background_pop(
        taxon_ids: &FxHashSet<TaxonID>, 
        dir: &str,
        categories: &Vec<EvidenceCategory>
    ) -> Result<Option<Self>> {
        
        let (taxon_protein_count, protein_to_go, go_term_count, go_term_to_protein_set) = taxon_ids
            .par_iter()
            .map(|&taxon_id| {
                let taxon_background_file_path = format!("{}/{}_background.txt", dir, taxon_id);
                
                match process_single_taxon(
                    &taxon_background_file_path,
                    categories
                ) {
                    Ok(Some(data)) => (taxon_id, Some(data)),
                    Ok(None) => {
                        eprintln!("No data found for taxon {}", taxon_id);
                        (taxon_id, None)
                    },
                    Err(err) => {
                        eprintln!("Error processing taxon {}: {}", taxon_id, err);
                        (taxon_id, None)
                    }
                }
            })
            .fold(
                || (FxHashMap::default(), FxHashMap::default(), FxHashMap::default(), FxHashMap::default()),
                |(mut pc, mut pg, mut gc, mut gp), (taxon_id, result)| {
                    if let Some((protein_count, protein_to_go_map, go_term_counts, go_term_to_protein)) = result {
                        pc.insert(taxon_id, protein_count);
                        pg.insert(taxon_id, protein_to_go_map);
                        gc.insert(taxon_id, go_term_counts);
                        gp.insert(taxon_id, go_term_to_protein);
                    }
                    (pc, pg, gc, gp)
                }
            )
            .reduce(
                || (FxHashMap::default(), FxHashMap::default(), FxHashMap::default(), FxHashMap::default()),
                |(mut pc1, mut pg1, mut gc1, mut gp1), (pc2, pg2, gc2, gp2)| {
                    pc1.extend(pc2);
                    pg1.extend(pg2);
                    gc1.extend(gc2);
                    gp1.extend(gp2);
                    (pc1, pg1, gc1, gp1)
                }
            );
    
        Ok(Some(Self {
            taxon_protein_count,
            protein_to_go,
            go_term_count,
            go_term_to_protein_set
        }))
    }

    pub fn filter_by_study_population(&mut self, taxon_ids: &FxHashSet<TaxonID>, study_pop: &StudyPop) {
        for &taxon_id in taxon_ids {
            let study_terms = match study_pop.go_term_count.get(&taxon_id) {
                Some(terms) => terms,
                None => continue,
            };
            
            let terms_to_remove_set: FxHashSet<GOTermID> = match self.go_term_count.get(&taxon_id) {
                Some(bg_terms) => bg_terms
                    .keys()
                    .filter(|&term_id| !study_terms.contains_key(term_id))
                    .copied()
                    .collect(),
                None => continue,
            };
            
            if terms_to_remove_set.is_empty() {
                continue;
            }
            
            if let Some(count_map) = self.go_term_count.get_mut(&taxon_id) {
                for &term_id in &terms_to_remove_set {
                    count_map.remove(&term_id);
                }
            }
            
            if let Some(term_map) = self.go_term_to_protein_set.get_mut(&taxon_id) {
                for &term_id in &terms_to_remove_set {
                    term_map.remove(&term_id);
                }
            }
            
            if let Some(protein_go_map) = self.protein_to_go.get_mut(&taxon_id) {
                protein_go_map.par_iter_mut().for_each(|(_, go_terms)| {
                    go_terms.retain(|go_term| !terms_to_remove_set.contains(go_term));
                });
            }
        }
    }   
}

fn process_single_taxon(
    taxon_background_path: impl AsRef<Path>,
    categories: &Vec<EvidenceCategory>
) -> Result<Option<(usize, ProteinToGO, GOTermCount, GOTermToProteinSet)>> {
    
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

    let mut protein_to_go_map: FxHashMap<CompactString, FxHashSet<GOTermID>> = FxHashMap::default();
    let mut go_term_counts: FxHashMap<GOTermID, usize> = FxHashMap::default();
    let mut go_term_to_protein_set: FxHashMap<GOTermID, FxHashSet<Protein>> = FxHashMap::default();

    for line_result in reader.lines() {
        let line = match line_result {
            Ok(l) => l,
            Err(e) => return Err(Error::new(
                ErrorKind::Other,
                format!("Error reading line: {}", e)
            ))
        };
        let parts: Vec<&str> = line.split('\t').collect();

        if parts.len() < 3 {
            continue;
        }

        let code = CompactString::new(parts[2]);
        let category = map_code_to_category(&code)?;

        if categories.contains(&category){
            let protein_str = CompactString::new(parts[0]);
            let protein_id = Arc::new(protein_str.clone());
            if let Some(go_str) = parts[1].strip_prefix("GO:") {
                if let Ok(go_id) = go_str.parse::<u32>() {
                    protein_to_go_map
                        .entry(protein_str)
                        .or_insert_with(FxHashSet::default)
                        .insert(go_id);
    
                    let is_new_association = go_term_to_protein_set
                        .entry(go_id)
                        .or_insert_with(FxHashSet::default)
                        .insert(Arc::clone(&protein_id));
                    
                    if is_new_association {
                        *go_term_counts.entry(go_id).or_insert(0) += 1;
                    }
    
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

