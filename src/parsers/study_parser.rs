use std::collections::{HashMap, HashSet};
use std::fs::{read_dir, File};
use std::io::{BufRead, BufReader, Error, ErrorKind, Result};
use std::path::Path;
use csv::Reader;  
use rayon::prelude::*;
use crate::parsers::background_parser::*;

#[derive(Debug, Default, Clone)]
pub struct StudyPop {
    pub taxon_map: HashMap<TaxonID, HashSet<Protein>>,
    pub taxon_protein_count: HashMap<TaxonID, usize>,
    pub go_term_count: HashMap<TaxonID, GOTermCount>,
    pub go_term_to_protein_set: HashMap<TaxonID, HashMap<GOTermID, HashSet<Protein>>>
}

impl StudyPop {

    pub fn from_csv_file(
        csv_file: impl AsRef<Path>,
        protein_to_go: &HashMap<TaxonID, ProteinToGO>,
    ) -> Result<Option<Self>> {
        
        let mut taxon_map: HashMap<TaxonID, HashSet<Protein>> = HashMap::new();
        
        let file = File::open(csv_file)?;
        let mut csv_reader = Reader::from_reader(file);
    
        let taxon_ids: Vec<TaxonID> = csv_reader
            .headers()?
            .clone()
            .iter()
            .filter_map(|id| id.parse::<u32>().ok())
            .collect();
            
        for result in csv_reader.records() {
            let record = result?;
            for (index, protein) in record.iter().enumerate() {
                if let Some(&taxon_id) = taxon_ids.get(index) {
                    if !protein.is_empty() {
                        taxon_map
                            .entry(taxon_id)
                            .or_insert_with(HashSet::new)
                            .insert(protein.to_string());
                    }
                }
            }
        }
        
        let mut taxon_protein_count = HashMap::with_capacity(taxon_map.len());
        let mut go_term_to_protein_set = HashMap::with_capacity(taxon_map.len());
        let mut go_term_count = HashMap::with_capacity(taxon_map.len());
        
        for (&taxon_id, proteins) in &taxon_map {
            taxon_protein_count.insert(taxon_id, proteins.len());
            
            if let Some(taxon_protein_to_go) = protein_to_go.get(&taxon_id) {
                let mut go_term_map = HashMap::new();
                
                for protein in proteins {
                    if let Some(go_terms) = taxon_protein_to_go.get(protein) {
                        for &go_term_id in go_terms {
                            go_term_map
                                .entry(go_term_id)
                                .or_insert_with(HashSet::new)
                                .insert(protein.to_string());
                        }
                    }
                }
                
                if !go_term_map.is_empty() {
                    let mut term_count = HashMap::new();
                    for (&go_term_id, go_proteins) in &go_term_map {
                        term_count.insert(go_term_id, go_proteins.len());
                    }
                    
                    go_term_count.insert(taxon_id, term_count);
                    go_term_to_protein_set.insert(taxon_id, go_term_map);
                }
            }
        }
        
        Ok(Some(Self {
            taxon_map,
            taxon_protein_count,
            go_term_count,
            go_term_to_protein_set
        }))
    }

    pub fn read_study_pop(
        study_data: impl AsRef<Path>,
        protein_to_go: &HashMap<TaxonID, ProteinToGO>
    ) -> Result<Option<Self>> {
        
        let path = study_data.as_ref();
        
        match path.extension().and_then(|s| s.to_str()) {
            Some("csv") if path.is_file() => {
                return StudyPop::from_csv_file(path, protein_to_go);
            }
            Some("fa") | Some("fasta") if path.is_file() => {
                if let Ok(Some((taxon_id, protein_set, go_term_count_map, go_term_to_proteins))) = parse_fasta_file(path, &protein_to_go) {

                    let mut taxon_map = HashMap::new();
                    let mut taxon_protein_count = HashMap::new();
                    let mut go_term_count = HashMap::new();
                    let mut go_term_to_protein_set = HashMap::new();

                    taxon_protein_count.insert(taxon_id, protein_set.len());
                    taxon_map.insert(taxon_id, protein_set);
                    go_term_count.insert(taxon_id, go_term_count_map);
                    go_term_to_protein_set.insert(taxon_id, go_term_to_proteins);

                    return Ok(Some(StudyPop {
                        taxon_map,
                        taxon_protein_count,
                        go_term_count,
                        go_term_to_protein_set,
                    }));
            }}
            _ => {}
        }

        let entries: Vec<_> = read_dir(path)?
            .filter_map(Result::ok)
            .collect();

        let results: Vec<_> = entries.par_iter()
        .filter_map(|entry| {
            parse_fasta_file(entry.path(), &protein_to_go)
                .ok()
                .flatten()
        })
        .collect();

        let mut taxon_map = HashMap::new();
        let mut taxon_protein_count = HashMap::new();
        let mut go_term_count = HashMap::new();
        let mut go_term_to_protein_set = HashMap::new();

        for (taxon_id, protein_set, go_term_count_map, go_term_to_proteins) in results {
        taxon_protein_count.insert(taxon_id, protein_set.len());
        taxon_map.insert(taxon_id, protein_set);
        go_term_count.insert(taxon_id, go_term_count_map);
        go_term_to_protein_set.insert(taxon_id, go_term_to_proteins);
        }
        
        Ok(Some(StudyPop {
            taxon_map,
            taxon_protein_count,
            go_term_count,
            go_term_to_protein_set
        }))
    }

    pub fn filter_by_threshold(
        &mut self, 
        taxon_ids: &HashSet<TaxonID>,
        threshold: usize) {  
        
        for taxon_id in taxon_ids {
            let terms_to_remove: Vec<GOTermID> = if let Some(term_count) = self.go_term_count.get(&taxon_id) {
                term_count
                    .iter()
                    .filter_map(|(term_id, &count)| {
                        if count <= threshold {
                            Some(term_id.clone())
                        } else {
                            None
                        }
                    })
                    .collect()
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
        }
    }
}

pub fn parse_fasta_file(
        fasta_file_path: impl AsRef<Path>,
        protein_to_go: &HashMap<TaxonID, ProteinToGO>,
        ) -> Result<Option<(
            TaxonID, 
            HashSet<Protein>, 
            GOTermCount,
            HashMap<GOTermID, HashSet<Protein>>)>> {

    let path = fasta_file_path.as_ref();
    
    if !path.is_file() {
        return Ok(None);
    }

    let file = File::open(path)?;
    let reader = BufReader::with_capacity(128 * 1024, file);

    let mut taxon_id: TaxonID = TaxonID::default();
    let mut protein_set: HashSet<Protein> = HashSet::new();
    
    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed_line = line.trim();
        if !trimmed_line.starts_with(">"){
            return Err(Error::new(
                ErrorKind::InvalidInput,
                format!("File does not have correct formatting: {:?}", path)
            ));
        }
        else if trimmed_line.starts_with(">") {
            let id_str = &trimmed_line[1..].trim();
            taxon_id = id_str.parse::<TaxonID>()
                .map_err(|e| Error::new(
                    ErrorKind::InvalidData,
                    format!("Failed to parse taxon ID: {}", e)
                ))?;
        }
        else if !trimmed_line.is_empty() {
            protein_set.insert(trimmed_line.to_owned());
        }
    }
    let mut go_term_count: GOTermCount = HashMap::new();
    let mut go_term_to_proteins: HashMap<GOTermID, HashSet<Protein>> = HashMap::new();

    let taxon_protein_to_go = protein_to_go.get(&taxon_id).unwrap();

    for protein in &protein_set {
        if let Some(protein_go_terms) = taxon_protein_to_go.get(protein) {
            for go_term in protein_go_terms {
                go_term_to_proteins
                    .entry(go_term.clone())
                    .or_insert_with(HashSet::new)
                    .insert(protein.to_string());
            }
        }
    }

    for (go_term, proteins) in &go_term_to_proteins {
        go_term_count.insert(go_term.clone(), proteins.len());
    }

    Ok(Some((
        taxon_id, 
        protein_set,
        go_term_count,
        go_term_to_proteins
    )))
}

pub fn collect_taxon_ids(
    study_data: impl AsRef<Path>
) -> Result<HashSet<TaxonID>> {

    let path = study_data.as_ref();
    let mut taxon_ids = HashSet::new();
    
    if let Some("csv") = path.extension().and_then(|s| s.to_str()) {
        if path.is_file() {
            let file = File::open(path)?;
            let mut csv_reader = Reader::from_reader(file);
            
            let header_taxons: Vec<TaxonID> = csv_reader
                .headers()?
                .clone()
                .iter()
                .filter_map(|id| id.parse::<u32>().ok())
                .collect();
                
            taxon_ids.extend(header_taxons);
            return Ok(taxon_ids);
        }
    }
    
    if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
        if (ext == "fa" || ext == "fasta") && path.is_file() {
            if let Some(id) = extract_taxon_id_from_fasta(path)? {
                taxon_ids.insert(id);
                return Ok(taxon_ids);
            }
        }
    }
    
    if path.is_dir() {
        let entries: Vec<_> = read_dir(path)?
            .filter_map(Result::ok)
            .collect();
            
        for entry in entries {
            let entry_path = entry.path();
            if let Some(ext) = entry_path.extension().and_then(|s| s.to_str()) {
                if ext == "fa" || ext == "fasta" {
                    if let Some(id) = extract_taxon_id_from_fasta(&entry_path)? {
                        taxon_ids.insert(id);
                    }
                }
            }
        }
    }
    
    Ok(taxon_ids)
}

fn extract_taxon_id_from_fasta(
    fasta_file_path: impl AsRef<Path>
) -> Result<Option<TaxonID>> {

    let file = File::open(fasta_file_path.as_ref())?;
    let reader = BufReader::new(file);
    
    for line_result in reader.lines() {
        let line = line_result?;
        let trimmed_line = line.trim();
        
        if trimmed_line.starts_with(">") {
            let id_str = &trimmed_line[1..].trim();
            return match id_str.parse::<TaxonID>() {
                Ok(id) => Ok(Some(id)),
                Err(e) => Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("Failed to parse taxon ID: {}", e)
                ))
            };
        } else {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "File does not have correct FASTA formatting (missing '>' header)"
            ));
        }
    }
    
    Ok(None)
}
