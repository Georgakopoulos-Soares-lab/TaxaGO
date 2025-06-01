use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::{read_dir, File};
use std::io::{BufRead, BufReader, ErrorKind};
use std::path::PathBuf;
use std::sync::Arc;
use csv::Reader;  
use rayon::prelude::*;
use crate::parsers::background_parser::*;
use compact_str::CompactString;
use thiserror::Error;

#[derive(Debug, Default, Clone)]
pub struct StudyPop {
    pub taxon_map: FxHashMap<TaxonID, FxHashSet<Protein>>,
    pub taxon_protein_count: FxHashMap<TaxonID, usize>,
    pub go_term_count: FxHashMap<TaxonID, GOTermCount>,
    pub go_term_to_protein_set: FxHashMap<TaxonID, FxHashMap<GOTermID, FxHashSet<Protein>>>
}

#[derive(Error, Debug)]
pub enum StudyPopError {
    #[error("File not found: {0}")]
    FileNotFound(PathBuf),

    #[error("Invalid file extension for file: {0}. Expected .csv, .fa, or .fasta.")]
    InvalidFileExtension(PathBuf),

    #[error("FASTA file ({0}) must start with a '>' header line.")]
    FastaMissingInitialHeader(PathBuf),

    #[error("FASTA file ({0}) contains multiple '>' taxon ID header lines. Only one is permitted.")]
    FastaMultipleHeaders(PathBuf),
}

type BoxedResult<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync + 'static>>;

impl StudyPop {
    pub fn from_csv_file(
        csv_file: &PathBuf,
        protein_to_go: &FxHashMap<TaxonID, ProteinToGO>,
    ) -> BoxedResult<Option<Self>> {
        let file = match File::open(csv_file) {
            Ok(f) => f,
            Err(e) if e.kind() == ErrorKind::NotFound => {
                return Err(Box::new(StudyPopError::FileNotFound(csv_file.clone())));
            }
            Err(e) => return Err(Box::new(e)),
        };
        let mut csv_reader = Reader::from_reader(file);

        let mut taxon_map: FxHashMap<TaxonID, FxHashSet<Protein>> = FxHashMap::default();

        let taxon_ids: Vec<TaxonID> = csv_reader
            .headers()
            .map_err(|e| Box::new(e) as Box<dyn std::error::Error + Send + Sync + 'static>)? 
            .clone()
            .iter()
            .filter_map(|id| id.parse::<u32>().ok())
            .collect();

        for result in csv_reader.records() {
            let record = result.map_err(|e| Box::new(e) as Box<dyn std::error::Error + Send + Sync + 'static>)?; 
            for (index, protein_str) in record.iter().enumerate() {
                if let Some(&taxon_id) = taxon_ids.get(index) {
                    if !protein_str.is_empty() {
                        taxon_map
                            .entry(taxon_id)
                            .or_insert_with(FxHashSet::default)
                            .insert(Arc::new(CompactString::new(protein_str)));
                    }
                }
            }
        }

        if taxon_map.is_empty() {
            return Ok(None);
        }

        let mut taxon_protein_count = FxHashMap::with_capacity_and_hasher(
            taxon_map.len(),
            rustc_hash::FxBuildHasher::default()
        );
        let mut go_term_to_protein_set = FxHashMap::with_capacity_and_hasher(
            taxon_map.len(),
            rustc_hash::FxBuildHasher::default()
        );
        let mut go_term_count = FxHashMap::with_capacity_and_hasher(
            taxon_map.len(),
            rustc_hash::FxBuildHasher::default()
        );

        for (&taxon_id, proteins_for_taxon) in &taxon_map {
            taxon_protein_count.insert(taxon_id, proteins_for_taxon.len());

            if let Some(specific_taxon_protein_to_go) = protein_to_go.get(&taxon_id) {
                if !proteins_for_taxon.is_empty() {
                    let (
                        calculated_go_term_count,
                        calculated_go_term_to_protein_set_map
                    ) = build_go_term_associations(proteins_for_taxon, specific_taxon_protein_to_go);

                    if !calculated_go_term_count.is_empty() {
                        go_term_count.insert(taxon_id, calculated_go_term_count);
                        go_term_to_protein_set.insert(taxon_id, calculated_go_term_to_protein_set_map);
                    }
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
        study_data_path: &PathBuf,
        protein_to_go: &FxHashMap<TaxonID, ProteinToGO>
    ) -> BoxedResult<Option<Self>> {
        if !study_data_path.exists() {
            return Err(Box::new(StudyPopError::FileNotFound(study_data_path.clone())));
        }

        if study_data_path.is_file() {
            match study_data_path.extension().and_then(|s| s.to_str()) {
                Some("csv") => {
                    return StudyPop::from_csv_file(study_data_path, protein_to_go);
                }
                Some("fa") | Some("fasta") => {
                    match parse_fasta_file(study_data_path, protein_to_go)? {
                        Some((taxon_id, protein_set, go_term_count_map, go_term_to_proteins)) => {
                            let mut taxon_map = FxHashMap::default();
                            let mut taxon_protein_count = FxHashMap::default();
                            let mut go_term_count = FxHashMap::default();
                            let mut go_term_to_protein_set = FxHashMap::default();

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
                        }
                        None => return Ok(None), 
                    }
                }
                _ => return Err(Box::new(StudyPopError::InvalidFileExtension(study_data_path.clone()))),
            }
        } else if study_data_path.is_dir() {
            let results: Vec<_> = read_dir(study_data_path)
                .map_err(|e| Box::new(e) as Box<dyn std::error::Error + Send + Sync + 'static>)?
                .filter_map(std::result::Result::ok)
                .par_bridge()
                .filter_map(|entry| {
                    let path = entry.path();
                    if let Some(ext) = path.extension().and_then(|s| s.to_str()) {
                        if ext == "fa" || ext == "fasta" {
                            return parse_fasta_file(&path, protein_to_go).ok().flatten();
                        }
                    }
                    None
                })
                .collect();

            if results.is_empty() {
                return Ok(None);
            }

            let mut taxon_map = FxHashMap::default();
            let mut taxon_protein_count = FxHashMap::default();
            let mut go_term_count = FxHashMap::default();
            let mut go_term_to_protein_set = FxHashMap::default();

            for (taxon_id, protein_set, go_term_count_map, go_term_to_proteins) in results {
                taxon_protein_count.insert(taxon_id, protein_set.len());
                taxon_map.insert(taxon_id, protein_set);
                go_term_count.insert(taxon_id, go_term_count_map);
                go_term_to_protein_set.insert(taxon_id, go_term_to_proteins);
            }

            return Ok(Some(StudyPop {
                taxon_map,
                taxon_protein_count,
                go_term_count,
                go_term_to_protein_set
            }));
        } else {
             return Err(Box::new(std::io::Error::new(ErrorKind::Other, format!("Path is not a file or directory: {}", study_data_path.display()))));
        }
    }

    pub fn filter_by_threshold(
        &mut self,
        taxon_ids: &FxHashSet<TaxonID>,
        threshold: usize
    ) {
        for taxon_id in taxon_ids {
            let terms_to_remove: Vec<GOTermID> = if let Some(term_count) = self.go_term_count.get(taxon_id) {
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

            if let Some(count_map) = self.go_term_count.get_mut(taxon_id) {
                for term_id in &terms_to_remove {
                    count_map.remove(term_id);
                }
            }

            if let Some(term_map) = self.go_term_to_protein_set.get_mut(taxon_id) {
                for term_id in &terms_to_remove {
                    term_map.remove(term_id);
                }
            }
        }
    }
}

pub type ParsedFastaData = (
    TaxonID,
    FxHashSet<Protein>,
    GOTermCount,
    FxHashMap<GOTermID, FxHashSet<Protein>>,
);

pub fn parse_fasta_file(
    fasta_file_path: &PathBuf,
    protein_to_go: &FxHashMap<TaxonID, ProteinToGO>,
) -> BoxedResult<Option<ParsedFastaData>> {
    let file = match File::open(fasta_file_path) {
        Ok(f) => f,
        Err(e) if e.kind() == ErrorKind::NotFound => {
            return Err(Box::new(StudyPopError::FileNotFound(fasta_file_path.clone())));
        }
        Err(e) => return Err(Box::new(e)), 
    };

    let reader = BufReader::with_capacity(128 * 1024, file);
    let mut lines_iter = reader.lines();

    let mut taxon_id_from_file: Option<TaxonID> = None;
    let mut protein_set: FxHashSet<Protein> = FxHashSet::default();
    let mut first_header_processed = false;

    while let Some(line_result) = lines_iter.next() {
        let line = line_result.map_err(Box::new)?; 
        let trimmed_line = line.trim();

        if trimmed_line.is_empty() {
            continue;
        }

        if !trimmed_line.starts_with('>') {
            return Err(Box::new(StudyPopError::FastaMissingInitialHeader(fasta_file_path.clone())));
        }

        let id_str_content = &trimmed_line[1..].trim();
        if id_str_content.is_empty() {
            return Err(Box::new(std::io::Error::new(
                ErrorKind::InvalidData,
                format!("Empty Taxon ID header (e.g., '>') found in file: {}", fasta_file_path.display())
            )));
        }

        taxon_id_from_file = Some(id_str_content.parse::<TaxonID>().map_err(|e| {
             Box::new(std::io::Error::new( 
                ErrorKind::InvalidData,
                format!("Failed to parse Taxon ID '{}' from header line '{}' in file {}: {}", id_str_content, trimmed_line, fasta_file_path.display(), e)
            ))
        })?);
        first_header_processed = true;
        break; 
    }

    if !first_header_processed {
        return Ok(None);
    }

    let taxon_id = taxon_id_from_file.unwrap(); 

    while let Some(line_result) = lines_iter.next() {
        let line = line_result.map_err(Box::new)?;
        let trimmed_line = line.trim();

        if trimmed_line.is_empty() {
            continue;
        }

        if trimmed_line.starts_with('>') {
            return Err(Box::new(StudyPopError::FastaMultipleHeaders(fasta_file_path.clone())));
        }

        protein_set.insert(Arc::new(CompactString::new(trimmed_line)));
    }

    let mut final_go_term_count: GOTermCount = FxHashMap::default();
    let mut final_go_term_to_proteins_map: FxHashMap<GOTermID, FxHashSet<Protein>> = FxHashMap::default();

    if let Some(specific_taxon_protein_to_go) = protein_to_go.get(&taxon_id) {
        if !protein_set.is_empty() {
            let (
                calculated_go_term_count,
                calculated_go_term_to_protein_set_map
            ) = build_go_term_associations(&protein_set, specific_taxon_protein_to_go);

            final_go_term_count = calculated_go_term_count;
            final_go_term_to_proteins_map = calculated_go_term_to_protein_set_map;
        }
    }

    Ok(Some((
        taxon_id,
        protein_set,
        final_go_term_count,
        final_go_term_to_proteins_map,
    )))
}


fn extract_taxon_id_from_fasta(
    fasta_file_path: &PathBuf,
) -> BoxedResult<Option<TaxonID>> {
    let file = match File::open(fasta_file_path) {
        Ok(f) => f,
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
            return Err(Box::new(StudyPopError::FileNotFound(fasta_file_path.clone())));
        }
        Err(e) => return Err(Box::new(e)),
    };
    let reader = BufReader::new(file);
    let mut lines_iter = reader.lines();
    let mut taxon_id_to_return: Option<TaxonID> = None;
    let mut first_header_found_and_processed = false;

    while let Some(line_result) = lines_iter.next() {
        let line = line_result.map_err(Box::new)?;
        let trimmed_line = line.trim();

        if trimmed_line.is_empty() {
            continue;
        }

        if trimmed_line.starts_with('>') {
            if !first_header_found_and_processed {
                let id_str = &trimmed_line[1..].trim();
                if id_str.is_empty() {
                    return Err(Box::new(std::io::Error::new(
                        ErrorKind::InvalidData,
                        format!("Empty Taxon ID header in FASTA file: {}", fasta_file_path.display())
                    )));
                }
                let id = id_str.parse::<TaxonID>().map_err(|e| {
                    Box::new(std::io::Error::new(
                        ErrorKind::InvalidData,
                        format!("Failed to parse taxon ID '{}' in FASTA file {}: {}", id_str, fasta_file_path.display(), e)
                    ))
                })?;
                taxon_id_to_return = Some(id);
                first_header_found_and_processed = true;
            } else {
                return Err(Box::new(StudyPopError::FastaMultipleHeaders(fasta_file_path.clone())));
            }
        } else {
            if !first_header_found_and_processed {
                return Err(Box::new(StudyPopError::FastaMissingInitialHeader(fasta_file_path.clone())));
            }

        }
    }

    Ok(taxon_id_to_return)
}

pub fn collect_taxon_ids(
    study_data: &PathBuf,
) -> BoxedResult<FxHashSet<TaxonID>> {
    if !study_data.exists() {
        return Err(Box::new(StudyPopError::FileNotFound(study_data.clone())));
    }

    let mut taxon_ids = FxHashSet::default();

    if study_data.is_file() {
        match study_data.extension().and_then(|s| s.to_str()) {
            Some("csv") => {
                let file = File::open(study_data).map_err(|e| {
                    if e.kind() == ErrorKind::NotFound { 
                        Box::new(StudyPopError::FileNotFound(study_data.clone()))
                    } else {
                        Box::new(e) as Box<dyn std::error::Error + Send + Sync + 'static>
                    }
                })?;
                let mut csv_reader = Reader::from_reader(file);
                let header_taxons: Vec<TaxonID> = csv_reader
                    .headers()
                    .map_err(|e| Box::new(e) as Box<dyn std::error::Error + Send + Sync + 'static>)?
                    .clone()
                    .iter()
                    .filter_map(|id| id.parse::<u32>().ok())
                    .collect();
                taxon_ids.extend(header_taxons);
            }
            Some("fa") | Some("fasta") => {
                if let Some(id) = extract_taxon_id_from_fasta(study_data)? {
                    taxon_ids.insert(id);
                }
            }
            _ => return Err(Box::new(StudyPopError::InvalidFileExtension(study_data.clone()))),
        }
    } else if study_data.is_dir() {
        let entries: Vec<_> = read_dir(study_data)
            .map_err(|e| Box::new(e) as Box<dyn std::error::Error + Send + Sync + 'static>)?
            .filter_map(std::result::Result::ok) 
            .collect();

        for entry in entries {
            let entry_path = entry.path();
            if entry_path.is_file() {
                if let Some(ext) = entry_path.extension().and_then(|s| s.to_str()) {
                    if ext == "fa" || ext == "fasta" {
                        match extract_taxon_id_from_fasta(&entry_path) {
                            Ok(Some(id)) => { taxon_ids.insert(id); }
                            Ok(None) => {  }
                            Err(_e) => {  }
                        }
                    }
                }
            }
        }
    }

    Ok(taxon_ids)
}

fn build_go_term_associations(
    proteins: &FxHashSet<Protein>,
    taxon_protein_to_go: &ProteinToGO,
) -> (GOTermCount, FxHashMap<GOTermID, FxHashSet<Protein>>) {
    let mut go_term_to_protein_set_map: FxHashMap<GOTermID, FxHashSet<Protein>> = FxHashMap::default();

    for protein_arc in proteins {
        if let Some(go_terms_for_protein) = taxon_protein_to_go.get(protein_arc.as_ref()) {
            for &go_term_id in go_terms_for_protein {
                go_term_to_protein_set_map
                    .entry(go_term_id)
                    .or_insert_with(FxHashSet::default)
                    .insert(Arc::clone(protein_arc));
            }
        }
    }

    let mut go_term_count_map: GOTermCount = FxHashMap::default();
    if !go_term_to_protein_set_map.is_empty() {
        for (&go_term_id, associated_proteins) in &go_term_to_protein_set_map {
            go_term_count_map.insert(go_term_id, associated_proteins.len());
        }
    }

    (go_term_count_map, go_term_to_protein_set_map)
}