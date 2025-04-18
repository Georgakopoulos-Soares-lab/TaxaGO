use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Result};
use std::path::Path;
use crate::analysis::enrichment_analysis::*;
use crate::parsers::background_parser::TaxonID;

pub fn read_lineage<P: AsRef<Path>>(path: P) -> Result<FxHashMap<TaxonID, Vec<String>>> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(32 * 1024, file);
    let mut taxonomy = FxHashMap::default();

    for line in reader.lines().skip(1) {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        
        if let Ok(taxon_id) = fields[0].parse::<u32>() {
            let lineage = fields[2..].iter()
                .map(|&s| s.to_string())
                .collect();
            taxonomy.insert(taxon_id, lineage);
        }
    }
    Ok(taxonomy)
}
pub fn taxid_to_species<P: AsRef<Path>>(path: P) -> Result<FxHashMap<TaxonID, String>> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(32 * 1024, file);
    let mut taxid_species_map = FxHashMap::default();

    for line in reader.lines().skip(1) {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();

        if let Ok(taxon_id) = fields[0].parse::<u32>() {
            let species_name = fields[1].to_string();
            taxid_species_map.insert(taxon_id, species_name);
        }
    }
    Ok(taxid_species_map)
}

pub fn taxid_to_level(
    significant_results: &FxHashMap<u32, FxHashMap<u32, GOTermResults>>,
    taxonomic_lineage: &FxHashMap<TaxonID, Vec<String>>,
    taxonomic_level: &str
) -> FxHashMap<String, Vec<TaxonID>> {
    let mut grouped_results: FxHashMap<String, Vec<u32>> = FxHashMap::default();
    
    let level_index = match taxonomic_level.to_lowercase().as_str() {
        "genus" => 0,
        "family" => 1,
        "order" => 2,
        "class" => 3,
        "phylum" => 4,
        "kingdom" => 5,
        "superkingdom" => 6,
        _ => return grouped_results,
    };

    for (taxon_id, _) in significant_results {
        if let Some(lineage) = taxonomic_lineage.get(taxon_id) {
            if let Some(level_name) = lineage.get(level_index) {
                grouped_results
                    .entry(level_name.to_string())
                    .or_default()
                    .push(*taxon_id);
            }
        }
    }

    grouped_results
}

pub fn read_taxon_organism_count <P: AsRef<Path>>(path: P) -> Result<FxHashMap<String, usize>> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(16 * 1024,file);
    let mut species_counts = FxHashMap::with_capacity_and_hasher(
        35000,
        rustc_hash::FxBuildHasher::default()
    );
    
    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        
        if parts.len() >= 2 {
            let species_name = parts[..parts.len()-1].join("\t").trim().to_string();
            if let Ok(count) = parts.last().unwrap().trim().parse::<usize>() {
                species_counts.insert(species_name, count);
            }
        }
    }
    
    Ok(species_counts)
}