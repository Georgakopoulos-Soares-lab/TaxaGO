use std::collections::{HashMap, HashSet};
use std::fs::{read_dir, File};
use std::io::{BufRead, BufReader, Error, ErrorKind, Result};
use std::path::Path;
use csv::Reader; 
use crate::parsers::{
    obo_parser::*,
    background_parser::*,
};
use daggy::{NodeIndex, Walker};
use::rayon::prelude::*;

#[derive(Debug, Default, Clone)]
pub struct StudyPop {
    pub taxon_map: HashMap<TaxonID, HashSet<String>>,
    pub taxon_protein_count: HashMap<TaxonID, usize>,
    pub go_term_count: HashMap<TaxonID, GOTermCount>,
}

impl StudyPop {
    pub fn map_go_terms(&mut self, protein_to_go: &HashMap<TaxonID, ProteinToGO>) -> Result<()> {
        
        self.go_term_count = HashMap::new();
        for (&taxon_id, proteins) in &self.taxon_map {
            let mut go_term_counts = HashMap::new();
            if let Some(taxon_proteins) = protein_to_go.get(&taxon_id) {
                for protein in proteins {
                    if let Some(go_terms) = taxon_proteins.get(protein) {
                        for &go_term in go_terms { 
                            *go_term_counts.entry(go_term).or_insert(0) += 1;
                        }
                    }
                }
            }
            if !go_term_counts.is_empty() {
                self.go_term_count.insert(taxon_id, go_term_counts);
            }
        }
        Ok(())
    }

    pub fn get_go_term_counts(&self, taxon_id: &TaxonID) -> Option<&GOTermCount> {
        self.go_term_count.get(taxon_id)
    }

    pub fn get_go_term_count(&self, taxon_id: &TaxonID, go_term: u32) -> usize {
        self.go_term_count
            .get(taxon_id)
            .and_then(|terms| terms.get(&go_term))
            .copied()
            .unwrap_or(0)
    }

    pub fn propagate_counts(
        &mut self,
        graph: &OntologyGraph,
        go_id_to_node_index: &HashMap<u32, NodeIndex>
    ) {
        let ancestor_cache: HashMap<NodeIndex, HashSet<NodeIndex>> = go_id_to_node_index
            .values()
            .par_bridge()
            .map(|&node_idx| (node_idx, get_unique_ancestors(node_idx, graph)))
            .collect();

        let node_index_to_go_id: HashMap<NodeIndex, u32> = go_id_to_node_index
            .iter()
            .map(|(&go_id, &node_idx)| (node_idx, go_id))
            .collect();

        let updated_counts: HashMap<_, _> = self.go_term_count
            .par_iter()
            .map(|(&taxon_id, go_term_counts)| {
                let mut updated_counts = HashMap::with_capacity(go_term_counts.len() * 2);
                
                for (&go_id, &direct_count) in go_term_counts.iter() {
                    if let Some(&node_idx) = go_id_to_node_index.get(&go_id) {
                        *updated_counts.entry(go_id).or_insert(0) += direct_count;
                        
                        if let Some(ancestors) = ancestor_cache.get(&node_idx) {
                            for &ancestor_idx in ancestors {
                                if let Some(&ancestor_id) = node_index_to_go_id.get(&ancestor_idx) {
                                    *updated_counts.entry(ancestor_id).or_insert(0) += direct_count;
                                }
                            }
                        }
                    }
                }
                
                (taxon_id, updated_counts)
            })
            .collect();
        
        self.go_term_count = updated_counts;
    }

    
    pub fn from_csv<P: AsRef<Path>>(csv_path: P) -> Result<Self> {
        let mut taxon_map = HashMap::new();
        let mut taxon_protein_count = HashMap::new();
        
        let file = File::open(csv_path)?;
        let mut rdr = Reader::from_reader(file);
        

        let headers = rdr.headers()?.clone();
        let taxon_ids: Vec<TaxonID> = headers
            .iter()
            .filter_map(|id| id.parse::<u32>().ok())
            .collect();
            

        for result in rdr.records() {
            let record = result?;
            for (idx, protein) in record.iter().enumerate() {
                if let Some(&taxon_id) = taxon_ids.get(idx) {
                    if !protein.is_empty() {
                        taxon_map
                            .entry(taxon_id)
                            .or_insert_with(HashSet::new)
                            .insert(protein.to_string());
                    }
                }
            }
        }
        

        for (&taxon_id, proteins) in &taxon_map {
            taxon_protein_count.insert(taxon_id, proteins.len());
        }
        
        Ok(StudyPop {
            taxon_map,
            taxon_protein_count,
            go_term_count: HashMap::new(),
        })
    }
}

fn get_unique_ancestors(
    node_idx: NodeIndex,
    graph: &OntologyGraph,
) -> HashSet<NodeIndex> {
    let mut ancestors = HashSet::new();
    let mut to_visit = vec![node_idx];
    let mut visited = HashSet::new();
    
    while let Some(current_idx) = to_visit.pop() {
        if !visited.insert(current_idx) {
            continue;
        }
        
        let mut parents = graph.parents(current_idx);
        while let Some((edge_idx, parent_idx)) = parents.walk_next(graph) {
            match graph.edge_weight(edge_idx).unwrap() {
                Relationship::IsA | Relationship::PartOf => {
                    ancestors.insert(parent_idx);
                    to_visit.push(parent_idx);
                },
                _ => continue,
            }
        }
    }
    
    ancestors
}

fn is_fasta_file(path: &Path) -> bool {
    matches!(
        path.extension().and_then(|s| s.to_str()),
        Some("fasta" | "fa")
    )
}

fn process_fasta_file(file_path: impl AsRef<Path>) -> Result<Option<(TaxonID, HashSet<String>)>> {
    let path = file_path.as_ref();
    if !path.is_file() || !is_fasta_file(path) {
        return Ok(None);
    }

    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(32 * 1024, file).lines();
    
    let taxon_id = match reader.next() {
        Some(Ok(first_line)) => parse_taxon_id(&first_line)?,
        _ => return Err(Error::new(ErrorKind::InvalidData, "Empty or invalid FASTA file")),
    };

    let proteins: HashSet<String> = reader
        .filter_map(|line| line.ok())
        .map(|line| line.trim().to_string())
        .filter(|line| !line.is_empty())
        .collect();

    Ok(Some((taxon_id, proteins)))
}

fn parse_taxon_id(line: &str) -> Result<TaxonID> {
    let taxon_id_str = line
        .trim()
        .strip_prefix('>')
        .ok_or_else(|| Error::new(ErrorKind::InvalidData, "Invalid Taxon ID format"))?;
    
    taxon_id_str
        .parse::<u32>()
        .map_err(|e| Error::new(ErrorKind::InvalidData, e))
}

pub fn read_study_pop<P: AsRef<Path>>(path: P) -> Result<StudyPop> {
    let path = path.as_ref();
    
    // If the path is a file and ends with .csv, use CSV parser
    if path.is_file() && path.extension().and_then(|s| s.to_str()) == Some("csv") {
        return StudyPop::from_csv(path);
    }
    
    // Otherwise, use existing FASTA directory processing
    let mut taxon_map = HashMap::new();
    let mut taxon_protein_count = HashMap::new();
    
    for entry in read_dir(path)? {
        if let Ok(entry) = entry {
            if let Some((taxon_id, proteins)) = process_fasta_file(entry.path())? {
                taxon_protein_count.insert(taxon_id, proteins.len());
                taxon_map.insert(taxon_id, proteins);
            }
        }
    }
    
    Ok(StudyPop {
        taxon_map,
        taxon_protein_count,
        go_term_count: HashMap::new(),
    })
}