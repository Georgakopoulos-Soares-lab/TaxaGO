use std::collections::VecDeque;
use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{BufReader, BufRead, Error};
use std::mem;
use std::path::PathBuf;
use ucfirst::ucfirst;
use daggy::{Dag, NodeIndex, Walker};
use regex::Regex;
use lazy_static::lazy_static;
use strum_macros::EnumIter; 
use thiserror::Error;

use super::background_parser::*;

pub type OboMap = FxHashMap<u32, OboTerm>;
pub type OntologyGraph = Dag<u32, Relationship, u32>;
pub type AncestryPath = Vec<(NodeIndex, Option<Relationship>)>;

pub type TermToLevel = FxHashMap<GOTermID, usize>;
pub type LevelToTerms = FxHashMap<usize, Vec<GOTermID>>;

#[derive(Default, Debug, Copy, Clone, Eq, Hash, PartialEq, EnumIter)]
pub enum NameSpace {
    #[default]
    BiologicalProcess,
    MolecularFunction,
    CellularComponent,
}

#[derive(Default, Debug, Clone)]
pub struct OboTerm {
    pub name: String,
    pub namespace: NameSpace,
    pub definition: String,
    pub is_obsolete: bool,
    pub relationships: FxHashMap<u32, Relationship>,
}
impl OboTerm {
    pub fn new() -> Self {
        OboTerm {
            name: String::with_capacity(90),
            namespace: NameSpace::BiologicalProcess,
            definition: String::with_capacity(350),
            is_obsolete: false,
            relationships: FxHashMap::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub enum Relationship {
    IsA,
    PartOf,
    OccursIn,
    Regulates,
    PositivelyRegulates,
    NegativelyRegulates,   
}

#[derive(Error, Debug)]
pub enum OboParserError {
    #[error("OBO file not found at path: '{filepath}'. Please ensure the file exists and the path is correct.")]
    FileNotFound { filepath: String }, 
    
    #[error("OBO I/O Error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Invalid file extension for file '{filename}'. Expected '.obo', but found '.{extension_found}'.")]
    InvalidFileExtension {
        filename: String, 
        extension_found: String 
    }
}

lazy_static! {
    static ref GO_ID_REGEX: Regex = Regex::new(r"GO:(\d{7})").unwrap();
}

lazy_static! {
    static ref RELATIONSHIP_REGEX: Regex = 
        Regex::new(r"relationship:\s+(\w+(?:_\w+)*)\s+GO:(\d+)(?:\s+!.*)?").unwrap();
}

lazy_static! {
    static ref INTERSECTION_REGEX: Regex = 
        Regex::new(r"intersection_of:\s+(\w+(?:_\w+)*)\s+GO:(\d+)(?:\s+!.*)?").unwrap();
}

pub fn parse_is_a(input: &str) -> Option<u32> {
    GO_ID_REGEX
        .captures(input)
        .and_then(|caps| caps.get(1))
        .and_then(|m| m.as_str().parse().ok())
}

fn parse_relationship_type(relation: &str) -> Option<Relationship> {
    match relation {
        "part_of" => Some(Relationship::PartOf),
        "occurs_in" => Some(Relationship::OccursIn),
        "regulates" => Some(Relationship::Regulates),
        "positively_regulates" => Some(Relationship::PositivelyRegulates),
        "negatively_regulates" => Some(Relationship::NegativelyRegulates),
        "has_part" | "happens_during" | "ends_during" | _ => None,
    }
}

pub fn parse_relationship(input: &str, regex: &Regex) -> Option<(u32, Relationship)> {
    let caps = regex.captures(input)?;
    
    let relation = caps.get(1)?.as_str();
    let relationship = parse_relationship_type(relation)?;
    
    let id_group = 2;
    let id = caps.get(id_group)?
        .as_str()
        .parse::<u32>()
        .ok()?;
    
    Some((id, relationship))
}

pub fn parse_obo_file(obo_file_path: &PathBuf) -> Result<OboMap, OboParserError> {

    if !obo_file_path.exists() {
        return Err(OboParserError::FileNotFound { 
            filepath: obo_file_path.display().to_string() 
        });
    }

    match obo_file_path.extension().and_then(std::ffi::OsStr::to_str) {
        Some(ext) if ext.eq_ignore_ascii_case("obo") => {
        }
        Some(ext) => {
            return Err(OboParserError::InvalidFileExtension {
                filename: obo_file_path.display().to_string(),
                extension_found: ext.to_string(),
            });
        }
        None => {
            return Err(OboParserError::InvalidFileExtension {
                filename: obo_file_path.display().to_string(),
                extension_found: "No extension".to_string(),
            });
        }
    }

    let mut obo_terms: FxHashMap<u32, OboTerm> = FxHashMap::with_capacity_and_hasher(
        41_000, 
        rustc_hash::FxBuildHasher::default()
    );

    let obo = File::open(obo_file_path)?;
    let reader = BufReader::with_capacity(3000 * 1024, obo);

    let mut new_term = false;
    let mut current_term = OboTerm::new();
    let mut current_id: u32 = 0;
    let mut obsolete_term: bool = false;

    for line in reader.lines() {
        let line = line?;
        
        if line == "[Term]" {
            new_term = true;
            current_term = OboTerm::new();
            obsolete_term = false; 
            current_id = 0; 

        } else if new_term {
            match line.trim() {
                line if line.starts_with("id: ") => {
                    let id: u32 = line.split("GO:")
                        .nth(1)
                        .expect("No GO: ID found")
                        .parse()
                        .expect("Invalid GO number");
                    current_id = id;
                },
                line if line.starts_with("name: ") => {
                    let name = ucfirst(line.split(": ")
                        .nth(1)
                        .expect("No name found"));
                    current_term.name = name;
                },
                line if line.starts_with("namespace: ") => {
                    let namespace = match line.split(": ")
                        .nth(1)
                        .expect("No namespace found") {
                            "biological_process" => NameSpace::BiologicalProcess,
                            "molecular_function" => NameSpace::MolecularFunction,
                            "cellular_component" => NameSpace::CellularComponent,
                            _ => panic!("Invalid namespace"),
                    };
                    current_term.namespace = namespace;
                },
                line if line.starts_with("def: ") => {
                    let definition =  line.split("\"")
                        .nth(1)
                        .expect("No definition found");

                    current_term.definition = definition.to_string();
                },
                line if line.starts_with("is_a: ") => {
                    let parent_id: u32 = parse_is_a(&line)
                        .expect("No parent ID found");
                    current_term.relationships.insert(parent_id, Relationship::IsA);
                },
                line if line.starts_with("is_obsolete: ") => {
                    let is_obsolete = match line.split(": ")
                        .nth(1)
                        .expect("Wrong is_obsolete format") {
                            "true" => true,
                            _ => false,
                    };
                    obsolete_term=is_obsolete;
                    current_term.is_obsolete = is_obsolete;
                },
                line if line.starts_with("relationship: ") => {
                    if let Some((id, relationship)) = parse_relationship(&line, &RELATIONSHIP_REGEX) {
                        current_term.relationships.insert(id, relationship);
                    }
                },
                line if line.starts_with("intersection_of: ") => {
                    if !line.contains("intersection_of: GO:") {
                        if let Some((id, relationship)) = parse_relationship(&line, &INTERSECTION_REGEX) {
                            current_term.relationships.insert(id, relationship);
                        }
                    } else {
                        continue;
                    }
                    
                }
                line if line.is_empty() => {
                    new_term = false;
                    if current_id != 0 && obsolete_term == false {  
                        obo_terms.insert(current_id, mem::take(&mut current_term));
                    }
                },
                _ => (), 
            }
        }
    }
    
    if new_term && current_id != 0 {
        obo_terms.insert(current_id, mem::take(&mut current_term));
    }
    
    Ok(obo_terms)
}

pub fn build_ontology_graph(obo_map: &OboMap) -> Result<(OntologyGraph, FxHashMap<u32, NodeIndex>), Error> {
    let mut ontology_graph: OntologyGraph = Dag::new();
    
    let mut go_id_to_node_index: FxHashMap<u32, NodeIndex> = FxHashMap::default();
    
    for node_id in obo_map.keys() {
        let node_index = ontology_graph.add_node(*node_id);
        go_id_to_node_index.insert(*node_id, node_index);
    }
    
    for (node_id, term) in obo_map.iter() {
        let source_index = go_id_to_node_index[node_id];
        
        for (parent_id, relationship_type) in term.relationships.iter() {
            let target_index = go_id_to_node_index[parent_id];
            ontology_graph
                .add_edge(target_index, source_index, relationship_type.clone())
                .unwrap();
        }
    }
    
    Ok((ontology_graph, go_id_to_node_index))
}
pub fn assign_levels_from_roots(
    graph: &OntologyGraph, 
    go_id_to_node_index: &FxHashMap<GOTermID, NodeIndex>,
    node_index_to_go_id: &FxHashMap<NodeIndex, GOTermID>,
    root_ids: &[u32]  
) -> (TermToLevel, LevelToTerms) {

    let mut term_to_level = FxHashMap::default();
    let mut level_to_terms = FxHashMap::default();

    for &term_id in go_id_to_node_index.keys() {
        term_to_level.insert(term_id, 0);
    }

    let mut queue = VecDeque::new();

    for &root_id in root_ids {
        if let Some(&node_index) = go_id_to_node_index.get(&root_id) {
            term_to_level.insert(root_id, 0);
            queue.push_back((node_index, 0));
        }
    }

    while let Some((node_index, current_level)) = queue.pop_front() {
        let term_id = *node_index_to_go_id.get(&node_index).unwrap();

        level_to_terms.entry(current_level)
            .or_insert_with(Vec::new)
            .push(term_id);

        let mut children = graph.children(node_index);
        while let Some((edge, child)) = children.walk_next(graph) {
            match graph.edge_weight(edge).unwrap() {
                Relationship::IsA | Relationship::PartOf => {
                    let child_id = *node_index_to_go_id.get(&child).unwrap();
                    let child_current_level = term_to_level[&child_id];
                    
                    if current_level + 1 > child_current_level {
                        term_to_level.insert(child_id, current_level + 1);
                        queue.push_back((child, current_level + 1));
                    }
                },
                _ => continue,
            }
        }
    }

    (term_to_level, level_to_terms)
}
