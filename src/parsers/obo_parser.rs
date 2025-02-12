use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufRead, Error};
use ucfirst::ucfirst;
use daggy::{Dag, NodeIndex};
use regex::Regex;
use lazy_static::lazy_static;

pub type OboMap = HashMap<u32, OboTerm>;
pub type OntologyGraph = Dag<u32, Relationship, u32>;
pub type AncestryPath = Vec<(NodeIndex, Option<Relationship>)>;

#[derive(Debug, Clone)]
pub enum NameSpace {
    BiologicalProcess,
    MolecularFunction,
    CellularComponent,
}

#[derive(Debug, Clone)]
pub struct OboTerm {
    pub name: String,
    pub namespace: NameSpace,
    pub definition: String,
    pub is_obsolete: bool,
    pub relationships: HashMap<u32, Relationship>,
}
impl OboTerm {
    pub fn new() -> Self {
        OboTerm {
            name: String::with_capacity(90),
            namespace: NameSpace::BiologicalProcess,
            definition: String::with_capacity(350),
            is_obsolete: false,
            relationships: HashMap::new(),
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

pub fn parse_obo_file(obo_file_path: &str) -> Result<OboMap, Error> {
    let mut obo_terms: HashMap<u32, OboTerm> = HashMap::with_capacity(41_000);

    let obo = File::open(obo_file_path)?;
    let reader = BufReader::with_capacity(3000 * 1024, obo);
    let lines = reader.lines();

    let mut new_term = false;
    let mut current_term = OboTerm::new();
    let mut current_id: u32 = 0;
    let mut obsolete_term: bool = false;

    for line in lines {
        let line = line?;
        
        if line == "[Term]" {
            new_term = true;
            current_term = OboTerm::new();
            obsolete_term = false; 
            current_id = 0; 

        } else if new_term {
            match line {
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
                        obo_terms.insert(current_id, current_term.clone());
                    }
                },
                _ => (), 
            }
        }
    }
    
    if new_term && current_id != 0 {
        obo_terms.insert(current_id, current_term);
    }
    
    Ok(obo_terms)
}

pub fn build_ontology_graph(obo_map: &OboMap) -> Result<(OntologyGraph, HashMap<u32, NodeIndex>), Error> {
    let mut ontology_graph: OntologyGraph = Dag::new();
    
    let mut node_indices: HashMap<u32, NodeIndex> = HashMap::new();
    
    for node_id in obo_map.keys() {
        let node_index = ontology_graph.add_node(*node_id);
        node_indices.insert(*node_id, node_index);
    }
    
    for (node_id, term) in obo_map.iter() {
        let source_index = node_indices[node_id];
        
        for (parent_id, relationship_type) in term.relationships.iter() {
            let target_index = node_indices[parent_id];
            ontology_graph
                .add_edge(target_index, source_index, relationship_type.clone())
                .unwrap();
        }
    }
    
    Ok((ontology_graph, node_indices))
}
