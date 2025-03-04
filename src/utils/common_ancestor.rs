use crate::parsers::obo_parser::*;
use daggy::{NodeIndex, Walker};
use std::collections::{HashMap, HashSet};

pub fn collect_ancestry_path(
    graph: &OntologyGraph,
    node_idx: NodeIndex,
) -> AncestryPath {
    let mut path = Vec::new();
    let mut to_visit = vec![(node_idx, None)];
    
    while let Some((current_idx, relationship)) = to_visit.pop() {
        path.push((current_idx, relationship));
        
        let mut parents = graph.parents(current_idx);
        while let Some((edge_idx, parent_idx)) = parents.walk_next(&graph) {
            let relationship = graph.edge_weight(edge_idx).unwrap().clone();
            to_visit.push((parent_idx, Some(relationship)));
        }
    }
    path
}

pub fn find_common_ancestors(
    paths: &[AncestryPath],
    node_index_to_go_id: &HashMap<NodeIndex, u32>
) -> Vec<u32> {
    if paths.is_empty() {
        return Vec::new();
    }

    let path_sets: Vec<HashSet<NodeIndex>> = paths
        .iter()
        .map(|path| path.iter().map(|(idx, _)| *idx).collect())
        .collect();
    
    let common_nodes = path_sets
        .iter()
        .skip(1)
        .fold(path_sets[0].clone(), |acc, set| {
            acc.intersection(set).cloned().collect()
        });
    
    common_nodes
        .into_iter()
        .map(|idx| node_index_to_go_id[&idx])
        .collect()
}

pub fn find_first_common_ancestor(
    paths: &[AncestryPath],
    node_index_to_go_id: &HashMap<NodeIndex, u32>,
    graph: &OntologyGraph
) -> Option<u32> {
    if paths.is_empty() {
        return None;
    }

    let path_sets: Vec<HashSet<NodeIndex>> = paths
        .iter()
        .map(|path| path.iter().map(|(idx, _)| *idx).collect())
        .collect();
    
    let common_nodes: HashSet<NodeIndex> = path_sets
        .iter()
        .skip(1)
        .fold(path_sets[0].clone(), |acc, set| {
            acc.intersection(set).cloned().collect()
        });

    let first_intersection = common_nodes.iter()
        .find(|&&node_idx| {
            let children: HashSet<NodeIndex> = graph.children(node_idx)
                .iter(graph)
                .map(|(_, child_idx)| child_idx)
                .collect();
            
            children.intersection(&common_nodes).next().is_none()
        });

    first_intersection.map(|&idx| node_index_to_go_id[&idx])
}

pub fn topological_sort(
    graph: &OntologyGraph,
    root_ids: &[u32],
    go_id_to_node_index: &HashMap<u32, NodeIndex>,
    node_index_to_go_id: &HashMap<NodeIndex, u32>,
) -> Vec<u32> {
    let mut sorted = Vec::new();
    let mut visited = HashSet::new();
    
    fn visit(
        node_id: u32,
        graph: &OntologyGraph,
        visited: &mut HashSet<u32>,
        sorted: &mut Vec<u32>,
        go_id_to_node_index: &HashMap<u32, NodeIndex>,
        node_index_to_go_id: &HashMap<NodeIndex, u32>,
    ) {
        if visited.insert(node_id) {
            let node_idx = go_id_to_node_index[&node_id];
            let mut parents = graph.parents(node_idx);
            
            while let Some((_, parent_idx)) = parents.walk_next(graph) {
                let parent_id = node_index_to_go_id[&parent_idx];
                visit(parent_id, graph, visited, sorted, go_id_to_node_index, node_index_to_go_id);
            }
            
            sorted.push(node_id);
        }
    }
    
    for &root_id in root_ids {
        visit(root_id, graph, &mut visited, &mut sorted, go_id_to_node_index, node_index_to_go_id);
    }
    
    sorted
}

pub fn generate_mermaid_chart(
    graph: &OntologyGraph,
    root_ids: &[u32],
    go_id_to_node_index: &HashMap<u32, NodeIndex>,
    node_index_to_go_id: &HashMap<NodeIndex, u32>,
    obo_map: &OboMap,
    first_common_ancestor: Option<u32>,
) -> String {
    let mut mermaid = String::from("\
    %%{init: {'flowchart': {'diagramPadding': 10,'nodeSpacing':20,'rankSpacing':20}}}%%\n\
    flowchart LR\n\
    %% Style declarations\n\
    classDef BP fill:#FFAB8350,stroke:black\n\
    classDef MF fill:#6CA9E050,stroke:black\n\
    classDef CC fill:#87E06C50,stroke:black\n\
    classDef inputNode fill:#D9B0D250,stroke:black\n\
    classDef intersectingNode fill:#E0B46F75,stroke:black\n\
    classDef legendNode fill:white, stroke:black\n\
    classDef bothNode fill:#E1C2B7, stroke:black\n\
    \n\
    %% Edge legend\n\
    subgraph Legends[ ]\n\
        direction LR\n\
        style Legends fill:white,stroke:#333,stroke-width:0.0px\n\
        \n\
        subgraph Edge_legend [Edge Types]\n\
            direction LR\n\
            style Edge_legend fill:white,stroke:#333,stroke-width:1.5px\n\
            \n\
            Start1([Child]) -->|Is a| End1([Parent])\n\
            Start2([Child]) -->|Part of| End2([Parent])\n\
            Start3([Child]) -->|Regulates| End3([Parent])\n\
            Start4([Child]) -->|Positively Regulates| End4([Parent])\n\
            Start5([Child]) -->|Negatively Regulates| End5([Parent])\n\
            Start6([Child]) -->|Occurs in| End6([Parent])\n\
            \n\
            linkStyle 0 stroke:#000000,stroke-width:2px\n\
            linkStyle 1 stroke:#FF00FF,stroke-width:2px\n\
            linkStyle 2 stroke:#006BFF,stroke-width:2px\n\
            linkStyle 3 stroke:#00AC00,stroke-width:2px\n\
            linkStyle 4 stroke:#EF0000,stroke-width:2px\n\
            linkStyle 5 stroke:#FFA500,stroke-width:2px\n\
            \n\
            class Start1,Start2,Start3,Start4,Start5,Start6,End1,End2,End3,End4,End5,End6 legendNode\n\
        end\n\
        \n\
        subgraph Node_legend [Node Types]\n\
            direction BT\n\
            style Node_legend fill:white,stroke:#333,stroke-width:1.5px\n\
            \n\
            BP_Node[\"Biological process\"]:::BP\n\
            MF_Node[\"Molecular function\"]:::MF\n\
            CC_Node[\"Cellular component\"]:::CC\n\
            Input_Node[\"Input nodes\"]:::inputNode\n\
            Intersecting_Node[\"First intersect\"]:::intersectingNode\n\
            Both_Node[\"Input & Intersect\"]:::bothNode\n\
            \n\
            BP_Node ~~~ MF_Node\n\
            MF_Node ~~~ CC_Node\n\
            \n\
            Input_Node ~~~ Intersecting_Node\n\
            Intersecting_Node ~~~ Both_Node\n\
        end\n\
    end\n\n\
    subgraph Ontology [ ]\n\
    direction BT\n\
    style Ontology fill:white, stroke:white\n\
    %% Nodes\n");

    let mut processed_nodes = HashSet::new();
    let mut edges = Vec::new();

    fn process_node(
        node_idx: NodeIndex,
        graph: &OntologyGraph,
        node_index_to_go_id: &HashMap<NodeIndex, u32>,
        processed_nodes: &mut HashSet<u32>,
        edges: &mut Vec<(u32, u32, Relationship)>,
    ) {
        let mut to_visit = vec![node_idx];
        
        while let Some(current_idx) = to_visit.pop() {
            let go_id = node_index_to_go_id[&current_idx];
            
            if processed_nodes.insert(go_id) {
                let mut parents = graph.parents(current_idx);
                while let Some((edge_idx, parent_idx)) = parents.walk_next(&graph) {
                    let parent_id = node_index_to_go_id[&parent_idx];
                    let relationship = graph.edge_weight(edge_idx).unwrap().clone();
                    edges.push((go_id, parent_id, relationship));
                    to_visit.push(parent_idx);
                }
            }
        }
    }

    for &root_id in root_ids {
        let root_idx = go_id_to_node_index[&root_id];
        process_node(
            root_idx,
            graph,
            node_index_to_go_id,
            &mut processed_nodes,
            &mut edges,
        );
    }

    let sorted_nodes = topological_sort(
        graph, 
        root_ids, 
        go_id_to_node_index, 
        node_index_to_go_id
    );

    for &go_id in &sorted_nodes {
        let term = &obo_map[&go_id];
        let mut class_name = match term.namespace {
            NameSpace::BiologicalProcess => ":::BP",
            NameSpace::MolecularFunction => ":::MF",
            NameSpace::CellularComponent => ":::CC",
        };

        if root_ids.contains(&go_id) && Some(go_id) == first_common_ancestor {
            class_name = ":::bothNode";
        } else if root_ids.contains(&go_id) {
            class_name = ":::inputNode";
        } else if Some(go_id) == first_common_ancestor {
            class_name = ":::intersectingNode";
        }
        
        mermaid.push_str(&format!(
            "    GO{:07}(\"**GO:{:07}**<br>{}\"){}\n",
            go_id, go_id, term.name.replace("_", " ").replace("\"", "&quot;"), class_name
        ));
    }

    mermaid.push_str("\n    %% Connections\n");
    
    let mut is_a_edges = Vec::new();
    let mut part_of_edges = Vec::new();
    let mut regulates_edges = Vec::new();
    let mut positively_regulates_edges = Vec::new();
    let mut negatively_regulates_edges = Vec::new();
    let mut occurs_in_edges = Vec::new();
    
    let mut edge_index = 10;
    for (child_id, parent_id, relationship) in edges {
        mermaid.push_str(&format!(
            "    GO{:07} --> GO{:07}\n",
            child_id, parent_id
        ));
        
        match relationship {
            Relationship::IsA => is_a_edges.push(edge_index),
            Relationship::PartOf => part_of_edges.push(edge_index),
            Relationship::Regulates => regulates_edges.push(edge_index),
            Relationship::PositivelyRegulates => positively_regulates_edges.push(edge_index),
            Relationship::NegativelyRegulates => negatively_regulates_edges.push(edge_index),
            Relationship::OccursIn => occurs_in_edges.push(edge_index),
        }
        edge_index += 1;
    }

    mermaid.push_str("\n    %% Link styles\n");
    
    if !is_a_edges.is_empty() {
        mermaid.push_str(&format!("    linkStyle {} stroke:#000000,stroke-width:2px\n", 
            is_a_edges.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(",")));
    }
    
    if !part_of_edges.is_empty() {
        mermaid.push_str(&format!("    linkStyle {} stroke:#FF00FF,stroke-width:2px\n", 
            part_of_edges.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(",")));
    }
    
    if !regulates_edges.is_empty() {
        mermaid.push_str(&format!("    linkStyle {} stroke:#006BFF,stroke-width:2px\n", 
            regulates_edges.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(",")));
    }
    
    if !positively_regulates_edges.is_empty() {
        mermaid.push_str(&format!("    linkStyle {} stroke:#00AC00,stroke-width:2px\n", 
            positively_regulates_edges.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(",")));
    }
    
    if !negatively_regulates_edges.is_empty() {
        mermaid.push_str(&format!("    linkStyle {} stroke:#EF0000,stroke-width:2px\n", 
            negatively_regulates_edges.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(",")));
    }
    
    if !occurs_in_edges.is_empty() {
        mermaid.push_str(&format!("    linkStyle {} stroke:#FFA500,stroke-width:2px\n", 
            occurs_in_edges.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(",")));
    }

    mermaid.push_str("\n    %% Click actions\n");
    for &go_id in &processed_nodes {
        mermaid.push_str(&format!(
            "    click GO{:07} href \"https://www.ebi.ac.uk/QuickGO/term/GO:{:07}\" \"Click to view GO:{:07}\"\n",
            go_id, go_id, go_id
        ));
    }
    mermaid.push_str(
        "
    end\n
    Ontology~~~Legends");
    mermaid
}