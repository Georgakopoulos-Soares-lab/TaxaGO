import os, argparse, textwrap, math, shutil, concurrent.futures
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from typing import *
import networkx as nx

def load_results(
        directory: Path) -> Dict[str, pd.DataFrame]:
    
    """
    Loads TaxaGO enrichment analysis results into a dictionary. 
    If analysis was performed on taxonomic level, the keys will be the taxonomy name, else the keys will be the species name. 
    Values are Pandas dataframes containing the corresponding results.
    """

    results_dict = dict() 

    for file in directory.glob("*_GOEA_results.txt"):
            goea_results_df = pd.read_csv(file, sep='\t')
            file_name = file.stem.replace("_GOEA_results", "").replace("_", " ")
            results_dict[file_name] = goea_results_df

    return results_dict

def get_default_path(
        file:str) -> Path:
    
    """
    Identifies the location where Rust's Cargo is installed and then creates the basic TaxaGO assets directory.
    """

    cargo_home = os.environ.get("CARGO_HOME")
    
    if cargo_home is None:
        home = Path(os.path.expanduser("~"))
        cargo_home = home / ".cargo"

    default_lineage = cargo_home / "taxago_assets" / file
    
    return default_lineage

def get_default_lineage() -> Path:

    """
    Helper function to get the default lineage file from the provided taxago_assets directory.
    """

    return get_default_path("full_lineage.txt")

def get_default_background() -> Path:

    """
    Helper function to get the default background population directory from the provided taxago_assets directory.
    """

    return get_default_path("background_pop")

def create_lineage(
        path: Path) -> Dict[str, Dict[str, str]]:
    
    """
    Creates a nested dictionary representing taxonomic lineage from NCBI TaxonID up to species superkingdom.
    """

    lineage_dict = dict()

    with open(path, 'r') as file:
                lines = file.readlines()
                for line in lines[1:]:
                    line = line.strip().split('\t')
                    lineage_dict[line[0]] = {
                        "Species": line[1],
                        "Genus": line[2],
                        "Family": line[3],
                        "Order": line[4],
                        "Class": line[5],
                        "Phylum": line[6],
                        "Kingdom": line[7],
                        "Superkingdom": line[8]}

    return lineage_dict

def get_taxon_ids(
        lineage: Dict[str, Dict[str, str]], 
        taxonomy_names: List[str]) -> Dict[str, List[str]]:

    """
    Takes a list of taxonomy names and returns a dictionary where:
    - Each key is a taxonomy name from the input list
    - Each value is a list of taxon IDs associated with that taxonomy name
    """

    results = {name: [] for name in taxonomy_names}
    
    for taxon_id, lineage_tree in lineage.items():
        for _, taxonomic_name in lineage_tree.items():
            if taxonomic_name in taxonomy_names:
                if taxon_id not in results[taxonomic_name]:
                    results[taxonomic_name].append(taxon_id)
                break
    
    return results

def parse_fasta_file(
          fasta_file: Path,
          study_dict: Dict[str, Set[str]],
          lineage: Dict[str, Dict[str, str]]) -> Tuple[str, Dict[str, Set[str]]]:

     """
     Parses a FASTA file containing the study population used in the analysis for a specific species, where NCBI TaxonID is the FASTA header followed by proteins used in the analysis.
     Stores the results into a dictionary with the species name as key and the proteins present in the study population as a set. Also, stores the parsed NCBI TaxonID as a string.
     Returns a tuple with the species NCBI TaxonID and the study population dictionary.
     """

     protein_set = set()
     with open(fasta_file, 'r') as file:
          lines = file.readlines()
          for line in lines:
               line = line.strip()
               if line.startswith(">"):
                   taxon_id = line.replace(">", "")
                   species_name = lineage[taxon_id].get("Species")
               else:
                    protein_set.add(line) 
     study_dict[species_name] = protein_set

     return (taxon_id, study_dict)

def parse_csv_file(
          csv_file: Path,
          study_dict: Dict[str, Set[str]],
          lineage: Dict[str, Dict[str, str]]) -> Tuple[List[str], Dict[str, Set[str]]]:

    """
    Parses a CSV file containing the study populations for all species, where NCBI TaxonIDs are column names and rows are the proteins used in the analysis.
    Stores the results into a dictionary with the species name as key and the proteins present in the study population as a set. Also, stores the parsed NCBI TaxonIDs in a list.
    Returns a tuple of the list with the NCBI TaxonIDs and the study population dictionary.
    """

    study_pop_df = pd.read_csv(csv_file, sep = ',', dtype=str)
    taxon_ids = study_pop_df.columns.tolist()

    for taxon_id in study_pop_df.columns:
        if taxon_id not in lineage.keys():
            print(f"Warning: Taxon ID {taxon_id} not found in taxonomic lineage")
            continue
        species_name = lineage[taxon_id].get("Species")
        protein_set= set(study_pop_df[taxon_id].dropna().tolist())
        study_dict[species_name]=protein_set
    
    return (taxon_ids, study_dict)

def parse_background_pop(
        species_background_file: Path,
        study_proteins: Set[str]):

    """
    Parses the background populations for each species, filtering to only keep those present in the study population.
    Stores the results into a Pandas DataFrame, with 2 columns one with the GO terms and one with the sets of proteins with only those present in the study population.
    """

    background = pd.read_csv(species_background_file, sep='\t', header=None, names=['Protein', "GO_term"])
    filtered_background = background[background['Protein'].isin(study_proteins)]
    grouped_background = filtered_background.groupby('GO_term')['Protein'].apply(set).reset_index()

    return grouped_background

def process_taxonomy_networks(
        taxonomy: str, 
        taxonomy_id_map: Dict[str, List[str]],
        background: Path, 
        lineage: Dict[str, Dict[str, str]], 
        study_dict: Dict[str, Set[str]], 
        results_dict: Dict[str, pd.DataFrame]) -> Tuple[str, Dict[str, Dict[str, Set[str]]]]:
    
    """
    Process a single taxonomy and return its network data.
    """

    taxonomy_data = {}
    
    for namespace in ['Biological Process', 'Molecular Function', 'Cellular Component']:
        taxonomy_data[namespace] = {}
    
    taxon_ids_for_taxonomy = taxonomy_id_map[taxonomy]

    for taxon_id in taxon_ids_for_taxonomy:
        species_background_path = background / f"{taxon_id}_background.txt"
        species_name = lineage[taxon_id].get("Species")
        if species_name not in study_dict.keys():
            continue
        study_proteins = study_dict[species_name]
        
        original_species_background = parse_background_pop(
            species_background_path, 
            study_proteins)
        
        for namespace in ['Biological Process', 'Molecular Function', 'Cellular Component']:
            namespace_results = results_dict[taxonomy][results_dict[taxonomy]["Namespace"] == namespace]
            enriched_terms = namespace_results['GO Term ID'].unique().tolist()
            
            species_background = original_species_background[original_species_background['GO_term'].isin(enriched_terms)]
            
            for _, row in species_background.iterrows():
                go_term = row['GO_term']
                proteins = row['Protein']
                
                if go_term not in taxonomy_data[namespace]:
                    taxonomy_data[namespace][go_term] = set()
                
                taxonomy_data[namespace][go_term].update(proteins)
    
    return taxonomy, taxonomy_data

def process_species_networks(
        taxon_id: str, 
        background: Path, 
        lineage: Dict[str, Dict[str, str]], 
        study_dict: Dict[str, Set[str]], 
        results_dict: Dict[str, pd.DataFrame]) -> Tuple[str, Dict[str, Dict[str, Set[str]]]]:
    
    """
    Process a single species and return its network data.
    """

    species_background_path = background / f"{taxon_id}_background.txt"
    species_name = lineage[taxon_id].get("Species")
    
    if species_name not in study_dict:
        return species_name, {}
    
    species_data = {}
    study_proteins = study_dict[species_name]
    
    if species_name not in results_dict:
        return species_name, {}
        
    species_results = results_dict[species_name]
    
    original_species_background = parse_background_pop(
        species_background_path, 
        study_proteins)
    
    for namespace in ['Biological Process', 'Molecular Function', 'Cellular Component']:
        namespace_results = species_results[species_results["Namespace"] == namespace]
        enriched_terms = namespace_results['GO Term ID'].unique().tolist()
        species_background = original_species_background[original_species_background['GO_term'].isin(enriched_terms)]

        species_background_dict = dict(zip(species_background['GO_term'], species_background['Protein']))

        if namespace not in species_data:
            species_data[namespace] = {}
        
        species_data[namespace] = species_background_dict
    
    return species_name, species_data

def create_network_data(
          study_data: Path,
          background: Path, 
          lineage: Dict[str, Dict[str, str]],
          taxonomic_level: str,
          taxonomy_id_map: Dict[str, List[str]], 
          results_dict: Dict[str, pd.DataFrame],
          ) -> Dict[str, Dict[str, Dict[str, Set[str]]]]:

    """
    Creates a structured network dataset linking proteins, GO terms, and taxonomic information for visualization and analysis of gene ontology enrichment results.
    This function processes study data (proteins) and background GO term associations to build a hierarchical network structure. 
    Depending on the taxonomic_level parameter, it organizes data either by individual species or by higher-level taxonomic groups.
    The final dictionary goes from taxonomy/species to namespace to GO terms (only those found enriched) to a set of proteins (only those present in the study population). When analyzing taxonomies then the protein set contains proteins across all species associated with the taxonomy.
    """

    study_data = Path(study_data)  
    study_dict = dict()
    taxon_ids = list()
    network_data = dict()

    if study_data.is_dir():
         fasta_files = list(study_data.glob("*.fa")) + list(study_data.glob("*.fasta"))
         if fasta_files:
              for fasta_file in fasta_files:
                taxon_id, study_dict = parse_fasta_file(
                    fasta_file,
                    study_dict,
                    lineage)
                taxon_ids.append(taxon_id)
         else:
            print(f"No files with .fa/.fasta extension found in {study_data}")
            exit()
    elif study_data.is_file():
        if study_data.suffix.lower() in [".fa", ".fasta"]:
            taxon_id, study_dict = parse_fasta_file(
                 study_data,
                 study_dict,
                 lineage)
            taxon_ids.append(taxon_id)
        elif study_data.suffix.lower() == ".csv":
            taxon_ids, study_dict = parse_csv_file(
                study_data,
                study_dict,
                lineage)
        else:
            print(f"{study_data.name} is not a supported file type (must be FASTA or CSV)")
            exit()
    else:
        print(f"No valid study data found in {study_data}. Must be FASTA, CSV or a directory with FASTA.")
        exit()

    if taxonomic_level == "taxonomy":
        for taxonomy in results_dict.keys():
            taxonomy, taxonomy_data = process_taxonomy_networks(
                taxonomy, 
                taxonomy_id_map,
                background,
                lineage,
                study_dict,
                results_dict
            )
            network_data[taxonomy] = taxonomy_data

    elif taxonomic_level == "species":
        for taxon_id in taxon_ids:
            species_name, species_data = process_species_networks(
                taxon_id,
                background,
                lineage,
                study_dict,
                results_dict
            )
            if species_data: 
                network_data[species_name] = species_data
    
    return network_data

def process_taxon_network(
        taxon_name: str, 
        taxon_network_data: Dict[str, Dict[str, Set[str]]], 
        taxon_enrichment_results: pd.DataFrame, 
        similarity_threshold: float)-> Tuple[str, Dict[str, nx.Graph]]:
    
    """
    Process network data for a single taxon.
    This helper function builds GO term networks for a specific taxon across all namespaces.
    """

    taxon_networks = {}
    
    namespace_term_data = dict()
    for namespace in ['Biological Process', 'Molecular Function', 'Cellular Component']:
        namespace_results = taxon_enrichment_results[taxon_enrichment_results['Namespace'] == namespace]
        namespace_term_data[namespace] = {
            row['GO Term ID']: row 
            for _, row in namespace_results.iterrows()
        }
    
    for namespace in ['Biological Process', 'Molecular Function', 'Cellular Component']:
        if namespace not in taxon_network_data:
            continue
            
        go_term_proteins = taxon_network_data[namespace]
        namespace_results_dict = namespace_term_data[namespace]

        term_network = nx.Graph()

        term_to_proteins = dict()
        for go_term_id, proteins in go_term_proteins.items():
            if go_term_id in namespace_results_dict:
                node_enrichment = {'log(Odds Ratio)': namespace_results_dict[go_term_id]['log(Odds Ratio)']}
                term_network.add_node(go_term_id, **node_enrichment)
                term_to_proteins[go_term_id] = proteins

        protein_to_terms = dict()
        term_sizes = dict()
        
        for term_id, proteins in term_to_proteins.items():
            term_sizes[term_id] = len(proteins)
            for protein in proteins:
                if protein not in protein_to_terms:
                    protein_to_terms[protein] = set()
                protein_to_terms[protein].add(term_id)
        
        candidate_pairs = set()
        for terms in protein_to_terms.values():
            terms = list(terms)
            for i in range(len(terms)):
                for j in range(i+1, len(terms)):
                    candidate_pairs.add(tuple(sorted([terms[i], terms[j]])))
        
        edges = []
        for term1_id, term2_id in candidate_pairs:
            size1 = term_sizes[term1_id]
            size2 = term_sizes[term2_id]
            
            max_possible_similarity = min(size1, size2) / max(size1, size2)
            if max_possible_similarity < similarity_threshold:
                continue
            
            proteins1 = term_to_proteins[term1_id]
            proteins2 = term_to_proteins[term2_id]
            
            intersection_size = len(proteins1.intersection(proteins2))
            if intersection_size == 0:
                continue
                
            union_size = size1 + size2 - intersection_size
            similarity = round(intersection_size / union_size, 4)
            
            if similarity >= similarity_threshold:
                edges.append((term1_id, term2_id, similarity))
        
        term_network.add_weighted_edges_from(edges)
        taxon_networks[namespace] = term_network
    
    return taxon_name, taxon_networks

def build_go_term_networks(
        network_data: Dict[str, Dict[str, Dict[str, Set[str]]]], 
        enrichment_results: Dict[str, pd.DataFrame], 
        similarity_threshold: float,
        executor: concurrent.futures.ProcessPoolExecutor) -> Dict[str, Dict[str, nx.Graph]]:
    
    """
    Constructs networks of GO terms based on protein overlap between terms.
    This function builds a graph-based representation of relationships between GO terms for each taxon (species or taxonomic group) and each GO namespace. 
    The relationships are determined by calculating the similarity between pairs of GO terms based on their
    shared proteins. 
    Only relationships that meet a minimum similarity threshold are included.
    """
    
    go_networks = dict()

    futures = []
    
    for taxon_name in network_data.keys():
        if taxon_name not in enrichment_results:
            continue
            
        futures.append(
            executor.submit(
                process_taxon_network,
                taxon_name,
                network_data[taxon_name],
                enrichment_results[taxon_name],
                similarity_threshold
            )
        )

    for future in concurrent.futures.as_completed(futures):
        try:
            taxon_name, taxon_networks = future.result()
            go_networks[taxon_name] = taxon_networks
        except Exception as e:
            print(f"Error processing network for taxon: {e}")

    
    return go_networks

def process_taxon_subnetworks(
        taxon_name: str, 
        taxon_networks: Dict[str, nx.Graph], 
        enrichment_threshold: float, 
        number_of_subnetworks: int)-> Tuple[str, Dict[str, List[nx.Graph]]]:
   
    """
    Process subnetworks for a single taxon across all namespaces.
    
    This function:
    1. Filters each network to keep only highly enriched nodes
    2. Finds connected components in the filtered networks
    3. Selects the largest and most informative subnetworks
    4. Returns subgraphs for visualization
    """
    taxon_final_subnetworks = {}
    
    for namespace in taxon_networks.keys():
        network = taxon_networks[namespace]
        if len(network.nodes) == 0:
            continue
        
        network_copy = network.copy()
        
        nodes_to_remove = [node for node, data in network_copy.nodes(data=True) 
                          if data.get('log(Odds Ratio)') <= enrichment_threshold]
        network_copy.remove_nodes_from(nodes_to_remove)
        
        if len(network_copy.nodes) == 0:
            continue
        
        subnetworks = sorted(
            nx.connected_components(network_copy), 
            key=len, 
            reverse=True)
        
        subnetworks = [subnetwork for subnetwork in subnetworks if len(subnetwork) > 1]
        
        if not subnetworks:
            continue
            
        most_informative_subnetworks = subnetworks[:min(number_of_subnetworks, len(subnetworks))]
        
        subgraphs = [network_copy.subgraph(component).copy() for component in most_informative_subnetworks]
        
        if not subgraphs:
            continue

        taxon_final_subnetworks[namespace] = subgraphs
    
    return taxon_name, taxon_final_subnetworks

def get_top_subnetworks(
        networks_dict: Dict[str, Dict[str, nx.Graph]],
        enrichment_threshold: float,
        number_of_subnetworks: int,
        executor: concurrent.futures.ProcessPoolExecutor
        ) -> Dict[str, Dict[str, List[nx.Graph]]]:

    """
    Extract most informative subnetworks from the enrichment networks created based on gene overlap.
    
    For each taxon and namespace combination, this function:
    1. Filters out nodes with log(Odds Ratio) below or equal to the enrichment threshold
    2. Identifies connected components in the filtered network
    3. Selects the largest subnetworks with more than one node
    4. Returns up to the specified number of most informative subnetworks
    """
    
    final_subnetworks = dict()

    futures = []
    
    for taxon_name in networks_dict.keys():
        futures.append(
            executor.submit(
                process_taxon_subnetworks,
                taxon_name,
                networks_dict[taxon_name],
                enrichment_threshold,
                number_of_subnetworks
            )
        )
    
    for future in concurrent.futures.as_completed(futures):
        try:
            taxon_name, taxon_final_subnetworks = future.result()
            if taxon_final_subnetworks:
                final_subnetworks[taxon_name] = taxon_final_subnetworks
        except Exception as e:
            print(f"Error processing subnetworks for taxon: {e}")

    return final_subnetworks

def create_single_network_plot(
        taxon_name: str, 
        taxon_results: pd.DataFrame, 
        namespace: str, 
        taxon_subnetworks: List[nx.Graph], 
        taxonomic_level: str, 
        output_dir: Path) -> Tuple[str, str]:

    """
    Creates a network visualization for a single taxon and namespace combination.
    
    This function:
    1. Processes GO term subnetworks for the specified taxon and namespace
    2. Arranges multiple subnetworks in a grid layout
    3. Creates an interactive visualization where:
       - Nodes represent GO terms
       - Edges represent similarity between terms
       - Node size corresponds to protein count
       - Node color represents enrichment score
    4. Saves the visualization as both HTML and SVG files

    """

    if not taxon_subnetworks:
        return taxon_name, namespace
        
    namespace_results = taxon_results[taxon_results['Namespace'] == namespace]
    
    go_term_data = {
        row['GO Term ID']: row for _, row in namespace_results.iterrows()
    }
    
    max_enrichment = namespace_results['log(Odds Ratio)'].max()
    min_enrichment = namespace_results['log(Odds Ratio)'].min()

    spacing = 2.5
    grid_cols = math.ceil(math.sqrt(len(taxon_subnetworks)))
    
    combined_subnetwork_layout = dict()
    subnetwork_node_groups = list()
    values_to_size = dict()
    node_sizes = dict()
    edge_data = list()

    for i, subnetwork in enumerate(taxon_subnetworks):
        nodes = list(subnetwork.nodes)
        subnetwork_node_groups.append(nodes)    
    
        subnetwork_layout = nx.kamada_kawai_layout(subnetwork)
        
        grid_x = i % grid_cols
        grid_y = i // grid_cols

        for node, pos in subnetwork_layout.items():
            combined_subnetwork_layout[node] = (pos[0] + grid_x * spacing, pos[1] + grid_y * spacing)
            
            if taxonomic_level == "taxonomy":
                values_to_size[node] = go_term_data[node]["N with GO term"]
            elif taxonomic_level == "species":
                values_to_size[node] = go_term_data[node]["N Study with term"]
        
        edge_data.extend([(u, v, data.get('weight')) for u, v, data in subnetwork.edges(data=True)])
    
    max_sizing_value = max(values_to_size.values())
    min_sizing_value = min(values_to_size.values())
    mid_sizing_value = sorted(values_to_size.values())[len(values_to_size.values()) // 2]
    
    for node in values_to_size.keys():
        node_sizes[node] = 7 + (28 * (math.sqrt(values_to_size[node] / math.sqrt(max_sizing_value))))

    edge_weights = [weight for _, _, weight in edge_data]
    min_edge_weight = min(edge_weights)
    max_edge_weight = max(edge_weights)
    mid_edge_weight = (min_edge_weight + max_edge_weight) / 2
    
    fig = go.Figure()

    for u, v, weight in edge_data:
        x0, y0 = combined_subnetwork_layout[u]
        x1, y1 = combined_subnetwork_layout[v]
        
        min_edge_width = 0.5
        max_edge_width = 3.0
        edge_width = min_edge_width + (max_edge_width - min_edge_width) * weight
        
        fig.add_trace(
            go.Scatter(
                x=[x0, x1, None],
                y=[y0, y1, None],
                mode='lines',
                line=dict(width=edge_width, color='rgba(0,0,0,0.08)'),
                showlegend=False,
            )
        )

    for nodes in subnetwork_node_groups:
        node_x = list()
        node_y = list()
        node_text = list()
        node_size_values = list()
        node_colors = list()
        hover_data_list = list()
        
        for node in nodes:
            row = go_term_data[node]
            x, y = combined_subnetwork_layout[node]
            node_x.append(x)
            node_y.append(y)
            node_text.append(node)
            node_size_values.append(node_sizes.get(node))
            node_colors.append(row['log(Odds Ratio)'])
        
            if taxonomic_level == "taxonomy":
                hover_text = "<b>{}</b><br>".format(row['Name']) + \
                "<b>GO ID:</b> {}<br>".format(row['GO Term ID']) + \
                "<b>log(Odds Ratio):</b> {:.2f}<br>".format(row['log(Odds Ratio)']) + \
                "<b>Odds Ratio:</b> {}<br>".format(row['Formatted Odds Ratio']) + \
                "<b>Stat. Sig.:</b> {}<br>".format(row['Formatted Statistical significance']) + \
                "<b>Heterogeneity:</b> {:.2f}".format(row['Heterogeneity'])
            
            elif taxonomic_level == "species":
                hover_text = "<b>{}</b><br>".format(row['Name']) + \
                "<b>GO ID:</b> {}<br>".format(row['GO Term ID']) + \
                "<b>log(Odds Ratio):</b> {:.2f}<br>".format(row['log(Odds Ratio)']) + \
                "<b>Odds Ratio:</b> {}<br>".format(row['Formatted Odds Ratio']) + \
                "<b>Stat. Sig.:</b> {}<br>".format(row['Formatted Statistical significance']) + \
                "<b>Number of Study proteins with term:</b> {}<br>".format(row['N Study with term']) + \
                "<b>Number of Study proteins without term:</b> {}<br>".format(row['N Study without term']) + \
                "<b>Number of Background proteins with term:</b> {}<br>".format(row['N Background with term']) + \
                "<b>Number of Background protein without term:</b> {}".format(row['N Background without term'])

            hover_data_list.append(hover_text)
        
        fig.add_trace(
            go.Scatter(
                x=node_x,
                y=node_y,
                text=node_text,
                mode='markers+text',
                textposition="top center",
                marker=dict(
                    size=node_size_values,
                    color=node_colors,
                    colorscale='matter',
                    opacity=1.0,
                    line=dict(width=1, color='black'),
                    colorbar=dict(
                        title="log(Odds Ratio)",
                        titleside="right",
                        thickness=15,
                        titlefont=dict(size=12),
                        len=0.4,
                        yanchor="middle",
                        y=0.8,
                        x=1
                    ),
                    cmin=min_enrichment,
                    cmax=max_enrichment
                ),
                textfont=dict(
                    family="Arial",
                    color='black',
                    weight='bold',
                    size=9,
                ),
                hoverinfo='text',
                hovertext=hover_data_list,
                showlegend=False
            )
        )
    
    fig.update_layout(
        showlegend=False,
        margin=dict(l=0, r=0, t=40, b=0),
        width=940, 
        height=460,
        plot_bgcolor='white',
        title=f"Enrichment network: {taxon_name} - {namespace}"
    )
    
    fig.update_xaxes(
        showticklabels=False, 
        showgrid=True, 
        zeroline=False,
        gridcolor='rgba(0,0,0,0.05)',
    )
    
    fig.update_yaxes(
        showticklabels=False, 
        showgrid=True, 
        zeroline=False,
        gridcolor='rgba(0,0,0,0.05)',
    )

    min_size = min(node_sizes.values()) if node_sizes else 0
    max_size = max(node_sizes.values()) if node_sizes else 0
    mid_size = min_size + (max_size - min_size) * 0.5
    
    size_to_study_count = {
        min_size: min_sizing_value,
        mid_size: mid_sizing_value,
        max_size: max_sizing_value
    }
    
    for size in [min_size, mid_size, max_size]:
        study_count = size_to_study_count[size]
        fig.add_trace(
            go.Scatter(
                x=[None], 
                y=[None],
                mode='markers',
                marker=dict(
                    size=size, 
                    color='black',
                    line=dict(width=1)
                ),
                name=f"{study_count}",
                legendgroup="node_sizes",
                legendgrouptitle_text="Study Count",
                legendgrouptitle_font=dict(size=12)
            )
        )
    
    min_edge_width = 0.5
    max_edge_width = 3.0
    mid_edge_width = (min_edge_width + max_edge_width) / 2
    
    edge_width_to_jaccard = {
        min_edge_width: min_edge_weight,
        mid_edge_width: mid_edge_weight,
        max_edge_width: max_edge_weight
    }
    
    for width in [min_edge_width, mid_edge_width, max_edge_width]:
        jaccard = edge_width_to_jaccard[width]
        fig.add_trace(
            go.Scatter(
                x=[None], 
                y=[None],
                mode='lines',
                line=dict(width=width, color='black'),
                name=f"{jaccard:.2f}",
                legendgroup="edge_widths",
                legendgrouptitle_text="Jaccard Index",
                legendgrouptitle_font=dict(size=12)
            )
        )
    
    fig.update_layout(
        showlegend=True,
        legend=dict(
            orientation="v",
            yanchor="bottom",
            y=0.2,
            xanchor="center",
            x=1.05,
            itemclick=False,
            itemdoubleclick=False,
        ),
        margin=dict(l=0, r=0, t=40, b=0)
    )
    
    namespace_dir = output_dir / namespace
    namespace_dir.mkdir(parents=True, exist_ok=True)
    
    html_file = f"{taxon_name}_network.html"
    svg_file = f"{taxon_name}_network.svg"

    fig.write_html(namespace_dir / html_file)
    fig.write_image(namespace_dir / svg_file)
    
    return taxon_name, namespace

def enrichment_network_plot(
        top_subnetworks: Dict[str, Dict[str, List[nx.Graph]]],
        results_dict: Dict[str, pd.DataFrame],
        taxonomic_level: str, 
        output_dir: Path,
        executor: concurrent.futures.ProcessPoolExecutor):
    
    """    
    This function creates network visualization plots for Gene Ontology (GO) term
    enrichment analysis results. It arranges multiple subnetworks in a grid layout,
    where nodes represent enriched GO terms and edges represent the Jaccard index between
    them. The Jaccard index is used here to represent the protein overlap between the terms.
    Node size corresponds to the number of proteins associated with the term,
    and node color represents the enrichment score (log odds ratio).
    """

    for taxon_name, taxon_results in results_dict.items():
        taxon_results['Formatted Statistical significance'] = taxon_results['Statistical significance'].apply(lambda x: f"{x:.2e}")
        taxon_results['Odds Ratio'] = np.exp(taxon_results['log(Odds Ratio)'])
        taxon_results['Formatted Odds Ratio'] = taxon_results['Odds Ratio'].apply(lambda x: f"{x:.2f}")

    plot_tasks = []
    
    for taxon_name in top_subnetworks.keys():
        if taxon_name not in results_dict:
            continue
            
        taxon_results = results_dict[taxon_name]
        
        for namespace in top_subnetworks[taxon_name].keys():
            taxon_subnetworks = top_subnetworks[taxon_name][namespace]
            plot_tasks.append((taxon_name, taxon_results, namespace, taxon_subnetworks, taxonomic_level, output_dir))
    

    futures = [
        executor.submit(create_single_network_plot, *task) 
        for task in plot_tasks
    ]
    
    for future in concurrent.futures.as_completed(futures):
        try:
            taxon_name, namespace = future.result()

        except Exception as e:
            print(f"Unexpected error in network plot task: {e}")
        

def wrap_text(
        text: str, 
        width: int) -> str:
    
    """
    Helper function to wrap text labels to fit better in the enrichment bar chart.
    """

    return "<br>".join(textwrap.wrap(text, width))

def create_single_barplot(
        taxon_name: str, 
        taxon_results: pd.DataFrame, 
        namespace: str, 
        taxonomic_level: str, 
        output_dir: Path) -> Tuple[str, str, Optional[str]]:
    
    """
    Creates a bar plot for a single taxon and namespace combination.
    
    This function:
    1. Processes the enrichment results for a specific taxon and namespace
    2. Identifies the top 20 most statistically significant GO terms
    3. Creates a horizontal bar chart with log(Odds Ratio) values
    4. Saves the plot as both HTML and SVG files
    """

    namespace_results = taxon_results[taxon_results['Namespace'] == namespace]
    
    if namespace_results.empty:
        return taxon_name, namespace, "skipped_empty"
        
    fig = go.Figure()

    sorted_namespace_results = namespace_results.sort_values(by="Statistical significance", ascending=True).head(20)
    sorted_namespace_results = sorted_namespace_results.sort_values(by="log(Odds Ratio)", ascending=True)
    
    sorted_rows = sorted_namespace_results.to_dict('records')
    
    if not sorted_rows:
        return taxon_name, namespace, "skipped_no_data"
    
    wrapped_names = [wrap_text(row['Name'], 30) for row in sorted_rows]
    
    hover_data_list = []
    
    for row in sorted_rows:
        if taxonomic_level == "taxonomy":
            hover_text = "<b>{}</b><br>".format(row['Name']) + \
            "<b>GO ID:</b> {}<br>".format(row['GO Term ID']) + \
            "<b>log(Odds Ratio):</b> {:.2f}<br>".format(row['log(Odds Ratio)']) + \
            "<b>Odds Ratio:</b> {}<br>".format(row['Formatted Odds Ratio']) + \
            "<b>Stat. Sig.:</b> {}<br>".format(row['Formatted Statistical significance']) + \
            "<b>-log10(Stat. Sig.):</b> {:.2f}<br>".format(row['-log10(Stat. Sig.)']) + \
            "<b>Heterogeneity:</b> {:.2f}".format(row['Heterogeneity'])
        
        elif taxonomic_level == "species":
            hover_text = "<b>{}</b><br>".format(row['Name']) + \
            "<b>GO ID:</b> {}<br>".format(row['GO Term ID']) + \
            "<b>log(Odds Ratio):</b> {:.2f}<br>".format(row['log(Odds Ratio)']) + \
            "<b>Stat. Sig.:</b> {}<br>".format(row['Formatted Statistical significance']) + \
            "<b>-log10(Stat. Sig.):</b> {:.2f}<br>".format(row['-log10(Stat. Sig.)']) + \
            "<b>Number of Study proteins with term:</b> {}<br>".format(row['N Study with term']) + \
            "<b>Number of Study proteins without term:</b> {}<br>".format(row['N Study without term']) + \
            "<b>Number of Background proteins with term:</b> {}<br>".format(row['N Background with term']) + \
            "<b>Number of Background protein without term:</b> {}".format(row['N Background without term'])
        
        else: 
            print("Incorrect taxonomic level specification")
            return taxon_name, namespace, "error_taxonomic_level"

        hover_data_list.append(hover_text)
    
    color_values = [row['-log10(Stat. Sig.)'] for row in sorted_rows]
    log_odds_values = [row['log(Odds Ratio)'] for row in sorted_rows]
    formatted_log_odds = [f"{val:.2f}" for val in log_odds_values]
    
    fig.add_trace(go.Bar(
        x=log_odds_values,
        y=wrapped_names,
        orientation='h',
        marker=dict(
            color=color_values,
            colorscale='Aggrnyl',
            colorbar=dict(
                title="-log10(Stat. Sig.)",
                titleside="right",
                titlefont=dict(size=12),
                tickfont=dict(size=10),
                len=0.2,
                thickness=15,
                x=1,
                y=0.9,
                yanchor="middle"
            )
        ),
        hovertext=hover_data_list,
        text=formatted_log_odds,
        textposition='auto'
    ))
    
    fig.update_layout(
        width=460,
        height=920,
        margin=dict(l=0, r=0, t=40, b=0),
        showlegend=False,
        dragmode=False,
        plot_bgcolor='white',
        xaxis=dict(
            title='log(Odds Ratio)',
            title_font=dict(size=12),
            tickfont=dict(size=10),
            showgrid=True,
            gridcolor='rgba(0,0,0,0.05)'
        ),
        yaxis=dict(
            title='',
            tickfont=dict(size=12),
            showgrid=True,
            gridcolor='rgba(0,0,0,0.05)'
        ),
        bargap=0.4,
        title=f"Bar chart: {taxon_name} - {namespace}"
    )

    namespace_dir = output_dir / namespace
    namespace_dir.mkdir(parents=True, exist_ok=True)
    
    html_file = f"{taxon_name}_bar.html"
    svg_file = f"{taxon_name}_bar.svg"

    fig.write_html(namespace_dir / html_file)
    fig.write_image(namespace_dir / svg_file)
    
    return taxon_name, namespace
    
def enrichment_bar_plot(
        results_dict: Dict[str, pd.DataFrame], 
        taxonomic_level: str,
        output_dir: Path,
        executor: concurrent.futures.ProcessPoolExecutor) -> go.Figure:
    
    """
    Creates a bar plot for the top 20 most statistically significant GO terms identified by TaxaGO. 
    Y-axis represents GO term names and X-axis represents log(Odds Ratio) values.
    Bars are colored based on -log10(Statistical significance).
    """ 

    for taxon_name, taxon_results in results_dict.items():
            min_value = np.nextafter(0, 1)
            taxon_results['-log10(Stat. Sig.)'] = -np.log10(
                taxon_results['Statistical significance'].replace(0, min_value)
            )
            taxon_results['Formatted Statistical significance'] = taxon_results['Statistical significance'].apply(lambda x: f"{x:.2e}")
            taxon_results['Odds Ratio'] = np.exp(taxon_results['log(Odds Ratio)'])
            taxon_results['Formatted Odds Ratio'] = taxon_results['Odds Ratio'].apply(lambda x: f"{x:.2f}")

    plot_tasks = []
    for taxon_name, taxon_results in results_dict.items():
        for namespace in taxon_results["Namespace"].unique():
            plot_tasks.append((taxon_name, taxon_results, namespace, taxonomic_level, output_dir))
    

    futures = [
        executor.submit(create_single_barplot, *task) 
        for task in plot_tasks
    ]
    
    for future in concurrent.futures.as_completed(futures):
        try:
            taxon_name, namespace = future.result()

        except Exception as e:
            print(f"Unexpected error in bar plot task: {e}")
        
def create_single_bubbleplot(
        taxon_name: str, 
        taxon_results: pd.DataFrame, 
        namespace: str, 
        taxonomic_level: str, 
        output_dir: Path) -> Tuple[str, str, Optional[str]]:

    """
    Creates a bubble plot for a single taxon and namespace combination.
    
    This function:
    1. Filters the enrichment results for a specific taxon and namespace
    2. Creates a bubble plot where:
       - X-axis represents log(Odds Ratio)
       - Y-axis represents -log10(Statistical Significance)
       - Bubble size represents either protein count or species percentage
    3. Labels the top 10 most significant GO terms on the plot
    4. Saves the visualization as both HTML and SVG files
    """

    namespace_results = taxon_results[taxon_results['Namespace'] == namespace]
   
    if namespace_results.empty:
        return taxon_name, namespace
    
    namespace_rows = namespace_results.to_dict('records')
    
    if not namespace_rows:
        return taxon_name, namespace
        
    row_lookup = {row['GO Term ID']: row for row in namespace_rows}
    
    fig = go.Figure()

    top_10_ids = namespace_results.nlargest(10, '-log10(Stat. Sig.)')['GO Term ID'].tolist()
    
    bubble_size_col = "Species Percentage" if taxonomic_level == "taxonomy" else "N Study with term"
    
    hover_data_list = []
    log_odds_values = []
    log_stat_sig_values = []
    sizes = []

    for row in namespace_rows:
        if taxonomic_level == "taxonomy":
            hover_text = "<b>{}</b><br>".format(row['Name']) + \
            "<b>GO ID:</b> {}<br>".format(row['GO Term ID']) + \
            "<b>log(Odds Ratio):</b> {:.2f}<br>".format(row['log(Odds Ratio)']) + \
            "<b>Stat. Sig.:</b> {}<br>".format(row['Formatted Statistical significance']) + \
            "<b>-log10(Stat. Sig.):</b> {:.2f}<br>".format(row['-log10(Stat. Sig.)']) + \
            "<b>Heterogeneity:</b> {:.2f}".format(row['Heterogeneity'])

        elif taxonomic_level == "species":
            hover_text =  "<b>{}</b><br>".format(row['Name']) + \
            "<b>GO ID:</b> {}<br>".format(row['GO Term ID']) + \
            "<b>log(Odds Ratio):</b> {:.2f}<br>".format(row['log(Odds Ratio)']) + \
            "<b>Stat. Sig.:</b> {}<br>".format(row['Formatted Statistical significance']) + \
            "<b>-log10(Stat. Sig.):</b> {:.2f}<br>".format(row['-log10(Stat. Sig.)']) + \
            "<b>Number of Study proteins with term:</b> {}<br>".format(row['N Study with term']) + \
            "<b>Number of Study proteins without term:</b> {}<br>".format(row['N Study without term']) + \
            "<b>Number of Background proteins with term:</b> {}<br>".format(row['N Background with term']) + \
            "<b>Number of Background protein without term:</b> {}".format(row['N Background without term'])
        
        else: 
            return taxon_name, namespace, "error_taxonomic_level"

        hover_data_list.append(hover_text)
        log_odds_values.append(row['log(Odds Ratio)'])
        log_stat_sig_values.append(row['-log10(Stat. Sig.)'])
        sizes.append(row[bubble_size_col])
    
    max_bubble_size = 30
    min_bubble_size = 5

    if len(sizes) > 0 and max(sizes) > 0:
        normalized_sizes = [
            min_bubble_size + (max_bubble_size - min_bubble_size) * (size / max(sizes)) 
            for size in sizes
        ]
    else:
        normalized_sizes = [min_bubble_size] * len(sizes)
    
    fig.add_trace(go.Scatter(
        x=log_odds_values,
        y=log_stat_sig_values,
        mode='markers',
        marker=dict(
            size=normalized_sizes,
            color="#FAB472",
            line=dict(width=1, color='#333333'),
            opacity=1
        ),
        hovertext=hover_data_list,
        showlegend=False
    ))
    
    fig.update_layout(
        width=940,
        height=460,
        margin=dict(l=0, r=0, t=40, b=0),
        xaxis=dict(
            title='log(Odds Ratio)',
            title_font=dict(size=12),
            tickfont=dict(size=10),
            showgrid=True,
            gridcolor='rgba(0,0,0,0.05)'
        ),
        yaxis=dict(
            title='-log10(Statistical Significance)',
            title_font=dict(size=12),
            tickfont=dict(size=10),
            showgrid=True,
            gridcolor='rgba(0,0,0,0.05)'
        ),
        plot_bgcolor='white',
        title=f"Bubble plot: {taxon_name} - {namespace}"
    )
    
    for term_id in top_10_ids:
        if term_id in row_lookup:
            row = row_lookup[term_id]
            fig.add_annotation(
                x=row['log(Odds Ratio)'],
                y=row['-log10(Stat. Sig.)'],
                text=row['GO Term ID'],
                font=dict(
                    size=8,
                    weight="bold",
                    family="Arial"
                ),
                showarrow=False,
                xanchor="left",
                yanchor="auto"
            )
    
    namespace_dir = output_dir / namespace
    namespace_dir.mkdir(parents=True, exist_ok=True)
    
    html_file = f"{taxon_name}_bubble.html"
    svg_file = f"{taxon_name}_bubble.svg"
    fig.write_html(namespace_dir / html_file)
    fig.write_image(namespace_dir / svg_file)
    
    return taxon_name, namespace

def enrichment_bubble_plot(
        results_dict: Dict[str, pd.DataFrame], 
        taxonomic_level: str,
        output_dir: Path,
        executor: concurrent.futures.ProcessPoolExecutor) -> go.Figure:
    
    """
    Creates a bubble plot for the enriched GO terms identified by TaxaGO. 
    Y-axis represents -log10(Statistical significane) and X-axis represents log(Odds Ratio) values.
    Bubble size is specified by number of proteins associated with the GO term in the study population for the species level plots or the percentage of species inside the taxonomy that have this GO term enriched.
    """

    for taxon_name, taxon_results in results_dict.items():
        min_value = np.nextafter(0, 1)
        taxon_results['-log10(Stat. Sig.)'] = -np.log10(
            taxon_results['Statistical significance'].replace(0, min_value)
        )
        taxon_results['Formatted Statistical significance'] = taxon_results['Statistical significance'].apply(lambda x: f"{x:.2e}")
        taxon_results['Odds Ratio'] = np.exp(taxon_results['log(Odds Ratio)'])
        taxon_results['Formatted Odds Ratio'] = taxon_results['Odds Ratio'].apply(lambda x: f"{x:.2f}")
    
    plot_tasks = []
    for taxon_name, taxon_results in results_dict.items():
        for namespace in taxon_results["Namespace"].unique():
            plot_tasks.append((taxon_name, taxon_results, namespace, taxonomic_level, output_dir))
    

    futures = [
        executor.submit(create_single_bubbleplot, *task) 
        for task in plot_tasks
    ]
    
    for future in concurrent.futures.as_completed(futures):
        try:
            taxon_name, namespace = future.result()
        except Exception as e:
            print(f"Unexpected error in bubble plot task: {e}")

def main():
    parser = argparse.ArgumentParser(description='TaxaGO enrichment results visualization script')
    parser.add_argument('-s', 
                        dest='study_data', 
                        required=True, 
                        help='Path where study data are stored either in FASTA format of in CSV file.')
    parser.add_argument('-r', 
                        dest='results_dir', 
                        required=True, 
                        help='Directory containing analysis results.')
    args = parser.parse_args()
    
    results_dir = Path(args.results_dir)
    
    if (results_dir / "combined_taxonomy_results").exists():
        results_dir = results_dir / "combined_taxonomy_results"
        taxonomic_level = "taxonomy"
        taxonomies = list()
        for file in results_dir.glob("*_GOEA_results.txt"):
            taxonomy = file.name.strip().split("_GOEA_results.txt")[0]
            taxonomies.append(taxonomy)
    elif (results_dir / "single_taxon_results").exists():
        results_dir = results_dir / "single_taxon_results"
        taxonomic_level = "species"
    else:
        print("No results directories found.")
        exit()
    
    plots_dir = results_dir / "plots"

    if plots_dir.exists():
        for item in plots_dir.iterdir():
            if item.is_file():
                item.unlink()
            elif item.is_dir():
                shutil.rmtree(item)
    else:
        plots_dir.mkdir(parents=True)
    
    for namespace in ['Biological Process', 'Molecular Function', 'Cellular Component']:
        namespace_dir = plots_dir / namespace
        namespace_dir.mkdir(parents=True)

    results_dict = load_results(
        results_dir)
    lineage = create_lineage(
        get_default_lineage())

    if taxonomic_level == 'taxonomy':
        taxonomy_id_map = get_taxon_ids(
            lineage, 
            taxonomies)
    
    # THIS IS A BIG TIME BOTTLENECK CURRENTLY
    network_data = create_network_data(
            args.study_data,
            get_default_background(),
            lineage,
            taxonomic_level,
            taxonomy_id_map,
            results_dict)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        networks = build_go_term_networks(
            network_data, 
            results_dict,
            0.5,
            executor)

        top_networks = get_top_subnetworks(
            networks,
            1.5,
            4,
            executor)

        enrichment_bar_plot(
            results_dict,
            taxonomic_level,
            plots_dir,
            executor)
    
        enrichment_bubble_plot(
            results_dict,
            taxonomic_level,
            plots_dir,
            executor)

        enrichment_network_plot(
            top_networks,
            results_dict,
            taxonomic_level,
            plots_dir,
            executor)

if __name__ == "__main__":
    main()