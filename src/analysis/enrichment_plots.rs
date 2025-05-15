use rustc_hash::{FxHashMap, FxHashSet};
use std::path::PathBuf;
use std::error::Error;
use std::fs;
use plotly::{
    Plot, Bar, Layout, Scatter,
    common::{
        Title, Font, HoverInfo,
        ColorScale, ColorScalePalette,
        Marker, ColorBar, Anchor, Side,
        ThicknessMode, Orientation, Mode,
    },
    layout::{
        Axis, Margin,
        DragMode, RangeMode,
    },
    color::NamedColor
};
use textwrap::wrap;
use std::cmp::Ordering::Equal;
use crate::{
    parsers::{
        background_parser::*,
        study_parser::*,
        obo_parser::*
    },
    analysis::{
        enrichment_analysis::*,
        phylogenetic_meta_analysis::*
    }
};
use petgraph::{
    Undirected,
    graph::NodeIndex,
    stable_graph::StableGraph,
    visit::EdgeRef
}; 
use rayon::prelude::*;
use itertools::Itertools;
use strum::IntoEnumIterator;
use std::collections::VecDeque;


#[derive(Debug, Copy, Clone)]
pub struct NetworkNode {
    pub go_id: GOTermID,
    pub lor: f64,
    pub stat_sig: f64,
    // pub size_statistic: u32
}
pub type JaccardIndex = f64;
pub type GoTermNetworkGraph = StableGraph<NetworkNode, JaccardIndex, Undirected>;

#[derive(Debug, Clone)]
pub struct TermPlotData {
    pub go_id: GOTermID,
    pub name: String,
    pub wrapped_name: String,
    pub lor: f64,
    pub stat_sig: f64,
    pub minus_log10_p_value: f64,
    pub namespace: NameSpace,
    pub hover_text: String,
}

#[derive(Debug, Clone, Default)]
pub struct PlotData {
    pub terms: Vec<TermPlotData>, 
}

pub trait EnrichmentResult {
    fn log_odds_ratio(&self) -> f64;
    fn p_value(&self) -> f64;
}

impl EnrichmentResult for GOTermResults {
    fn log_odds_ratio(&self) -> f64 {
        self.log_odds_ratio
    }

    fn p_value(&self) -> f64 {
        self.p_value
    }
}

impl EnrichmentResult for TaxonomyGOResult {
    fn log_odds_ratio(&self) -> f64 {
        self.log_odds_ratio
    }

    fn p_value(&self) -> f64 {
        self.p_value
    }
}

#[derive(Clone)]
pub enum ProteinDataProvider<'a> {
    Species(
        &'a FxHashMap<String, GOTermToProteinSet>
    ), 
    Taxonomy {
        species_data_by_id: &'a FxHashMap<TaxonID, GOTermToProteinSet>,
        taxonomy_to_species_ids: &'a FxHashMap<String, Vec<TaxonID>>,
    },
}

impl<'a> ProteinDataProvider<'a> {
    fn get_proteins_for_taxon(
        &self,
        taxon_name: &str,
        relevant_go_ids: &FxHashSet<GOTermID>,
    ) -> GOTermToProteinSet {
        match self {
            ProteinDataProvider::Species(species_map) => {
                species_map.get(taxon_name).cloned().unwrap_or_default()
            }
            ProteinDataProvider::Taxonomy {
                species_data_by_id,
                taxonomy_to_species_ids,
            } => {
                let species_ids = taxonomy_to_species_ids.get(taxon_name).unwrap();
                let mut aggregated_go_to_proteins: GOTermToProteinSet = FxHashMap::default();
                for go_id in relevant_go_ids {
                    let mut proteins_for_go: FxHashSet<Protein> = FxHashSet::default();
                    for species_id in species_ids {
                        if let Some(species_go_map) = species_data_by_id.get(species_id) {
                            if let Some(proteins) = species_go_map.get(go_id) {
                                proteins_for_go.extend(proteins.iter().cloned());
                            }
                        }
                    }
                    if !proteins_for_go.is_empty() {
                        aggregated_go_to_proteins.insert(*go_id, proteins_for_go);
                    }
                }
                aggregated_go_to_proteins
            }
        }
    }
}


fn wrap_text(
    text: &str, 
    width: usize
) -> String {
    wrap(text, width).join("<br>")
}

pub fn process_species_data(
    mut significant_species_results: FxHashMap<TaxonID, FxHashMap<GOTermID, GOTermResults>>,
    study_pop: &StudyPop,
    taxon_id_to_name: &FxHashMap<TaxonID, String>,
) -> (
    FxHashMap<String, FxHashMap<GOTermID, GOTermResults>>,
    FxHashMap<String, FxHashMap<GOTermID, FxHashSet<Protein>>>,
) {

    let capacity = significant_species_results.len();
    let mut significant_results_by_name: FxHashMap<String, FxHashMap<GOTermID, GOTermResults>> =
        FxHashMap::with_capacity_and_hasher(capacity, Default::default());

    let original_go_term_to_protein_set = &study_pop.go_term_to_protein_set;
    let mut new_go_term_to_protein_set: FxHashMap<String, FxHashMap<GOTermID, FxHashSet<Protein>>> =
        FxHashMap::with_capacity_and_hasher(original_go_term_to_protein_set.len(), Default::default());

    for (taxon_u32_id, go_term_map_value) in significant_species_results.drain() {
        let taxon_name = taxon_id_to_name.get(&taxon_u32_id).unwrap();
        significant_results_by_name.insert(taxon_name.clone(), go_term_map_value);
    }

    for (taxon_u32_id, go_term_data_for_taxon) in original_go_term_to_protein_set.iter() {
        let taxon_name = taxon_id_to_name.get(taxon_u32_id).unwrap();
        new_go_term_to_protein_set.insert(taxon_name.clone(), go_term_data_for_taxon.clone());
    }

    (significant_results_by_name, new_go_term_to_protein_set)
}

pub fn prepare_plot_data<R>(
    significant_results: &FxHashMap<String, FxHashMap<GOTermID, R>>,
    ontology: &OboMap
) -> FxHashMap<String, FxHashMap<NameSpace, PlotData>>
where
    R: EnrichmentResult + Clone + Send + Sync
{
    significant_results
        .par_iter()
        .filter_map(|(species_name, go_term_results_map)| {
            let mut species_terms_by_namespace: FxHashMap<NameSpace, Vec<TermPlotData>> =
                FxHashMap::default();

            for (go_id, results) in go_term_results_map {
                if let Some(obo_term) = ontology.get(go_id) {
                    let current_p_value = results.p_value();
                    let current_lor = results.log_odds_ratio();
                    let minus_log_10_p = -current_p_value.log10();
                    let go_id_string = format!("GO:{:07}", go_id);
                    let original_name = obo_term.name.clone();
                    let wrapped_display_name = wrap_text(&original_name, 30);

                    let hover_html_content = format!(
                        "<b>Term Name:</b> {}<br><b>Term ID:</b> {}<br><b>log(Odds Ratio):</b> {:.3}<br><b>-log10(Stat. Sig.):</b> {:.3}",
                        original_name,
                        go_id_string,
                        current_lor,
                        minus_log_10_p,
                    );

                    let rich_term = TermPlotData {
                        go_id: *go_id,
                        name: original_name,
                        wrapped_name: wrapped_display_name,
                        lor: current_lor,
                        stat_sig: current_p_value,
                        minus_log10_p_value: minus_log_10_p,
                        namespace: obo_term.namespace.clone(),
                        hover_text: hover_html_content,
                    };

                    species_terms_by_namespace
                        .entry(obo_term.namespace.clone())
                        .or_default()
                        .push(rich_term);
                }
            }

            let data_for_plotting: FxHashMap<NameSpace, PlotData> =
                species_terms_by_namespace
                    .into_iter()
                    .map(|(ns_enum, mut terms_in_ns_group)| { 
                        terms_in_ns_group.sort_by(|a, b| {
                            a.stat_sig.partial_cmp(&b.stat_sig).unwrap_or(Equal)
                        });

                        (
                            ns_enum,
                            PlotData {
                                terms: terms_in_ns_group,
                            },
                        )
                    })
                    .collect();
            
            if data_for_plotting.is_empty() {
                None
            } else {
                Some((species_name.clone(), data_for_plotting))
            }
        })
        .collect()
}

pub fn bar_plot(
    plot_data_map: &FxHashMap<String, FxHashMap<NameSpace, PlotData>>,
    plots_dir: &PathBuf,
) -> Result<(), Box<dyn Error + Send + Sync>> {

    plot_data_map
        .into_iter()
        .flat_map(|(taxon_name, namespace_map)| {
            namespace_map
                .into_iter()
                .map(move |(namespace, current_plot_data)| {
                    (taxon_name.clone(), namespace, current_plot_data)
                })
        })
        .par_bridge()
        .try_for_each(|(taxon_name, namespace, namespace_plot_data)| -> Result<(), Box<dyn Error + Send + Sync>> {

            let namespace_str: String = match namespace {
                NameSpace::BiologicalProcess => "Biological_Process".to_string(),
                NameSpace::MolecularFunction => "Molecular_Function".to_string(),
                NameSpace::CellularComponent => "Cellular_Component".to_string(),
            };
            let namespace_subdir = plots_dir.join(&namespace_str);
            fs::create_dir_all(&namespace_subdir)?;
            
            let mut top_20_terms: Vec<TermPlotData> = namespace_plot_data.terms
                .iter()
                .take(20)
                .cloned()
                .collect();
            
            top_20_terms.sort_by(|a, b| {
                a.lor.partial_cmp(&b.lor).unwrap_or(Equal)
            });
            
            let capacity = top_20_terms.len();
            let mut term_names_display: Vec<String> = Vec::with_capacity(capacity);
            let mut log_odds_ratios_values: Vec<f64> = Vec::with_capacity(capacity);
            let mut minus_log_10_stat_sigs_for_color: Vec<f64> = Vec::with_capacity(capacity);
            let mut html_hover_texts_vec: Vec<String> = Vec::with_capacity(capacity);

            for term_data in top_20_terms {
                term_names_display.push(term_data.wrapped_name.clone());
                log_odds_ratios_values.push(term_data.lor);
                minus_log_10_stat_sigs_for_color.push(term_data.minus_log10_p_value);
                html_hover_texts_vec.push(term_data.hover_text.clone());
            }

            let color_bar = ColorBar::new()
                .title(
                    Title::from("-log10(Stat. Sig.)")
                        .side(Side::Right)
                        .font(Font::new().size(12)),
                )
                .tick_font(Font::new().size(10))
                .len_mode(ThicknessMode::Pixels)
                .len(200)
                .thickness(15)
                .x(1.0)
                .y(0.9)
                .y_anchor(Anchor::Middle);

            let marker = Marker::new()
                .color_array(minus_log_10_stat_sigs_for_color)
                .color_scale(ColorScale::Palette(ColorScalePalette::Cividis))
                .color_bar(color_bar)
                .show_scale(true);

            let bar_trace = Bar::new(log_odds_ratios_values, term_names_display)
                .orientation(Orientation::Horizontal)
                .marker(marker)
                .hover_text_array(html_hover_texts_vec)
                .hover_info(HoverInfo::Text)
                .show_legend(false);

            let mut plot = Plot::new();
            plot.add_trace(bar_trace);

            let layout = Layout::new()
                .width(460)
                .height(920)
                .margin(Margin::new()
                    .left(30)
                    .right(0)
                    .top(30)
                    .bottom(0))
                .x_axis(
                    Axis::new()
                        .title(Title::with_text("log(Odds Ratio)").font(Font::new().size(12)))
                        .tick_font(Font::new().size(10))
                        .show_line(true)
                        .line_color(NamedColor::Black)
                        .show_grid(true)
                        .grid_color("rgba(0,0,0,0.05)")
                        .show_tick_labels(true)
                        .auto_margin(true),
                )
                .y_axis(
                    Axis::new() 
                        .title(Title::with_text(""))
                        .tick_font(Font::new().size(12))
                        .show_line(true)
                        .line_color(NamedColor::Black)
                        .show_grid(true)
                        .grid_color("rgba(0,0,0,0.05)")
                        .show_tick_labels(true)
                        .auto_margin(true),
                )
                .drag_mode(DragMode::False)
                .bar_gap(0.4);
            plot.set_layout(layout);

            let plot_name_path = namespace_subdir.join(format!("{}_barplot.html", taxon_name));
            plot.write_html(plot_name_path); 
                                             
            Ok(())
        })?; 

    Ok(())
}

pub fn bubble_plot(
    plot_data_map: FxHashMap<String, FxHashMap<NameSpace, PlotData>>,
    plots_dir: &PathBuf,
) -> Result<(), Box<dyn Error + Send + Sync>> {

    plot_data_map
        .into_iter()
        .flat_map(|(taxon_name, namespace_map)| {
            namespace_map
                .into_iter()
                .map(move |(namespace, current_plot_data)| {
                    (taxon_name.clone(), namespace, current_plot_data)
                })
        })
        .par_bridge()
        .try_for_each(|(taxon_name, namespace, namespace_plot_data)| -> Result<(), Box<dyn Error + Send + Sync>> {

            let namespace_str: String = match namespace {
                NameSpace::BiologicalProcess => "Biological_Process".to_string(),
                NameSpace::MolecularFunction => "Molecular_Function".to_string(),
                NameSpace::CellularComponent => "Cellular_Component".to_string(),
            };
            let namespace_subdir = plots_dir.join(&namespace_str);
            fs::create_dir_all(&namespace_subdir)?;

            let enrichment_values: Vec<f64> = namespace_plot_data.terms.iter().map(|t| t.lor).collect();
            let stat_sig_values: Vec<f64> = namespace_plot_data.terms.iter().map(|t| t.minus_log10_p_value).collect();
            let hover_texts: Vec<String> = namespace_plot_data.terms.iter().map(|t| t.hover_text.clone()).collect();
        
            let scatter_trace = Scatter::new(enrichment_values, stat_sig_values)
                .mode(Mode::Markers)
                .marker(Marker::new().color(NamedColor::DimGray))
                .hover_text_array(hover_texts) 
                .hover_info(HoverInfo::Text) 
                .show_legend(false);

            let mut plot = Plot::new();
            plot.add_trace(scatter_trace);

            let layout = Layout::new()
                .width(940)
                .height(460)
                .margin(Margin::new()
                    .left(60)
                    .right(0)
                    .top(30)
                    .bottom(40))
                .x_axis(
                    Axis::new()
                        .title(Title::with_text("log(Odds Ratio)").font(Font::new().size(12)))
                        .tick_font(Font::new().size(10))
                        .show_line(true)
                        .line_color(NamedColor::Black)
                        .show_grid(true)
                        .grid_color("rgba(0,0,0,0.05)")
                        .show_tick_labels(true)
                        .auto_margin(true)
                )
                .y_axis(
                    Axis::new() 
                        .title(Title::with_text("-log10(Stat. Sig.)").font(Font::new().size(12)))
                        .tick_font(Font::new().size(10))
                        .show_line(true)
                        .line_color(NamedColor::Black)
                        .show_grid(true)
                        .grid_color("rgba(0,0,0,0.05)")
                        .show_tick_labels(true)
                        .auto_margin(true)
                        .range_mode(RangeMode::ToZero)
                );
            plot.set_layout(layout);

            let plot_name_path = namespace_subdir.join(format!("{}_bubble.html", taxon_name));
            plot.write_html(plot_name_path); 
                                             
            Ok(())
        })?; 

    Ok(())
}

pub fn prepare_network_data<R>(
    significant_results: &FxHashMap<String, FxHashMap<GOTermID, R>>,
    protein_provider: &ProteinDataProvider,
    ontology: &OboMap,
) -> FxHashMap<String, FxHashMap<NameSpace, GOTermToProteinSet>>
where
    R: EnrichmentResult + Clone + Send + Sync, 
{
    significant_results
        .par_iter() 
        .map(|(taxon_name, enriched_go_terms_map)| {

            let relevant_go_ids_for_taxon: FxHashSet<GOTermID> = enriched_go_terms_map.keys().cloned().collect();
            let current_taxon_go_to_proteins = protein_provider
                .get_proteins_for_taxon(taxon_name, &relevant_go_ids_for_taxon);

            let network_data_by_namespace = enriched_go_terms_map
                .iter()
                .filter_map(|(go_id, _result)| { 
                    ontology.get(go_id).and_then(|obo_term| { 
                        let namespace = obo_term.namespace.clone();
                        current_taxon_go_to_proteins.get(go_id).map(|protein_set_for_go_term| {
                            (*go_id, namespace, protein_set_for_go_term.clone())
                        })
                    })
                })
                .fold(
                    FxHashMap::<NameSpace, GOTermToProteinSet>::default(),
                    |mut acc, (go_id, namespace, protein_set)| {
                        acc.entry(namespace)
                           .or_insert_with(FxHashMap::default)
                           .insert(go_id, protein_set);
                        acc
                    },
                );
            (taxon_name.clone(), network_data_by_namespace)
        })
        .collect()
}


pub fn build_networks<R>(
    network_data: &FxHashMap<String, FxHashMap<NameSpace, GOTermToProteinSet>>,
    enrichment_results: &FxHashMap<String, FxHashMap<GOTermID, R>>,
    k_communities: usize,
) -> FxHashMap<String, FxHashMap<NameSpace, Vec<GoTermNetworkGraph>>> 
where 
    R: EnrichmentResult + Clone + Send + Sync
{
    network_data
        .par_iter()
        .filter_map(|(taxon_name, taxon_specific_network_data)| {
            enrichment_results
                .get(taxon_name)
                .map(|taxon_specific_enrichment_results| {
                    let mut taxon_networks_graphs: FxHashMap<NameSpace, Vec<GoTermNetworkGraph>> =
                        FxHashMap::default();

                    for current_namespace in NameSpace::iter() {
                        let go_term_proteins_in_namespace = match taxon_specific_network_data
                            .get(&current_namespace)
                        {
                            Some(data) => data,
                            None => {
                                taxon_networks_graphs
                                    .insert(current_namespace.clone(), Vec::new());
                                continue;
                            }
                        };

                        let mut current_term_network: GoTermNetworkGraph = StableGraph::default();
                        let mut term_to_node_index_map: FxHashMap<GOTermID, NodeIndex> = FxHashMap::default();
                        let mut term_to_proteins_map_for_nodes: FxHashMap<GOTermID, &FxHashSet<Protein>> = FxHashMap::default();

                        for (go_term_id, protein_set) in go_term_proteins_in_namespace {
                            let enrichment_detail = taxon_specific_enrichment_results.get(go_term_id).unwrap();
                            let node_data = NetworkNode {
                                go_id: *go_term_id,
                                lor: enrichment_detail.log_odds_ratio(),
                                stat_sig: enrichment_detail.p_value(),
                            };

                            let node_idx = current_term_network.add_node(node_data);
                            term_to_node_index_map.insert(*go_term_id, node_idx);
                            term_to_proteins_map_for_nodes.insert(*go_term_id, protein_set);
                        }

                        let mut protein_to_terms_map: FxHashMap<&Protein, FxHashSet<GOTermID>> = FxHashMap::default();
                        let mut term_node_sizes: FxHashMap<GOTermID, usize> = FxHashMap::default();

                        for (term_id, protein_set) in &term_to_proteins_map_for_nodes {
                            term_node_sizes.insert(*term_id, protein_set.len());
                            for protein in protein_set.iter() {
                                protein_to_terms_map
                                    .entry(protein)
                                    .or_insert_with(FxHashSet::default)
                                    .insert(*term_id);
                            }
                        }

                        let mut candidate_go_pairs: FxHashSet<(GOTermID, GOTermID)> = FxHashSet::default();
                        for go_terms_sharing_protein in protein_to_terms_map.values() {
                            if go_terms_sharing_protein.len() >= 2 {
                                for combo in go_terms_sharing_protein.iter().combinations(2) {
                                    let term1 = *combo[0];
                                    let term2 = *combo[1];
                                    let pair = if term1 < term2 {
                                        (term1, term2)
                                    } else {
                                        (term2, term1)
                                    };
                                    candidate_go_pairs.insert(pair);
                                }
                            }
                        }

                        for (term1_id, term2_id) in candidate_go_pairs {
                            let size1 = *term_node_sizes.get(&term1_id).unwrap_or(&0);
                            let size2 = *term_node_sizes.get(&term2_id).unwrap_or(&0);

                            let proteins1 = match term_to_proteins_map_for_nodes.get(&term1_id) {
                                Some(p_set) => p_set,
                                None => continue, 
                            };
                            let proteins2 = match term_to_proteins_map_for_nodes.get(&term2_id) {
                                Some(p_set) => p_set,
                                None => continue,
                            };

                            let intersection_size = proteins1.intersection(proteins2).count();
                            let union_size = size1 + size2 - intersection_size;

                            let jaccard_similarity: JaccardIndex =
                                (intersection_size as f64) / (union_size as f64);

                            if jaccard_similarity >= 0.25 {
                                if let (Some(&node_idx1), Some(&node_idx2)) = (
                                    term_to_node_index_map.get(&term1_id),
                                    term_to_node_index_map.get(&term2_id),
                                ) {
                                    current_term_network.add_edge(node_idx1, node_idx2, jaccard_similarity);
                                }
                            }
                        }
                        
                        let top_k_subgraphs = extract_top_k_communities(&current_term_network, k_communities);
                        taxon_networks_graphs.insert(current_namespace.clone(), top_k_subgraphs);
                    }
                    (taxon_name.clone(), taxon_networks_graphs)
                })
        })
        .collect()
}

fn get_all_connected_components(
    graph: &GoTermNetworkGraph,
) -> Vec<Vec<NodeIndex>> {
    let mut visited_nodes = FxHashSet::default();
    let mut all_components = Vec::new();

    for node_idx in graph.node_indices() {
        if !visited_nodes.contains(&node_idx) {
            let mut current_component_nodes = Vec::new();
            let mut queue = VecDeque::new();

            queue.push_back(node_idx);
            visited_nodes.insert(node_idx);

            while let Some(u_idx) = queue.pop_front() {
                current_component_nodes.push(u_idx);
                for neighbor_idx in graph.neighbors(u_idx) {
                    if !visited_nodes.contains(&neighbor_idx) {
                        visited_nodes.insert(neighbor_idx);
                        queue.push_back(neighbor_idx);
                    }
                }
            }
            if !current_component_nodes.is_empty() {
                all_components.push(current_component_nodes);
            }
        }
    }
    all_components
}

fn extract_top_k_communities(
    graph: &GoTermNetworkGraph,
    k: usize,
) -> Vec<GoTermNetworkGraph> {
    let mut components_node_indices = get_all_connected_components(graph);

    components_node_indices.sort_by_key(|comp| std::cmp::Reverse(comp.len()));

    let mut top_k_graphs = Vec::new();

    for component_nodes in components_node_indices.into_iter().take(k) {
        if component_nodes.is_empty() {
            continue;
        }

        let mut subgraph: GoTermNetworkGraph = StableGraph::default();
        let mut old_to_new_node_map = FxHashMap::default();

        for &old_node_idx in &component_nodes {
            if let Some(node_weight) = graph.node_weight(old_node_idx) {
                let new_node_idx = subgraph.add_node(node_weight.clone());
                old_to_new_node_map.insert(old_node_idx, new_node_idx);
            }
        }

        let component_node_set: FxHashSet<NodeIndex> = component_nodes.iter().cloned().collect();

        for &old_u_idx in &component_nodes {
            let new_u_idx = old_to_new_node_map[&old_u_idx];

            for edge_ref in graph.edges(old_u_idx) {
                let old_v_idx = edge_ref.target();
                if component_node_set.contains(&old_v_idx) {
                    if old_u_idx <= old_v_idx {
                        if let Some(&new_v_idx) = old_to_new_node_map.get(&old_v_idx) {
                            if let Some(edge_weight) = graph.edge_weight(edge_ref.id()) {
                                subgraph.add_edge(new_u_idx, new_v_idx, *edge_weight);
                            }
                        }
                    }
                }
            }
        }
        top_k_graphs.push(subgraph);
    }

    top_k_graphs
}

