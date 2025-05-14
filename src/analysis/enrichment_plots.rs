use rustc_hash::{FxHashMap, FxHashSet};
use std::path::PathBuf;
use std::error::Error;
use std::fs;
use plotly::{
    Plot, Bar, Layout, Scatter,
    common::{
        Title, Font,
        ColorScale, ColorScalePalette,
        Marker, ColorBar, Anchor, Side,
        ThicknessMode, Orientation, Mode
    },
    layout::{
        Axis, Margin,
        DragMode
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
    stable_graph::StableGraph
}; 
use rayon::prelude::*;
use itertools::Itertools;
use strum::IntoEnumIterator;


#[derive(Debug, Clone)]
pub struct NetworkNode {
    pub go_id: GOTermID,
    pub lor: f64,
    pub stat_sig: f64,
    // pub size_statistic: u32
}
pub type JaccardIndex = f64;
pub type GoTermNetworkGraph = StableGraph<NetworkNode, JaccardIndex, Undirected>;

#[derive(Debug, Clone)]
pub struct TermToPlot {
    pub go_id: GOTermID,
    pub name: String,
    pub lor: f64,
    pub stat_sig: f64,
    pub namespace: NameSpace
}
#[derive(Debug, Clone, Default)]
pub struct PlotData {
    term_names: Vec<String>,
    lor_vec: Vec<f64>,
    stat_sig_vec: Vec<f64>
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

fn wrap_text(
    text: &str, 
    width: usize
) -> String {
    wrap(text, width).join("<br>")
}

pub fn prepare_plot_data<R>(
    significant_results: &FxHashMap<TaxonID, FxHashMap<GOTermID, R>>,
    ontology: &OboMap,
    taxid_species_map: &FxHashMap<TaxonID, String>,

) -> FxHashMap<String, FxHashMap<NameSpace, PlotData>> 
where
    R: EnrichmentResult + Clone + Send + Sync
{
    significant_results
        .par_iter()
        .filter_map(|(taxon_id, go_term_results_map)| {
            let species_name = taxid_species_map.get(taxon_id).cloned().unwrap();

            let mut species_terms_by_namespace: FxHashMap<NameSpace, Vec<TermToPlot>> =
                FxHashMap::default();

            for (go_id, results) in go_term_results_map {
                if let Some(obo_term) = ontology.get(go_id) {
                    let namespace_clone = obo_term.namespace.clone();
                    let term_data = TermToPlot {
                        go_id: *go_id,
                        name: obo_term.name.clone(),
                        lor: results.log_odds_ratio(),
                        stat_sig: results.p_value(),
                        namespace: namespace_clone.clone(),
                    };

                    species_terms_by_namespace
                        .entry(namespace_clone)
                        .or_default()
                        .push(term_data);
                }
            }

            let data_for_plotting: FxHashMap<NameSpace, PlotData> =
                species_terms_by_namespace
                    .into_iter()
                    .map(|(ns_enum, mut terms_in_ns_group)| {
                        terms_in_ns_group.sort_by(|a, b| {
                            a.stat_sig.partial_cmp(&b.stat_sig).unwrap_or(Equal)
                        });
                        
                        let mut top_terms: Vec<TermToPlot> =
                            terms_in_ns_group.into_iter().take(20).collect();
                        
                        top_terms.sort_by(|a, b| {
                            a.lor.partial_cmp(&b.lor).unwrap_or(Equal)
                        });

                        let capacity = top_terms.len();
                        let mut term_names_display: Vec<String> = Vec::with_capacity(capacity);
                        let mut log_odds_ratios_values: Vec<f64> = Vec::with_capacity(capacity);
                        let mut minus_log_10_stat_sigs: Vec<f64> = Vec::with_capacity(capacity);

                        for term_data in top_terms { 
                            term_names_display.push(
                                wrap_text(&term_data.name, 30)
                            );
                            log_odds_ratios_values.push(term_data.lor);
                            minus_log_10_stat_sigs.push(-term_data.stat_sig.log10());
                        }

                        (
                            ns_enum,
                            PlotData {
                                term_names: term_names_display,
                                lor_vec: log_odds_ratios_values,
                                stat_sig_vec: minus_log_10_stat_sigs,
                            },
                        )
                    })
                    .collect();
            
            Some((species_name, data_for_plotting))
            
        })
        .collect()
}

pub fn barplot(
    plot_data: FxHashMap<String, FxHashMap<NameSpace, PlotData>>,
    plots_dir: &PathBuf,
) -> Result<(), Box<dyn Error + Send + Sync>> {

    plot_data
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

            let enrichment_values = namespace_plot_data.lor_vec;
            let go_term_names = namespace_plot_data.term_names;

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
                .color_array(namespace_plot_data.stat_sig_vec)
                .color_scale(ColorScale::Palette(ColorScalePalette::Cividis))
                .color_bar(color_bar)
                .show_scale(true);

            let bar_trace = Bar::new(enrichment_values, go_term_names)
                .orientation(Orientation::Horizontal)
                .marker(marker)
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
    plot_data: FxHashMap<String, FxHashMap<NameSpace, PlotData>>,
    plots_dir: &PathBuf,
) -> Result<(), Box<dyn Error + Send + Sync>> {

    plot_data
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

            let stat_sig = namespace_plot_data.stat_sig_vec;
            let enrichment_values = namespace_plot_data.lor_vec;
        
            let scatter_trace = Scatter::new(enrichment_values, stat_sig)
                .mode(Mode::Markers)
                .marker(Marker::new().color(NamedColor::DimGray))
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
                        .auto_margin(true),
                )
                .y_axis(
                    Axis::new() 
                        .title(Title::with_text("-log10(Stat. Sig.").font(Font::new().size(12)))
                        .tick_font(Font::new().size(10))
                        .show_line(true)
                        .line_color(NamedColor::Black)
                        .show_grid(true)
                        .grid_color("rgba(0,0,0,0.05)")
                        .show_tick_labels(true)
                        .auto_margin(true),
                );
            plot.set_layout(layout);

            let plot_name_path = namespace_subdir.join(format!("{}_bubble.html", taxon_name));
            plot.write_html(plot_name_path); 
                                             
            Ok(())
        })?; 

    Ok(())
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

pub fn prepare_network_data<R>(
    significant_results: &FxHashMap<String, FxHashMap<GOTermID, R>>,
    go_term_to_protein_set: &FxHashMap<String, FxHashMap<GOTermID, FxHashSet<Protein>>>,
    ontology: &OboMap,
) -> FxHashMap<String, FxHashMap<NameSpace, GOTermToProteinSet>>
where
    R: EnrichmentResult + Clone + Send + Sync,
{
    significant_results
        .par_iter()
        .map(|(taxon_name, enriched_go_terms_map)| {
            let taxon_go_to_proteins = go_term_to_protein_set.get(taxon_name).unwrap();

            let network_data_by_namespace = enriched_go_terms_map
                .iter()
                .filter_map(|(go_id, _result)| {
                    let obo_term = ontology.get(go_id).unwrap();
                    let namespace = obo_term.namespace.clone();

                    taxon_go_to_proteins.get(go_id).map(|protein_set_for_go_term| {
                        (*go_id, namespace, protein_set_for_go_term.clone())
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

pub fn build_networks(
    network_data: &FxHashMap<String, FxHashMap<NameSpace, GOTermToProteinSet>>,
    enrichment_results: &FxHashMap<String, FxHashMap<GOTermID, GOTermResults>>,
) -> FxHashMap<String, FxHashMap<NameSpace, GoTermNetworkGraph>> {
    network_data
        .par_iter()
        .filter_map(|(taxon_name, taxon_specific_network_data)| {
            enrichment_results
                .get(taxon_name)
                .map(|taxon_specific_enrichment_results| {
                    let mut taxon_networks_graphs: FxHashMap<NameSpace, GoTermNetworkGraph> =
                        FxHashMap::default();

                    for current_namespace in NameSpace::iter() {
                        let go_term_proteins_in_namespace = match taxon_specific_network_data
                            .get(&current_namespace)
                        {
                            Some(data) => data,
                            None => {
                                taxon_networks_graphs
                                    .insert(current_namespace.clone(), StableGraph::default());
                                continue;
                            }
                        };

                        let mut current_term_network: GoTermNetworkGraph = StableGraph::default();
                        let mut term_to_node_index_map: FxHashMap<GOTermID, NodeIndex> =
                            FxHashMap::default();
                        let mut term_to_proteins_map_for_nodes: FxHashMap<
                            GOTermID,
                            &FxHashSet<Protein>,
                        > = FxHashMap::default();

                        for (go_term_id, protein_set) in go_term_proteins_in_namespace {
                            if let Some(enrichment_detail) =
                                taxon_specific_enrichment_results.get(go_term_id)
                            {
                                let node_data = NetworkNode {
                                    go_id: *go_term_id,
                                    lor: enrichment_detail.log_odds_ratio,
                                    stat_sig: enrichment_detail.p_value,
                                };

                                let node_idx = current_term_network.add_node(node_data);
                                term_to_node_index_map.insert(*go_term_id, node_idx);
                                term_to_proteins_map_for_nodes.insert(*go_term_id, protein_set);
                            }
                        }

                        if term_to_proteins_map_for_nodes.is_empty() {
                            taxon_networks_graphs
                                .insert(current_namespace.clone(), current_term_network);
                            continue;
                        }

                        let mut protein_to_terms_map: FxHashMap<&Protein, FxHashSet<GOTermID>> =
                            FxHashMap::default();
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

                        let mut candidate_go_pairs: FxHashSet<(GOTermID, GOTermID)> =
                            FxHashSet::default();
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

                            if size1 == 0 || size2 == 0 {
                                continue;
                            }
                            
                            if (size1.min(size2) as f64) / (size1.max(size2) as f64) <= 0.5 {
                                continue;
                            }

                            let proteins1 = match term_to_proteins_map_for_nodes.get(&term1_id) {
                                Some(p_set) => p_set,
                                None => continue, 
                            };
                            let proteins2 = match term_to_proteins_map_for_nodes.get(&term2_id) {
                                Some(p_set) => p_set,
                                None => continue,
                            };

                            let intersection_size = proteins1.intersection(proteins2).count();
                            if intersection_size == 0 {
                                continue;
                            }

                            let union_size = size1 + size2 - intersection_size;
                            if union_size == 0 {
                                continue;
                            } 

                            let jaccard_similarity: JaccardIndex =
                                (intersection_size as f64) / (union_size as f64);

                            if jaccard_similarity >= 0.5 {
                                if let (Some(&node_idx1), Some(&node_idx2)) = (
                                    term_to_node_index_map.get(&term1_id),
                                    term_to_node_index_map.get(&term2_id),
                                ) {
                                    current_term_network.add_edge(node_idx1, node_idx2, jaccard_similarity);
                                }
                            }
                        }
                        taxon_networks_graphs.insert(current_namespace.clone(), current_term_network);
                    }
                    (taxon_name.clone(), taxon_networks_graphs)
                })
        })
        .collect()
}
