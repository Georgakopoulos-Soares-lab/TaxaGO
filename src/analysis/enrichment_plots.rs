use rustc_hash::{FxHashMap, FxHashSet};
use std::path::PathBuf;
use std::error::Error;
use std::fs;
use plotly::{
    Plot, Bar, Layout, Scatter, Trace,
    common::{
        Title, Font, HoverInfo,
        ColorScale, ColorScalePalette,
        Marker, ColorBar, Anchor, Side,
        ThicknessMode, Orientation, Mode,
        Line
    },
    layout::{
        Axis, Margin, Legend,
        DragMode, RangeMode, Annotation,
        ItemClick
    },
    color::{Rgb, NamedColor, Rgba},
    // ImageFormat
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
    Directed,
    graph::NodeIndex,
    stable_graph::StableGraph,
    visit::{EdgeRef, IntoEdgeReferences}
}; 
use rayon::prelude::*;
use itertools::Itertools;
use strum::IntoEnumIterator;
use std::collections::VecDeque;
use fdg::{
    init_force_graph_uniform,
    ForceGraph,
    fruchterman_reingold::{
        FruchtermanReingold, 
        FruchtermanReingoldConfiguration
    },
    Force,
    simple::Center,
    
};

const PLOT_WIDTH: f32 = 10.0;
const PLOT_HEIGHT: f32 = 6.0;
const COLS: usize = 2;
const ROWS: usize = 2;
const QUADRANTS: usize = COLS * ROWS;

const QUADRANT_WIDTH: f32 = PLOT_WIDTH / COLS as f32;
const QUADRANT_HEIGHT: f32 = PLOT_HEIGHT / ROWS as f32;

const QUADRANT_DEFINITIONS: [(f32, f32); QUADRANTS] = [
    (0.0 * QUADRANT_WIDTH, 1.0 * QUADRANT_HEIGHT),
    (1.0 * QUADRANT_WIDTH, 1.0 * QUADRANT_HEIGHT),
    (0.0 * QUADRANT_WIDTH, 0.0 * QUADRANT_HEIGHT),
    (1.0 * QUADRANT_WIDTH, 0.0 * QUADRANT_HEIGHT),
];

pub type JaccardIndex = f32;
pub type GoTermNetworkGraph = StableGraph<GOTermPlotData, JaccardIndex, Directed>;
pub type LayoutGraph = ForceGraph<f32, 2, GOTermPlotData, JaccardIndex, Directed>;

#[derive(Debug, Clone)]
pub struct GOTermPlotData {
    pub go_id: GOTermID,
    pub name: String,
    pub wrapped_name: String,
    pub lor: f64,
    pub stat_sig: f64,
    pub minus_log10_p_value: f64,
    pub size_statistic: usize,
    pub namespace: NameSpace,
    pub hover_text: String,
}

pub trait EnrichmentResult {
    fn log_odds_ratio(&self) -> f64;
    fn p_value(&self) -> f64;
    fn size(&self) -> usize;
}

impl EnrichmentResult for GOTermResults {
    fn log_odds_ratio(&self) -> f64 {
        self.log_odds_ratio
    }

    fn p_value(&self) -> f64 {
        self.p_value
    }
    fn size(&self) -> usize {
        self.contingency_table[0]
    }
}

impl EnrichmentResult for TaxonomyGOResult {
    fn log_odds_ratio(&self) -> f64 {
        self.log_odds_ratio
    }

    fn p_value(&self) -> f64 {
        self.p_value
    }

    fn size(&self) -> usize {
        self.species_number
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
                species_map.get(taxon_name).cloned().unwrap()
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
fn create_legend_trace(
    dummy_x: Vec<std::option::Option<f64>>,
    dummy_y: Vec<std::option::Option<f64>>,
    max_bubble_size: usize,
) -> Box<dyn Trace> {
    Scatter::new(dummy_x, dummy_y)
        .mode(Mode::Markers)
        .marker(Marker::new().size(max_bubble_size).color(NamedColor::Black))
        .name(format!("{}", max_bubble_size))
        .text_font(Font::new().size(8))
        .legend_group("sizes")
        .show_legend(true)
}

fn get_namespace_subdir(namespace: &NameSpace, plots_dir: &PathBuf) -> Result<PathBuf, Box<dyn Error + Send + Sync>> {
    let namespace_str: String = match namespace {
        NameSpace::BiologicalProcess => "Biological_Process".to_string(),
        NameSpace::MolecularFunction => "Molecular_Function".to_string(),
        NameSpace::CellularComponent => "Cellular_Component".to_string(),
    };
    let namespace_subdir: PathBuf = plots_dir.join(&namespace_str);
    fs::create_dir_all(&namespace_subdir)?;
    Ok(namespace_subdir)
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
) -> FxHashMap<String, FxHashMap<NameSpace, Vec<GOTermPlotData>>>
where
    R: EnrichmentResult + Clone + Send + Sync
{
    significant_results
        .par_iter()
        .filter_map(|(species_name, go_term_results_map)| {
            let mut terms_by_namespace: FxHashMap<NameSpace, Vec<GOTermPlotData>> = FxHashMap::default();

            for (go_id, results) in go_term_results_map {
                if let Some(obo_term) = ontology.get(go_id) {
                    let current_p_value = results.p_value();
                    let current_lor = results.log_odds_ratio();

                    let minus_log_10_p = if current_p_value > 0.0 {-current_p_value.log10()} else {0.0};
                    let go_id_string = format!("GO:{:07}", go_id);
                    let original_name = obo_term.name.clone();
                    let wrapped_display_name = wrap_text(&original_name, 30);
                    let size_stat = results.size();
                    let term_namespace = obo_term.namespace.clone();

                    let hover_html_content = format!(
                        "<b>Term Name:</b> {}<br><b>Term ID:</b> {}<br><b>log(Odds Ratio):</b> {:.3}<br><b>-log10(Stat. Sig.):</b> {:.3}",
                        original_name,
                        go_id_string,
                        current_lor,
                        minus_log_10_p,
                    );

                    let rich_term = GOTermPlotData {
                        go_id: *go_id,
                        name: original_name,
                        wrapped_name: wrapped_display_name,
                        lor: current_lor,
                        stat_sig: current_p_value,
                        minus_log10_p_value: minus_log_10_p,
                        size_statistic: size_stat,
                        namespace: term_namespace.clone(),
                        hover_text: hover_html_content,
                    };

                    terms_by_namespace
                        .entry(term_namespace)
                        .or_insert_with(Vec::new)
                        .push(rich_term);
                }
            }

            for terms_vec in terms_by_namespace.values_mut() {
                terms_vec.sort_by(|a, b| {
                    a.stat_sig.partial_cmp(&b.stat_sig).unwrap_or(Equal)
                });
            }
            Some((species_name.clone(), terms_by_namespace))
        })
        .collect()
}

pub fn bar_plot(
    plot_data_map: &FxHashMap<String, FxHashMap<NameSpace, Vec<GOTermPlotData>>>,
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

            let namespace_subdir=  get_namespace_subdir(namespace, plots_dir)?;
            let mut top_20_terms: Vec<GOTermPlotData> = namespace_plot_data
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
                    .left(50)
                    .right(0)
                    .top(30)
                    .bottom(0))
                .x_axis(
                    Axis::new()
                        .title(Title::with_text("log(Odds Ratio)").font(Font::new().size(14)))
                        .tick_font(Font::new().size(12))
                        .show_line(true)
                        .line_color(NamedColor::Black)
                        .show_grid(true)
                        .grid_color(Rgba::new(0, 0, 0, 0.05))
                        .show_tick_labels(true)
                        .auto_margin(true),
                )
                .y_axis(
                    Axis::new() 
                        .title(Title::with_text(""))
                        .tick_font(Font::new().size(14))
                        .show_line(true)
                        .line_color(NamedColor::Black)
                        .show_grid(true)
                        .grid_color(Rgba::new(0, 0, 0, 0.05))
                        .show_tick_labels(true)
                        .auto_margin(true),
                )
                .drag_mode(DragMode::False)
                .bar_gap(0.4);
            plot.set_layout(layout);

            let html_file = namespace_subdir.join(format!("{}_bar_plot.html", taxon_name));
            plot.write_html(html_file); 
            
            // let svg_file = namespace_subdir.join(format!("{}_bar_plot.svg", taxon_name));
            // plot.write_image(svg_file, ImageFormat::SVG, 940, 460, 1.0);
                                             
            Ok(())
        })?; 

    Ok(())
}

pub fn bubble_plot(
    plot_data_map: FxHashMap<String, FxHashMap<NameSpace, Vec<GOTermPlotData>>>,
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
            let namespace_subdir = get_namespace_subdir(&namespace, plots_dir)?;

            let enrichment_values: Vec<f64> = namespace_plot_data.iter().map(|t| t.lor).collect();
            let stat_sig_values: Vec<f64> = namespace_plot_data.iter().map(|t| t.minus_log10_p_value).collect();
            let hover_texts: Vec<String> = namespace_plot_data.iter().map(|t| t.hover_text.clone()).collect();
            let size_statistics: Vec<usize> = namespace_plot_data.iter().map(|t| t.size_statistic).collect();

            let min_bubble_size: f64 = 10.0;
            let max_bubble_size: f64 = 25.0;

            let min_stat: f64 = *size_statistics.iter().min().unwrap() as f64;
            let max_stat: f64 = *size_statistics.iter().max().unwrap() as f64;

            let bubble_sizes: Vec<usize> = size_statistics
                .iter()
                .map(|&stat| {
                    let scaled_size_f64 = if max_stat == min_stat {
                        min_bubble_size + (max_bubble_size - min_bubble_size) / 2.0
                    } else {
                        let normalized_size = (stat as f64 - min_stat) / (max_stat - min_stat);
                        min_bubble_size + (normalized_size * (max_bubble_size - min_bubble_size))
                    };
                    scaled_size_f64.round() as usize
                })
                .collect();

            let min_bubble_size: usize = *bubble_sizes.iter().min().unwrap();
            let max_bubble_size: usize = *bubble_sizes.iter().max().unwrap();
            let mid_bubble_size_float: f64 = (min_bubble_size as f64 + max_bubble_size as f64) / 2.0;
            let mid_bubble_size: usize = mid_bubble_size_float as usize;

            let dummy_x = vec![None::<f64>];
            let dummy_y = vec![None::<f64>];

            let legend_trace_small = create_legend_trace(dummy_x.clone(), dummy_y.clone(), min_bubble_size);
            let legend_trace_medium = create_legend_trace(dummy_x.clone(), dummy_y.clone(), mid_bubble_size);
            let legend_trace_large = create_legend_trace(dummy_x.clone(), dummy_y.clone(), max_bubble_size);

            let scatter_trace = Scatter::new(enrichment_values, stat_sig_values)
                .mode(Mode::Markers)
                .marker(
                    plotly::common::Marker::new()
                        .color(Rgb::new(156, 148, 120))
                        .size_array(bubble_sizes)
                        .opacity(0.9),
                )
                .hover_text_array(hover_texts)
                .hover_info(HoverInfo::Text)
                .show_legend(false);

            let mut terms_for_sorting = namespace_plot_data.clone();
            terms_for_sorting.sort_by(|a, b| {
                b.minus_log10_p_value
                    .partial_cmp(&a.minus_log10_p_value)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            let mut top_10_significant_terms: Vec<GOTermPlotData> = terms_for_sorting
                .iter()
                .take(10)
                .cloned()
                .collect();

            top_10_significant_terms.sort_by(|a, b| {
                a.lor.partial_cmp(&b.lor)
                    .unwrap_or(Equal)
            });

            let mut annotations: Vec<Annotation> = Vec::new();
            let text_positions_cycle = vec![
                (-30, 20),  // top left
                (30, 10),   // top right
                (-30, -20), // bottom left
                (30, -10),  // bottom right
            ];

            for (i, term) in top_10_significant_terms.iter().enumerate() {
                let (ax_offset, ay_offset) = text_positions_cycle[i % text_positions_cycle.len()];
                annotations.push(
                    Annotation::new()
                        .x(term.lor)
                        .y(term.minus_log10_p_value)
                        .text(format!("GO:{:07}", term.go_id))
                        .show_arrow(true)
                        .font(
                            Font::new()
                                .size(10)
                                .color(NamedColor::Black
                            ))
                        .arrow_head(2)
                        .arrow_size(1.0)
                        .arrow_width(1.1)
                        .arrow_color(NamedColor::DimGray)
                        .ax(ax_offset)
                        .ay(ay_offset)
                        .opacity(0.9)
                );
            }

            let mut plot = Plot::new();
            plot.add_trace(scatter_trace);
            plot.add_trace(legend_trace_small);
            plot.add_trace(legend_trace_medium);
            plot.add_trace(legend_trace_large);

            let layout = Layout::new()
                .width(940)
                .height(460)
                .margin(Margin::new()
                    .left(50)
                    .right(0)
                    .top(30)
                    .bottom(0))
                .x_axis(
                    Axis::new()
                        .title(Title::with_text("log(Odds Ratio)").font(Font::new().size(14)))
                        .tick_font(Font::new().size(12))
                        .show_line(true)
                        .line_color(NamedColor::Black)
                        .show_grid(true)
                        .grid_color(Rgba::new(0, 0, 0, 0.05))
                        .show_tick_labels(true)
                        .auto_margin(true)
                        .range_mode(RangeMode::ToZero),
                )
                .y_axis(
                    Axis::new()
                        .title(Title::with_text("-log10(Stat. Sig.)").font(Font::new().size(14)))
                        .tick_font(Font::new().size(12))
                        .show_line(true)
                        .line_color(NamedColor::Black)
                        .show_grid(true)
                        .grid_color(Rgba::new(0, 0, 0, 0.05))
                        .show_tick_labels(true)
                        .auto_margin(true)
                        .range_mode(RangeMode::ToZero),
                )
                .legend(
                    Legend::new()
                        .x(1.0)
                        .y(1.0)
                        .trace_group_gap(10)
                        .title(Title::with_text("GO Term size").font(Font::new().size(12)))
                        .item_click(ItemClick::False)
                        .item_double_click(ItemClick::False)
                )
                .annotations(annotations);

            plot.set_layout(layout);

            let html_file = namespace_subdir.join(format!("{}_bubble_plot.html", taxon_name));
            plot.write_html(&html_file);
            // let svg_file = namespace_subdir.join(format!("{}_bubble_plot.svg", taxon_name));
            // plot.write_image(&svg_file, ImageFormat::SVG, 940, 460, 1.0)

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
    ontology: &OboMap,
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
                        if let Some(proteins) = taxon_specific_network_data.get(&current_namespace) {
                            let go_term_proteins_in_namespace = proteins;
                        

                        let mut current_namespace_network: GoTermNetworkGraph = StableGraph::default();
                        let mut term_to_node_index_map: FxHashMap<GOTermID, NodeIndex> = FxHashMap::default();
                        let mut term_to_proteins_map_for_nodes: FxHashMap<GOTermID, &FxHashSet<Protein>> = FxHashMap::default();

                        for (go_term_id, protein_set) in go_term_proteins_in_namespace {
                            if let (Some(enrichment_detail), Some(obo_term)) = (
                                taxon_specific_enrichment_results.get(go_term_id),
                                ontology.get(go_term_id)
                            ) {
                                let p_value = enrichment_detail.p_value();
                                let lor = enrichment_detail.log_odds_ratio();
                                let minus_log10_p = if p_value > 0.0 { -p_value.log10() } else { 0.0 };
                                let go_id_str = format!("GO:{:07}", go_term_id);
                                let name = obo_term.name.clone();
                                let wrapped_name = wrap_text(&name, 30);
                                let size = enrichment_detail.size();

                                let hover_text = format!(
                                    "<b>Term Name:</b> {}<br><b>Term ID:</b> {}<br><b>log(Odds Ratio):</b> {:.3}<br><b>-log10(Stat. Sig.):</b> {:.3}",
                                    name, go_id_str, lor, minus_log10_p
                                );

                                let node_data = GOTermPlotData {
                                    go_id: *go_term_id,
                                    name: name.to_string(),
                                    wrapped_name,
                                    lor,
                                    stat_sig: p_value,
                                    minus_log10_p_value: minus_log10_p,
                                    size_statistic: size,
                                    namespace: obo_term.namespace.clone(), 
                                    hover_text,
                                };

                            let node_idx = current_namespace_network.add_node(node_data);
                            term_to_node_index_map.insert(*go_term_id, node_idx);
                            term_to_proteins_map_for_nodes.insert(*go_term_id, protein_set);
                            }
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
                            let size1 = *term_node_sizes.get(&term1_id).unwrap();
                            let size2 = *term_node_sizes.get(&term2_id).unwrap();

                            let proteins1 = term_to_proteins_map_for_nodes.get(&term1_id).unwrap();
                            let proteins2 = term_to_proteins_map_for_nodes.get(&term2_id).unwrap();

                            let intersection_size = proteins1.intersection(proteins2).count();
                            let union_size = size1 + size2 - intersection_size;

                            let jaccard_similarity: JaccardIndex =
                                (intersection_size as f32) / (union_size as f32);

                            if jaccard_similarity >= 0.25 {
                                let &node_idx1= term_to_node_index_map.get(&term1_id).unwrap();
                                let &node_idx2= term_to_node_index_map.get(&term2_id).unwrap();
                                
                                current_namespace_network.add_edge(node_idx1, node_idx2, jaccard_similarity);
                            }

                        }
                        
                        let top_k_subgraphs = extract_top_k_communities(&current_namespace_network, 4);
                        taxon_networks_graphs.insert(current_namespace.clone(), top_k_subgraphs);
                    }
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

        let mut nodes_to_remove = Vec::new();
        for node_index in subgraph.node_indices() {
            if subgraph.neighbors_undirected(node_index).count() == 0 { 
                nodes_to_remove.push(node_index);
            }
        }

        for node_index in nodes_to_remove {
            subgraph.remove_node(node_index);
        }

        top_k_graphs.push(subgraph);
    }

    top_k_graphs
}

fn apply_fruchterman_reingold_layout(
    original_graph: &GoTermNetworkGraph,
    iterations: usize,
) -> ForceGraph<f32, 2, GOTermPlotData, JaccardIndex, Directed> {

    let mut force_layout_graph=
        init_force_graph_uniform(
            original_graph.clone(),
            10.0,
    );

    let mut fr_force = FruchtermanReingold {
        conf: FruchtermanReingoldConfiguration {
            dt: 0.02,
            cooloff_factor: 0.975,
            scale: 50.0,
        },
        ..Default::default()
    };

    fr_force.apply_many(&mut force_layout_graph, iterations);
    Center::default().apply(&mut force_layout_graph);

    force_layout_graph
}

pub fn network_plot(
    top_networks_map: &FxHashMap<String, FxHashMap<NameSpace, Vec<GoTermNetworkGraph>>>,
    plots_dir: &PathBuf
) -> Result<(), Box<dyn Error + Send + Sync>> {

    let mut network_layouts_map: FxHashMap<String, FxHashMap<NameSpace, Vec<LayoutGraph>>> =
        top_networks_map
            .par_iter() 
            .map(|(taxon_name, namespace_to_networks_map)| {
                let layouts_for_taxon: FxHashMap<NameSpace, Vec<LayoutGraph>> =
                    namespace_to_networks_map
                        .iter() 
                        .map(|(namespace, networks_vec)| {
                            let layouts_for_namespace: Vec<LayoutGraph> = networks_vec
                                .iter()
                                .map(|network_graph| {
                                    apply_fruchterman_reingold_layout(network_graph, 1000)
                                })
                                .collect();
                            (*namespace, layouts_for_namespace)
                        })
                        .collect();
                (taxon_name.clone(), layouts_for_taxon)
            })
            .collect();

    network_layouts_map
        .par_iter_mut() 
        .try_for_each(|(taxon_name, namespace_map)| {
            namespace_map
                .iter_mut()
                .try_for_each(|(namespace, layouts_vec)| {
                    layouts_vec
                        .iter_mut()
                        .enumerate()
                        .for_each(|(i, graph)| {
                            let (quadrant_min_x, quadrant_min_y) = QUADRANT_DEFINITIONS[i];

                            let padding_abs_x = QUADRANT_WIDTH * 0.1;
                            let padding_abs_y = QUADRANT_HEIGHT * 0.1;

                            let drawable_origin_x = quadrant_min_x + padding_abs_x;
                            let drawable_origin_y = quadrant_min_y + padding_abs_y;
                            let mut drawable_width = QUADRANT_WIDTH - 2.0 * padding_abs_x;
                            let mut drawable_height = QUADRANT_HEIGHT - 2.0 * padding_abs_y;

                            drawable_width = drawable_width.max(0.0);
                            drawable_height = drawable_height.max(0.0);

                            let mut min_graph_x = f32::MAX;
                            let mut max_graph_x = f32::MIN;
                            let mut min_graph_y = f32::MAX;
                            let mut max_graph_y = f32::MIN;

                            graph.node_weights().for_each(|(_node_data, location)| {
                                min_graph_x = min_graph_x.min(location.x);
                                max_graph_x = max_graph_x.max(location.x);
                                min_graph_y = min_graph_y.min(location.y);
                                max_graph_y = max_graph_y.max(location.y);
                            });

                            let current_graph_width = max_graph_x - min_graph_x;
                            let current_graph_height = max_graph_y - min_graph_y;

                            let scale_ratio_x = drawable_width / current_graph_width;
                            let scale_ratio_y = drawable_height / current_graph_height;
                            let mut scale_factor = scale_ratio_x.min(scale_ratio_y);

                            scale_factor = scale_factor.max(0.0);

                            let scaled_graph_width = current_graph_width * scale_factor;
                            let scaled_graph_height = current_graph_height * scale_factor;

                            let offset_x_in_drawable = (drawable_width - scaled_graph_width) / 2.0;
                            let offset_y_in_drawable = (drawable_height - scaled_graph_height) / 2.0;

                            let final_translation_x = drawable_origin_x + offset_x_in_drawable;
                            let final_translation_y = drawable_origin_y + offset_y_in_drawable;

                            graph.node_weights_mut().for_each(|(_node_data, location)| {
                                let original_relative_x = location.x - min_graph_x;
                                let original_relative_y = location.y - min_graph_y;

                                location.x = original_relative_x * scale_factor + final_translation_x;
                                location.y = original_relative_y * scale_factor + final_translation_y;
                            });
                        }); 
                    
                    let mut plot = Plot::new();
                    const MIN_EDGE_WIDTH: f64 = 2.0;
                    const MAX_EDGE_WIDTH: f64 = 8.0; 

                    let mut all_jaccard_indices_for_this_plot: Vec<JaccardIndex> = Vec::new();
                    for graph in layouts_vec.iter() {
                        for edge_ref in graph.edge_references() {
                            all_jaccard_indices_for_this_plot.push(*edge_ref.weight());
                        }
                    }

                    let min_jaccard_opt = all_jaccard_indices_for_this_plot.iter().copied().reduce(f32::min);
                    let max_jaccard_opt = all_jaccard_indices_for_this_plot.iter().copied().reduce(f32::max);

                    if let (Some(min_j), Some(max_j)) = (min_jaccard_opt, max_jaccard_opt) {
                        let min_jaccard = min_j;
                        let max_jaccard = max_j;

                        for graph in layouts_vec.iter() {
                            for edge_ref in graph.edge_references() {
                                let jaccard_index_val = *edge_ref.weight();
                                let source_idx = edge_ref.source();
                                let target_idx = edge_ref.target();

                                if let (Some(source_node_info), Some(target_node_info)) =
                                    (graph.node_weight(source_idx), graph.node_weight(target_idx))
                                {
                                    let x_start = source_node_info.1.x as f64;
                                    let y_start = source_node_info.1.y as f64;
                                    let x_end = target_node_info.1.x as f64;
                                    let y_end = target_node_info.1.y as f64;

                                    let scaled_width = if max_jaccard <= min_jaccard { 
                                        MIN_EDGE_WIDTH + (MAX_EDGE_WIDTH - MIN_EDGE_WIDTH) / 2.0
                                    } else {
                                        let normalized_val = (jaccard_index_val - min_jaccard) as f64 / (max_jaccard - min_jaccard) as f64;
                                        MIN_EDGE_WIDTH + normalized_val * (MAX_EDGE_WIDTH - MIN_EDGE_WIDTH)
                                    };
                                    let final_edge_width = scaled_width.max(MIN_EDGE_WIDTH).min(MAX_EDGE_WIDTH);

                                    let edge_segment_trace = Scatter::new(vec![x_start, x_end], vec![y_start, y_end])
                                        .mode(Mode::Lines)
                                        .line(
                                            Line::new()
                                                .width(final_edge_width)
                                                .color(Rgba::new(200, 200, 200, 0.5))
                                                
                                        )
                                        .show_legend(false);
                                    plot.add_trace(edge_segment_trace);
                                }
                            }
                        }
                    }
                    let mut all_plot_annotations: Vec<Annotation> = Vec::new();
                    let text_positions_cycle = vec![
                        (30, -30),
                        (-30, 15),
                        (30, 30),
                        (-30, -15)
                    ];
                    let mut annotation_offset_idx_counter = 0;

                    let mut all_nodes_x: Vec<f32> = Vec::new();
                    let mut all_nodes_y: Vec<f32> = Vec::new();
                    let mut all_nodes_hover_text: Vec<String> = Vec::new();
                    let mut all_nodes_color_values: Vec<f64> = Vec::new(); 
                    let mut all_nodes_sizes: Vec<f64> = Vec::new();

                    let mut all_edges_coordinates: Vec<((f64, f64), (f64, f64))> = Vec::new();

                    for graph in layouts_vec.iter() {
                        for (node_plot_data, location) in graph.node_weights() {
                            all_nodes_x.push(location.x);
                            all_nodes_y.push(location.y);
                            all_nodes_hover_text.push(node_plot_data.hover_text.clone());
                            all_nodes_color_values.push(node_plot_data.minus_log10_p_value);
                            all_nodes_sizes.push(node_plot_data.size_statistic as f64);

                            let (x_shift, y_shift) = &text_positions_cycle[annotation_offset_idx_counter % text_positions_cycle.len()];
                            annotation_offset_idx_counter += 1;

                            let annotation = Annotation::new()
                                .x(location.x as f64)
                                .y(location.y as f64)
                                .text(format!("GO:{:07}", node_plot_data.go_id))
                                .show_arrow(true)
                                .font(Font::new().size(10).color(NamedColor::Black))
                                .ax(*x_shift)
                                .ay(*y_shift)
                                .opacity(0.9);
                            all_plot_annotations.push(annotation);
                        }

                        for edge_ref in graph.edge_references() {
                            let source_idx = edge_ref.source();
                            let target_idx = edge_ref.target();

                            if let (Some(source_node_info), Some(target_node_info)) =
                                (graph.node_weight(source_idx), graph.node_weight(target_idx))
                            {
                                let source_coords = (source_node_info.1.x as f64, source_node_info.1.y as f64);
                                let target_coords = (target_node_info.1.x as f64, target_node_info.1.y as f64);
                                all_edges_coordinates.push((source_coords, target_coords));
                            }
                        }
                    }
                    
                    let min_size: f64 = 15.0;
                    let max_size: f64 = 35.0;

                    let min_stat: f64 = all_nodes_sizes
                        .iter()
                        .copied()
                        .reduce(f64::min)
                        .unwrap();

                    let max_stat: f64 = all_nodes_sizes
                        .iter()
                        .copied()
                        .reduce(f64::max)
                        .unwrap();

                    let node_sizes: Vec<usize> = all_nodes_sizes
                        .iter()
                        .map(|&stat| {
                            let scaled_size_f64 = if max_stat == min_stat {
                                min_size + (max_size - min_size) / 2.0
                            } else {
                                let normalized_size = (stat as f64 - min_stat) / (max_stat - min_stat);
                                min_size + (normalized_size * (max_size - min_size))
                            };
                            scaled_size_f64.round() as usize
                        })
                        
                        .collect();

                    let color_bar = ColorBar::new()
                        .title(
                            Title::from("log(Odds Ratio)")
                                .side(Side::Right)
                                .font(Font::new().size(12)),
                        )
                        .tick_font(Font::new().size(10))
                        .len_mode(ThicknessMode::Pixels)
                        .len(200)
                        .thickness(15)
                        .x(1.0)
                        .y(0.8)
                        .y_anchor(Anchor::Middle);

                    let node_trace = Scatter::new(all_nodes_x, all_nodes_y)
                        .mode(Mode::Markers)
                        .marker(
                            Marker::new()
                                .color_array(all_nodes_color_values)
                                .color_scale(ColorScale::Palette(ColorScalePalette::Viridis))
                                .color_bar(color_bar)
                                .size_array(node_sizes)
                                .show_scale(true)
                                .opacity(1.0)
                            )
                        .hover_text_array(all_nodes_hover_text) 
                        .hover_info(HoverInfo::Text) 
                        .show_legend(false);

                    let mut edge_x_coords: Vec<Option<f64>> = Vec::new();
                    let mut edge_y_coords: Vec<Option<f64>> = Vec::new();

                    for (i, edge) in all_edges_coordinates.iter().enumerate() {
                        let ((x_start, y_start), (x_end, y_end)) = edge;

                        edge_x_coords.push(Some(*x_start));
                        edge_y_coords.push(Some(*y_start));
                        edge_x_coords.push(Some(*x_end));
                        edge_y_coords.push(Some(*y_end));

                        if i < all_edges_coordinates.len() - 1 {
                            edge_x_coords.push(None);
                            edge_y_coords.push(None);
                        }
                    }

                    plot.add_trace(node_trace);

                    let layout = Layout::new()
                        .width(940)
                        .height(460)
                        .margin(Margin::new()
                            .left(50)
                            .right(0)
                            .top(30)
                            .bottom(0))
                        .x_axis(
                            Axis::new()
                                .show_line(false)
                                .zero_line(false)
                                .show_grid(true)
                                .show_tick_labels(false)
                                .auto_margin(true)
                        )
                        .y_axis(
                            Axis::new() 
                                .show_line(false)
                                .zero_line(false)
                                .show_grid(true)
                                .show_tick_labels(false)
                                .auto_margin(true)
                        )
                        .annotations(all_plot_annotations);
                    plot.set_layout(layout);

                    let namespace_subdir = get_namespace_subdir(namespace, plots_dir)?;
                    let html_file = namespace_subdir.join(format!("{}_network_plot.html", taxon_name));
                    plot.write_html(html_file); 
                    // let svg_file = namespace_subdir.join(format!("{}_network_plot.svg", taxon_name));
                    // plot.write_image(svg_file, ImageFormat::SVG, 940, 460, 1.0);

                    Ok::<(), Box<dyn Error + Send + Sync>>(())
                })
        })?;

Ok(())
}
