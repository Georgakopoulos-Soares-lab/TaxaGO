use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use clap::{Parser, ValueEnum, ArgGroup};
use std::fs;
use std::env::var;
use std::process::ExitCode;
use rustc_hash::{FxHashMap, FxHashSet};
use daggy::NodeIndex;
use std::path::PathBuf;
use dirs::home_dir;

use TaxaGO::parsers::{
    background_parser::*, obo_parser::*, study_parser::*
};
use TaxaGO::analysis::{
    enrichment_analysis::*, 
    multiple_testing_correction::*, 
    write_results::*,
    handle_lineage::*,
    result_combination::*,
    count_propagation::*,
    phylogenetic_meta_analysis::*,
    enrichment_plots::*
};

fn get_default_asset_path(filename: &str) -> String {
    let cargo_home = var("CARGO_HOME")
        .unwrap_or_else(|_| {
            home_dir()
                .expect("Could not determine home directory")
                .join(".cargo")
                .to_string_lossy()
                .into_owned()
        });
    PathBuf::from(cargo_home)
        .join("taxago_assets")
        .join(filename)
        .to_string_lossy()
        .into_owned()
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum PropagationMethod {
    None,
    Classic,
    Elim,
    Weight
}

#[derive(Parser, Debug)]
#[command(name = "taxago", about, version, author)]
#[command(group(
    ArgGroup::new("meta_analysis")
        .requires("vcv_matrix")
        .requires("combine_results")
))]
struct CliArgs {
    #[arg(
        short = 'o',
        long = "obo",
        value_name = "FILE",
        help = "Path to the Gene Ontology file in OBO format.",
    )]
    obo_file: Option<String>,
    
    #[arg(
        short = 's',
        long = "study",
        value_name = "FILE_OR_DIR",
        help = "Directory containing study population for each taxon in FASTA format or CSV file with the study population for each species.",
        required = true
    )]
    study_pop: String,
    
    #[arg(
        short = 'b',
        long = "background",
        value_name = "DIRECTORY",
        help = "Directory containing background populations.",
    )]
    background_pop: Option<String>,

    #[arg(
        short = 'e',
        long = "evidence",
        value_name = "CATEGORY",
        help = "Evidence code categories to parse in background associations. [possible values: all, experimental, phylogenetic, computational, author, curator, electronic]",
        default_value = "all"
    )]
    evidence_categories: String,   
    
    #[arg(
        short = 'd',
        long = "dir",
        value_name = "DIRECTORY",
        help = "Directory to write results for each taxon and the combined results.",
        required = true
    )]
    output_dir: PathBuf,
    
    #[arg(
        short = 'p',
        long = "propagate-counts",
        value_enum,
        help = "Method to propagate GO term counts up the Ontology graph.",
        default_value_t = PropagationMethod::None
    )]
    propagate_counts: PropagationMethod,

    #[arg(
        short = 't',
        long = "test",
        value_enum,
        help = "Statistical test to use.",
        default_value_t = StatisticalTest::Fishers
    )]
    statistical_test: StatisticalTest,

    #[arg(
        short = 'm',
        long = "min-prot",
        value_name = "COUNT",
        help = "Minimum protein count a GO Term must have to be processed.",
        default_value_t = 5
    )]
    min_protein_count: usize,
    
    #[arg(
        short = 'r',
        long = "min-score",
        value_name = "SCORE",
        help = "Minimum log(odds ratio) a GO Term must have to be included in results.",
        default_value_t = 0.2
    )]
    min_odds_ratio: f64,

    #[arg(
        short = 'a',
        long = "alpha",
        value_name = "THRESHOLD",
        help = "P-value threshold to determine significant results.",
        default_value_t = 0.05
    )]
    significance_threshold: f64,

    #[arg(
        short = 'c',
        long = "correction-method",
        value_enum,
        help = "Method to adjust p-values for multiple test correction.",
        default_value_t = AdjustmentMethod::BenjaminiHochberg
    )]
    correction_method: AdjustmentMethod,

    #[arg(
        short = 'g',
        long = "group-results",
        value_name = "TAXONOMIC_LEVEL",
        help = "Combine results based on the specified taxonomic level.",
    )]
    combine_results: Option<String>,
    
    #[arg(
        short = 'l',
        long = "lineage-percentage",
        value_name = "PERCENTAGE",
        help = "Percentage of species inside the taxonomic level where GO term must be found.",
        default_value_t = 0.25
    )]
    lineage_percentage: f64,

    #[arg(
        long = "vcv-matrix",
        value_name = "FILE",
        help = "Path to a custom variance-covariance matrix file for phylogenetic meta-analysis.",
    )]
    vcv_matrix: Option<PathBuf>,

    #[arg(
        long = "permutations",
        value_name = "COUNT",
        help = "Number of permutations for phylogenetic meta-analysis.",
        default_value_t = 10000
    )]
    permutations: u32,

    #[arg(
        long = "cores",
        value_name = "NUMBER",
        help = "Number of cores to use for the analysis. Uses all available by default.",
        default_value_t = num_cpus::get()
    )]
    num_cores: usize,

    #[arg(
        long = "save-plots",
        help = "If specified, TaxaGO will save the enrichment plots.",
        value_enum,
        default_value_t = PlotType::Interactive
    )]
    save_plots: PlotType,
}

fn main() -> ExitCode{
    let cli_args: CliArgs = CliArgs::parse();
    
    let default_obo_path = get_default_asset_path("go.obo");
    let default_background_path = get_default_asset_path("background_pop");
    let default_lineage_path = get_default_asset_path("lineage.txt");
    
    let obo_file = cli_args.obo_file.unwrap_or(default_obo_path);
    let background_pop = cli_args.background_pop.unwrap_or(default_background_path);
    
    let cargo_home = var("CARGO_HOME")
            .unwrap_or_else(|_| {
                home_dir()
                    .expect("Could not determine home directory")
                    .join(".cargo")
                    .to_string_lossy()
                    .into_owned()
            });

    println!("\nAnalysis will be performed with {} core(s)", &cli_args.num_cores);

    if let Err(e) = rayon::ThreadPoolBuilder::new()
        .num_threads(cli_args.num_cores)
        .build_global() {
        eprintln!("Failed to initialize Rayon global thread pool: {:?}", e);
    }; 
    println!("\nCleaning previous results");
    clean_directory(&cli_args.output_dir).unwrap_or_else(|e| {
        eprintln!("Error cleaning output directory: {}", e);
    });

    fs::create_dir_all(&cli_args.output_dir).unwrap_or_else(|e| {
        eprintln!("Error creating output directory: {}", e);
    });
    
    let obo_file_path = PathBuf::from(&obo_file);

    println!("\nReading ontology information from: {}", &obo_file_path.to_string_lossy());

     let ontology = match parse_obo_file(&obo_file_path) {
        Ok(parsed_ontology) => parsed_ontology,
        Err(e) => {
            eprintln!("\nError processing OBO file '{}':", obo_file);
            eprintln!("{}", e);
            return ExitCode::FAILURE; 
        }
    };
    
    println!("Successfully parsed OBO file with {} GO terms\n", ontology.len());

    let (ontology_graph, go_id_to_node_index) = match build_ontology_graph(&ontology) {
        Ok(graph_data) => graph_data,
        Err(e) => {
            eprintln!("\nError building ontology graph from OBO file '{}':", obo_file);
            eprintln!("{}", e);
            return ExitCode::FAILURE;
        }
    };

    let node_index_to_go_id: FxHashMap<NodeIndex, GOTermID> = go_id_to_node_index.iter()
            .map(|(go_term, node_index)| (*node_index, go_term.clone()))
            .collect();

    let root_go_ids: Vec<u32> = vec![8150, 3674, 5575];

    let (_, level_to_go_term) = assign_levels_from_roots(
        &ontology_graph,
        &go_id_to_node_index,
        &node_index_to_go_id,
        &root_go_ids
    );
    
    println!("Reading background populations from: {}\n", &background_pop);
    let taxon_ids: FxHashSet<TaxonID> = match collect_taxon_ids(&PathBuf::from(&cli_args.study_pop)) {
        Ok(taxon_ids) => taxon_ids,
        Err(e) => {
            eprintln!("\nError collecting taxon IDs from study population '{}':", cli_args.study_pop);
            eprintln!("{}", e);
            return ExitCode::FAILURE;
        }
    };
    let categories: Vec<EvidenceCategory> = match map_input_to_category(cli_args.evidence_categories.clone()) {
        Ok(categories) => categories,
        Err(e) => {
            eprintln!("\nError parsing evidence categories '{}':", cli_args.evidence_categories);
            eprintln!("{}", e);
            return ExitCode::FAILURE;
        }
    };

    let mut background_population = match BackgroundPop::read_background_pop(
        &taxon_ids, 
        &background_pop,
        &categories
    ) {
        Ok(Some(background_pop)) => {
            println!("Successfully loaded background population for {} taxa\n", &taxon_ids.len());
            background_pop
        },
        Ok(None) => {
            eprintln!("Error: No background population data could be loaded from directory '{}'", background_pop);
            return ExitCode::FAILURE;
        },
        Err(e) => {
            eprintln!("Error reading background population: {}", e);
            return ExitCode::FAILURE;
        }
    };

    println!("Reading study populations from: {}\n", &cli_args.study_pop.clone());
    
    let mut study_population = match StudyPop::read_study_pop(
        &PathBuf::from(&cli_args.study_pop), 
        &background_population.protein_to_go
    ){
        Ok(Some(study_pop)) => {
            println!("Successfully loaded study population for {} taxa\n", &taxon_ids.len());
            study_pop
        },
        Ok(None) => {
            eprintln!("Error: No study population data could be loaded from '{}'", cli_args.study_pop);
            return ExitCode::FAILURE;
        },
        Err(e) => {
            eprintln!("Error reading study population: {}", e);
            return ExitCode::FAILURE;
        }
    };

    let should_propagate = match cli_args.propagate_counts {
        PropagationMethod::None => false,
        PropagationMethod::Classic | PropagationMethod::Elim | PropagationMethod::Weight => true,
    };

    if should_propagate{
        println!("Propagating counts up the Ontology graph\n");
        
        let ancestor_cache: GOAncestorCache = match GOAncestorCache::new(
            &ontology_graph, 
            &ontology, 
            &go_id_to_node_index,
            &node_index_to_go_id) {
                Ok(cache) => cache,
                Err(e) => {
                    eprintln!("\nError creating GO ancestor cache: {}", e);
                    return ExitCode::FAILURE;
                }};

        study_population.propagate_counts(
            &taxon_ids, 
            &ancestor_cache
        );
        background_population.propagate_counts(
            &taxon_ids, 
            &ancestor_cache
        );
        
    }
    
    study_population.filter_by_threshold(
        &taxon_ids,
        cli_args.min_protein_count
    );
    
    background_population.filter_by_study_population(
        &taxon_ids, 
        &study_population
    );
    
    println!("Starting Gene Ontology (GO) term enrichment analysis\n");
    
    let analysis = EnrichmentAnalysis::new(cli_args.statistical_test);

    let enrichment_results = match cli_args.propagate_counts {
        PropagationMethod::Elim => {
            println!("Performing elim algorithm on propagated counts\n");
            analysis.elim_analysis(
                &taxon_ids,
                cli_args.significance_threshold,
                &study_population,
                &background_population,
                &level_to_go_term
            )
        },
        PropagationMethod::Classic => {
            println!("Performing classic analysis with propagated counts\n");
            analysis.classic(
                &taxon_ids,          
                &background_population.go_term_count,
                &study_population.go_term_count,
                &background_population.taxon_protein_count,
                &study_population.taxon_protein_count,
            )
        },
        PropagationMethod::Weight => {
            println!("Performing weight algorithm with propagated counts\n");
            analysis.weight(
                &taxon_ids, 
                &study_population, 
                &background_population, 
                &level_to_go_term, 
                &ontology_graph, 
                &go_id_to_node_index, 
                &node_index_to_go_id)
        },
        PropagationMethod::None => {
            println!("Performing classic analysis without count propagation\n");
            analysis.classic(
                &taxon_ids,          
                &background_population.go_term_count,
                &study_population.go_term_count,
                &background_population.taxon_protein_count,
                &study_population.taxon_protein_count,
            )
        }
    };
    
    let significant_species_results = adjust_species_p_values(
        &enrichment_results, 
        cli_args.correction_method, 
        Some(cli_args.significance_threshold),
        cli_args.min_odds_ratio
    );
        
    let taxid_species_map = match taxid_to_species(default_lineage_path.clone()) {
        Ok(map) => map,
        Err(e) => {
            eprintln!("\nError reading taxonomic lineage information from '{}':", default_lineage_path);
            eprintln!("{}", e);
            return ExitCode::FAILURE;
        }
    };

    match write_single_taxon_results(
        &significant_species_results,
        &ontology,
        &taxid_species_map,
        &cli_args.output_dir,
    ) {
        Ok(_) => {
        }
        Err(e) => {
            eprintln!(
                "Error: Failed to write single taxon results to directory '{}': {}",
                cli_args.output_dir.display(), 
                e 
            );
            return ExitCode::FAILURE; 
        }
    }
    if cli_args.save_plots != PlotType::None {
        println!("Generating enrichment plots\n");
        let species_plots_subdir = cli_args.output_dir.join("single_taxon_results").join("plots");
        fs::create_dir_all(&species_plots_subdir).unwrap_or_else(|e| {
            eprintln!("Error creating species plot  directory: {}", e);
        });

        let (processed_species_data, go_term_to_protein_set) = process_species_data(
            significant_species_results,
            &study_population,
            &taxid_species_map
        );

        let species_plot_data = prepare_plot_data(
            &processed_species_data, 
            &ontology);
        
        let _species_bar_plots = bar_plot(
            &species_plot_data, 
            &species_plots_subdir,
            cli_args.save_plots
        );
        
        let _species_bubble_plots = bubble_plot(
            species_plot_data, 
            &species_plots_subdir,
            cli_args.save_plots);
        
        let species_protein_provider = ProteinDataProvider::Species(&go_term_to_protein_set);
        let species_network_data = prepare_network_data(
            &processed_species_data,
            &species_protein_provider, 
            &ontology,
        );

        let species_networks = build_networks(
            &species_network_data,
            &processed_species_data,
            &ontology
        );
        
        let _species_network_plots = network_plot(
            &species_networks, 
            &species_plots_subdir,
            cli_args.save_plots);
    }  
    
    if let Some(level_to_combine) = &cli_args.combine_results {

        println!("Reading taxonomic lineage information from: {}\n", default_lineage_path);
        
        let lineage = match read_lineage(default_lineage_path.clone()) {
            Ok(lineage) => lineage,
            Err(e) => {
                eprintln!("\nError reading taxonomic lineage information from '{}':", default_lineage_path);
                eprintln!("{}", e);
                return ExitCode::FAILURE;
            }
        };

        let superkingdom = match get_superkingdom(&taxon_ids, &lineage) {
            Ok(superkingdom) => superkingdom,
            Err(e) => {
                eprintln!("\nError determining superkingdom for taxon IDs: {}", e);
                return ExitCode::FAILURE;
            }
        };
        
        println!("Grouping species based on {}\n", level_to_combine);
        
        let grouped_species = taxid_to_level(
            &enrichment_results,
            &lineage,
            level_to_combine
        );

        let lineage_organized_results= group_results_by_taxonomy(
            &grouped_species, 
            &enrichment_results, 
            cli_args.lineage_percentage
        );

        let matrix_path = if let Some(custom_path) = &cli_args.vcv_matrix {
            println!("Using custom VCV matrix from: {:?} \n", custom_path);
            custom_path.clone()
        } else {
            let matrix_filename = format!("{}.dmat", &superkingdom);
            let default_path = PathBuf::from(&cargo_home)
                .join("taxago_assets")
                .join(matrix_filename);
            println!("Reading {} VCV matrix from: {:?} \n", &superkingdom, default_path);
            default_path
        };

        let vcv_matrix = read_vcv_matrix(
            matrix_path
        ).unwrap();

        println!("Performing phylogenetic meta-analysis with {} permutations", &cli_args.permutations);

        let phylogenetic_results = phylogenetic_meta_analysis(
            &taxon_ids,
            lineage_organized_results, 
            vcv_matrix,
            cli_args.permutations
        );

        let significant_taxonomy_results = adjust_taxonomy_p_values(
            &phylogenetic_results, 
            cli_args.correction_method, 
            Some(cli_args.significance_threshold),
            cli_args.min_odds_ratio,
            level_to_combine);
        
        match write_taxonomy_results(
        &significant_taxonomy_results,
            &ontology,
            &cli_args.output_dir,
            level_to_combine
        ) {
            Ok(_) => {
            }
            Err(e) => {
                eprintln!(
                    "Error: Failed to write taxonomy results to directory '{}': {}",
                    cli_args.output_dir.display(), 
                    e 
                );
                return ExitCode::FAILURE; 
            }
        }

        if cli_args.save_plots != PlotType::None{
            let taxonomy_plots_subdir = cli_args.output_dir.join("combined_taxonomy_results").join("plots");
            fs::create_dir_all(&taxonomy_plots_subdir).unwrap_or_else(|e| {
                eprintln!("Error creating taxonomy plot directory: {}", e);
            });

            let taxonomy_plot_data = prepare_plot_data(
                &significant_taxonomy_results, 
                &ontology);

            let _taxonomy_bar_plots = bar_plot(
                &taxonomy_plot_data, 
                &taxonomy_plots_subdir,
                cli_args.save_plots);

            let _taxonomy_bubble_plots = bubble_plot(
                taxonomy_plot_data, 
                &taxonomy_plots_subdir,
                cli_args.save_plots);

            let taxonomy_protein_provider = ProteinDataProvider::Taxonomy {
                species_data_by_id: &study_population.go_term_to_protein_set,
                taxonomy_to_species_ids: &grouped_species,
            };
            let taxon_network_data = prepare_network_data(
                &significant_taxonomy_results,
                &taxonomy_protein_provider,
                &ontology,
            );
            
            let taxon_networks = build_networks(
                &taxon_network_data,
                &significant_taxonomy_results,
                &ontology
            );
            
            let _taxon_network_plots = network_plot(
                &taxon_networks, 
                &taxonomy_plots_subdir,
                cli_args.save_plots);
            }

        }
    
    println!("Finished analysis\n");
    ExitCode::SUCCESS
}