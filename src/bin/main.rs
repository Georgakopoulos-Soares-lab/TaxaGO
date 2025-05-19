use actix_web::{web, App, HttpResponse, HttpServer, Result, get, post, Responder};
use std::sync::Arc;
use std::path::Path;
use std::fs;
use std::io::{BufReader, BufRead, Write};
use serde::{Deserialize};
use log::{info, warn, error, debug};
use std::fs::{OpenOptions, File};
use std::process::Command;
use serde_json::json;
use ctrlc;
use actix_multipart::Multipart;
use futures::{StreamExt, TryStreamExt};
use std::collections::HashMap;
use TaxaGO::utils::download;

#[derive(Deserialize)]
struct QueryParams {
    odds_ratio_threshold: Option<f64>,
    n_protein_threshold: Option<i32>,
    selected_background: Option<String>,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
struct CategoryUpdate {
    category: String,
    level: String,
    superkingdom: String,
    #[serde(rename = "isChecked")]
    is_checked: bool,
    #[serde(rename = "isUserAction")]
    is_user_action: bool
}

#[derive(Debug, serde::Serialize)]
struct SavedCategory {
    category: String,
    level: String,
    superkingdom: String,
}

// AppState will hold our shared data
struct AppState {
}

async fn index() -> HttpResponse {
    debug!("Serving index page");
    let html_content = include_str!("templates/index.html");
    HttpResponse::Ok()
        .content_type("text/html")
        .body(html_content)
}

async fn loading(query: web::Query<QueryParams>) -> HttpResponse {
    debug!("Serving loading page");
    let html_content = include_str!("templates/loading.html");
    
    let odds_ratio = query.odds_ratio_threshold.unwrap_or(1.0);
    let n_protein = query.n_protein_threshold.unwrap_or(5);
    let html_with_params = html_content.replace(
        "window.location.href = '/results';",
        &format!("window.location.href = '/results?odds_ratio_threshold={}&n_protein_threshold={}';", odds_ratio, n_protein)
    );
    
    HttpResponse::Ok()
        .content_type("text/html")
        .body(html_with_params)
}

async fn results(query: web::Query<QueryParams>, _state: web::Data<Arc<AppState>>) -> HttpResponse {
    info!("Processing analysis results request");
    
    // Log the selected background if present
    if let Some(background) = &query.selected_background {
        debug!("Selected background: {}", background);
    }
    
    // Read the original results file
    let results_path = Path::new("results/results.txt");
    if let Ok(content) = fs::read_to_string(results_path) {
        let mut biological_process_lines = Vec::new();
        let mut molecular_function_lines = Vec::new();
        let mut cellular_component_lines = Vec::new();
        let mut is_header = true;
        let mut header_line = String::new();

        for line in content.lines() {
            if is_header {
                header_line = line.to_string();
                is_header = false;
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 {
                // Extract namespace directly from the results file (3rd column)
                if let Some(namespace) = fields.get(2) {
                    match namespace.trim().to_lowercase().as_str() {
                        "biological_process" => biological_process_lines.push(line.to_string()),
                        "molecular_function" => molecular_function_lines.push(line.to_string()),
                        "cellular_component" => cellular_component_lines.push(line.to_string()),
                        _ => {}
                    }
                }
            }
        }

        // Write namespace-specific files
        let write_namespace_file = |lines: Vec<String>, filename: &str| -> Result<(), std::io::Error> {
            let mut content = vec![header_line.clone()];
            content.extend(lines);
            fs::write(format!("results/{}", filename), content.join("\n"))
        };

        if let Err(e) = write_namespace_file(biological_process_lines, "biological_process_results.txt") {
            error!("Failed to write biological process results: {}", e);
            return HttpResponse::InternalServerError().body("Failed to write biological process results file");
        }

        if let Err(e) = write_namespace_file(molecular_function_lines, "molecular_function_results.txt") {
            error!("Failed to write molecular function results: {}", e);
            return HttpResponse::InternalServerError().body("Failed to write molecular function results file");
        }

        if let Err(e) = write_namespace_file(cellular_component_lines, "cellular_component_results.txt") {
            error!("Failed to write cellular component results: {}", e);
            return HttpResponse::InternalServerError().body("Failed to write cellular component results file");
        }
        
        let html_content = include_str!("templates/results.html");
        HttpResponse::Ok()
            .content_type("text/html")
            .body(html_content)
    } else {
        error!("Failed to read original results file");
        
        // Create a custom HTML page with an error popup
        let error_html = r#"<!DOCTYPE html>
<html>
<head>
    <title>TaxaGO - Error</title>
    <style>
        :root {
            --primary-color: #4f46e5;
            --primary-hover: #4338ca;
            --background-color: #f8fafc;
            --card-background: #ffffff;
            --text-primary: #1e293b;
            --text-secondary: #64748b;
            --border-color: #e2e8f0;
            --error-color: #ef4444;
        }

        body {
            margin: 0;
            padding: 20px;
            font-family: system-ui, -apple-system, sans-serif;
            background-color: var(--background-color);
            color: var(--text-primary);
            line-height: 1.5;
            display: flex;
            justify-content: center;
            align-items: center;
            min-height: 100vh;
        }

        .error-container {
            text-align: center;
            padding: 2rem;
            background: var(--card-background);
            border-radius: 1rem;
            box-shadow: 0 4px 6px -1px rgb(0 0 0 / 0.1);
            max-width: 600px;
            width: 100%;
        }

        .error-icon {
            color: var(--error-color);
            font-size: 3rem;
            margin-bottom: 1rem;
        }

        .error-title {
            font-size: 1.5rem;
            font-weight: 600;
            margin-bottom: 1rem;
            color: var(--text-primary);
        }

        .error-message {
            font-size: 1rem;
            margin-bottom: 2rem;
            color: var(--text-secondary);
        }

        .new-analysis-btn {
            background-color: var(--primary-color);
            color: white;
            border: none;
            padding: 0.75rem 1.5rem;
            font-size: 1rem;
            font-weight: 600;
            border-radius: 0.5rem;
            cursor: pointer;
            transition: all 0.3s ease;
        }

        .new-analysis-btn:hover {
            background-color: var(--primary-hover);
            transform: translateY(-2px);
        }
    </style>
</head>
<body>
    <div class="error-container">
        <div class="error-title">Analysis Error</div>
        <div class="error-message">The background and study population do not match. Please try again with a compatible background.</div>
        <button class="new-analysis-btn" id="newAnalysisBtn">New Analysis</button>
    </div>

    <script>
        document.getElementById('newAnalysisBtn').addEventListener('click', async () => {
            try {
                // Clear results first
                await fetch('/clear-results', {
                    method: 'POST'
                });
                
                // Also clear data folder
                await fetch('/clear-data-folder', {
                    method: 'POST'
                });
                
                // Redirect to index page
                window.location.href = '/';
            } catch (error) {
                console.error('Error clearing data:', error);
                // Redirect anyway
                window.location.href = '/';
            }
        });
    </script>
</body>
</html>"#;
        
        HttpResponse::Ok()
            .content_type("text/html")
            .body(error_html)
    }
}

async fn get_results(_state: web::Data<Arc<AppState>>) -> HttpResponse {
    debug!("Handling results data request");
    let results_path = Path::new("results/results.txt");
    match fs::read_to_string(results_path) {
        Ok(content) => {
            debug!("Successfully retrieved results data");
            HttpResponse::Ok()
                .content_type("text/plain")
                .body(content)
        },
        Err(e) => {
            warn!("Failed to read results file: {}", e);
            HttpResponse::NotFound()
                .json(json!({
                    "error": "Results file not found",
                    "message": "The background and study population do not match. Please try again with a compatible background."
                }))
        }
    }
}

async fn serve_plot(name: web::Path<String>, _state: web::Data<Arc<AppState>>) -> HttpResponse {
    debug!("Serving plot: {}", name);
    let plot_path = format!("results/plots/{}", name);
    
    // Check if the file exists
    if !Path::new(&plot_path).exists() {
        warn!("Plot file not found: {}", name);
        
        // Determine content type based on file extension for proper error handling
        let content_type = if name.ends_with(".html") {
            // For HTML plots, return an HTML error page
            return HttpResponse::NotFound()
                .content_type("text/html")
                .body(format!(
                    r#"<!DOCTYPE html>
                    <html>
                    <head>
                        <title>Error: Plot Not Found</title>
                        <style>
                            body {{ font-family: system-ui, sans-serif; text-align: center; padding: 2rem; }}
                            .error-message {{ color: #64748b; }}
                        </style>
                    </head>
                    <body class="plot-error-page">
                        <h1>Plot Not Found</h1>
                        <p class="error-message">The requested plot "{}" could not be found.</p>
                    </body>
                    </html>"#,
                    name
                ));
        } else {
            "application/octet-stream"
        };
        
        return HttpResponse::NotFound()
            .content_type(content_type)
            .body(format!("Plot {} not found", name));
    }
    
    // Determine content type based on file extension
    let content_type = if name.ends_with(".html") {
        "text/html"
    } else {
        "application/octet-stream"
    };
    
    // Read the file as bytes (works for both text and binary files)
    match fs::read(&plot_path) {
        Ok(content) => {
            debug!("Successfully retrieved plot data");
            HttpResponse::Ok()
                .content_type(content_type)
                .body(content)
        },
        Err(e) => {
            warn!("Failed to read plot file {}: {}", name, e);
            
            // Return appropriate error based on file type
            if content_type == "text/html" {
                HttpResponse::InternalServerError()
                    .content_type("text/html")
                    .body(format!(
                        r#"<!DOCTYPE html>
                        <html>
                        <head>
                            <title>Error: Failed to Read Plot</title>
                            <style>
                                body {{ font-family: system-ui, sans-serif; text-align: center; padding: 2rem; }}
                                .error-message {{ color: #64748b; }}
                            </style>
                        </head>
                        <body>
                            <h1>Error Reading Plot</h1>
                            <p class="error-message">Failed to read plot "{}": {}</p>
                        </body>
                        </html>"#,
                        name, e
                    ))
            } else {
                HttpResponse::InternalServerError()
                    .content_type(content_type)
                    .body(format!("Failed to read plot {}: {}", name, e))
            }
        }
    }
}

async fn download_results(format: web::Path<String>, _state: web::Data<Arc<AppState>>) -> HttpResponse {
    info!("Processing download request for format: {}", format);
    let download_manager = download::DownloadManager::new("results/results.txt".to_string());
    
    match download_manager.generate_file_content(&format) {
        Ok((filename, content)) => {
            let content_type = match format.as_str() {
                "csv" => "text/csv",
                "tsv" => "text/tab-separated-values",
                "json" => "application/json",
                _ => "application/octet-stream",
            };
            
            debug!("Successfully generated {} file for download", format);
            HttpResponse::Ok()
                .content_type(content_type)
                .append_header(("Content-Disposition", format!("attachment; filename=\"{}\"", filename)))
                .body(content)
        },
        Err(e) => {
            error!("Failed to generate download file: {}", e);
            HttpResponse::InternalServerError()
                .body(format!("Error generating file: {}", e))
        }
    }
}

async fn get_namespace_results(namespace: web::Path<String>) -> Result<HttpResponse> {
    let file_path = match namespace.as_str() {
        "biological-process" => "results/biological_process_results.txt",
        "molecular-function" => "results/molecular_function_results.txt",
        "cellular-component" => "results/cellular_component_results.txt",
        _ => return Ok(HttpResponse::BadRequest().body("Invalid namespace")),
    };

    if let Ok(content) = fs::read_to_string(file_path) {
        Ok(HttpResponse::Ok()
            .content_type("text/plain")
            .body(content))
    } else {
        error!("Failed to read namespace file: {}", file_path);
        Ok(HttpResponse::NotFound().body("Results file not found"))
    }
}

#[derive(serde::Serialize)]
struct Species {
    taxon_id: String,
    name: String,
}

async fn get_species() -> Result<HttpResponse> {
    // Return an empty list since full_lineage.txt is no longer used
    let species_list: Vec<Species> = Vec::new();
    Ok(HttpResponse::Ok().json(species_list))
}

async fn get_lineage_data() -> Result<HttpResponse> {
    // Return an empty response since full_lineage.txt is no longer used
    Ok(HttpResponse::Ok()
        .content_type("text/plain")
        .body(""))
}

#[post("/update-categories")]
async fn update_categories(data: web::Json<CategoryUpdate>) -> impl Responder {
    debug!("Received category update: {:?}", data);
    let file_path = "data/user_selected_categories.txt";
    
    if data.is_user_action {
        // Create data directory if it doesn't exist
        if let Err(e) = std::fs::create_dir_all("data") {
            error!("Failed to create data directory: {}", e);
            return HttpResponse::InternalServerError().body("Failed to create data directory");
        }

        let category_line = format!("{}\t{}\t{}", data.superkingdom, data.level, data.category);
        debug!("Category line: {}", category_line);
        
        if data.is_checked {
            // Add category
            match OpenOptions::new()
                .create(true)
                .append(true)
                .open(file_path) 
            {
                Ok(mut file) => {
                    if let Err(e) = writeln!(file, "{}", category_line) {
                        error!("Failed to write to file: {}", e);
                        return HttpResponse::InternalServerError().body("Failed to write to file");
                    }
                    debug!("Successfully added category to file");
                },
                Err(e) => {
                    error!("Failed to open file for writing: {}", e);
                    return HttpResponse::InternalServerError().body("Failed to open file for writing");
                }
            }
        } else {
            // Remove category
            if Path::new(file_path).exists() {
                match File::open(file_path) {
                    Ok(file) => {
                        let reader = BufReader::new(file);
                        let lines: Vec<String> = reader
                            .lines()
                            .filter_map(Result::ok)
                            .filter(|line| line != &category_line)
                            .collect();

                        match File::create(file_path) {
                            Ok(mut file) => {
                                for line in lines {
                                    if let Err(e) = writeln!(file, "{}", line) {
                                        error!("Failed to write to file: {}", e);
                                        return HttpResponse::InternalServerError().body("Failed to write to file");
                                    }
                                }
                                debug!("Successfully removed category from file");
                            },
                            Err(e) => {
                                error!("Failed to create file: {}", e);
                                return HttpResponse::InternalServerError().body("Failed to create file");
                            }
                        }
                    },
                    Err(e) => {
                        error!("Failed to open file for reading: {}", e);
                        return HttpResponse::InternalServerError().body("Failed to open file for reading");
                    }
                }
            }
        }
    }
    
    HttpResponse::Ok().finish()
}

async fn get_saved_categories() -> impl Responder {
    let file_path = "data/user_selected_categories.txt";
    
    if !Path::new(file_path).exists() {
        return HttpResponse::Ok().json(Vec::<SavedCategory>::new());
    }

    match File::open(file_path) {
        Ok(file) => {
            let reader = BufReader::new(file);
            let categories: Vec<SavedCategory> = reader
                .lines()
                .filter_map(Result::ok)
                .filter_map(|line| {
                    let parts: Vec<&str> = line.split('\t').collect();
                    if parts.len() >= 3 {
                        Some(SavedCategory {
                            superkingdom: parts[0].to_string(),
                            level: parts[1].to_string(),
                            category: parts[2].to_string(),
                        })
                    } else {
                        None
                    }
                })
                .collect();

            HttpResponse::Ok().json(categories)
        },
        Err(e) => {
            error!("Failed to read categories file: {}", e);
            HttpResponse::InternalServerError().body("Failed to read categories")
        }
    }
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "snake_case")]
struct AnalysisParams {
    study_pop_path: String,
    obo_path: Option<String>,
    background_path: Option<String>,
    vcv_matrix_path: Option<String>,
    odds_ratio_threshold: f64,
    n_protein_threshold: i32,
    lineage_percentage: f64,
    pm_tolerance: f64,
    statistical_test: String,
    fdr_enabled: bool,
    fdr_method: Option<String>,
    fdr_threshold: Option<f64>,
    propagate_enabled: bool,
    propagate_algorithm: Option<String>,
    combine_enabled: bool,
    taxonomic_level: Option<String>,
    multiple_testing: Option<bool>,
    evidence_codes: Option<String>,
}

async fn execute_analysis(params: web::Json<AnalysisParams>) -> HttpResponse {
    info!("Starting analysis with parameters: {:?}", params);
    
    // Create results directory if it doesn't exist
    if let Err(e) = fs::create_dir_all("results") {
        error!("Failed to create results directory: {}", e);
        return HttpResponse::InternalServerError().body("Failed to create results directory");
    }
    
    // Create plots directory if it doesn't exist
    if let Err(e) = fs::create_dir_all("results/plots") {
        error!("Failed to create plots directory: {}", e);
        return HttpResponse::InternalServerError().body("Failed to create plots directory");
    }
    
    // Create a vector to store all arguments
    let mut args = Vec::new();
    
    // Add all arguments to the vector
    args.push("-s".to_string());
    
    // Use the file from the data directory instead of the original path
    let study_pop_path = format!("data/{}", params.study_pop_path);
    args.push(study_pop_path.clone());
    
    // Print the study population path for debugging
    println!("\n==================================================");
    println!("STUDY POPULATION PATH: {}", study_pop_path);
    println!("Checking if file exists: {}", Path::new(&study_pop_path).exists());
    
    // Print the contents of the data directory
    println!("Contents of data directory:");
    let mut found_obo_file = None;
    let mut found_background_file = None;
    let mut found_vcv_file = None;
    
    if let Ok(entries) = fs::read_dir("data") {
        for entry in entries.filter_map(Result::ok) {
            if let Ok(metadata) = entry.metadata() {
                let file_type = if metadata.is_dir() { "Directory" } else { "File" };
                let path = entry.path();
                println!("  {} - {}", file_type, path.display());
                
                // Check for OBO files
                if let Some(extension) = path.extension() {
                    if extension == "obo" {
                        found_obo_file = Some(path.file_name().unwrap().to_string_lossy().to_string());
                        println!("  Found OBO file: {}", found_obo_file.as_ref().unwrap());
                    }
                }
                
                // Check for background files in the data directory (not in background_population)
                if metadata.is_file() {
                    if let Some(file_name) = path.file_name() {
                        let file_name_str = file_name.to_string_lossy();
                        if file_name_str.ends_with("_background.txt") {
                            found_background_file = Some(file_name_str.to_string());
                            println!("  Found background file in data directory: {}", found_background_file.as_ref().unwrap());
                        } else if file_name_str.ends_with("_dmat.csv") {
                            found_vcv_file = Some(file_name_str.to_string());
                            println!("  Found VCV matrix file: {}", found_vcv_file.as_ref().unwrap());
                        }
                    }
                }
            }
        }
    } else {
        println!("  Failed to read data directory");
    }
    
    // Check for the background_population directory
    let background_population_dir = Path::new("data/background_population");
    let has_background_population = background_population_dir.exists() && background_population_dir.is_dir();
    
    if has_background_population {
        println!("  Found background_population directory");
        
        // Check if the directory contains any background files
        if let Ok(entries) = fs::read_dir(background_population_dir) {
            let background_files: Vec<String> = entries
                .filter_map(Result::ok)
                .filter_map(|entry| {
                    if let Ok(metadata) = entry.metadata() {
                        if metadata.is_file() {
                            if let Some(file_name) = entry.path().file_name() {
                                let file_name_str = file_name.to_string_lossy().to_string();
                                if file_name_str.ends_with("_background.txt") {
                                    return Some(file_name_str);
                                }
                            }
                        }
                    }
                    None
                })
                .collect();
            
            if !background_files.is_empty() {
                println!("  Found background files in background_population directory:");
                for file in &background_files {
                    println!("    {}", file);
                }
            } else {
                println!("  No background files found in background_population directory");
            }
        }
    } else {
        println!("  No background_population directory found");
    }
    
    println!("==================================================\n");

    // Only add the OBO path if it's explicitly provided in params
    if let Some(obo_path) = &params.obo_path {
        if !obo_path.is_empty() {
            println!("Custom OBO file specified: data/{}", obo_path);
            args.push("-o".to_string());
            args.push(format!("data/{}", obo_path));
        }
    } else if let Some(obo_file) = found_obo_file {
        // Use the OBO file found in the data directory
        println!("Found OBO file in data directory: {}", obo_file);
        args.push("-o".to_string());
        args.push(format!("data/{}", obo_file));
    }
    // If no OBO file is provided or found, don't add the -o option at all

    // Check for background population folder and add it to the command if it exists
    let background_population_dir = Path::new("data/background_population");
    if background_population_dir.exists() && background_population_dir.is_dir() {
        // Check if the directory contains any background files
        let has_background_files = if let Ok(entries) = fs::read_dir(background_population_dir) {
            entries.filter_map(Result::ok)
                .any(|entry| {
                    if let Ok(metadata) = entry.metadata() {
                        if metadata.is_file() {
                            if let Some(file_name) = entry.path().file_name() {
                                let file_name_str = file_name.to_string_lossy();
                                return file_name_str.ends_with("_background.txt");
                            }
                        }
                    }
                    false
                })
        } else {
            false
        };

        if has_background_files {
            args.push("-b".to_string());
            args.push("data/background_population".to_string());
            println!("Using background_population directory for background files");
        }
    } else if let Some(background_path) = &params.background_path {
        // Fallback to the specified background path if provided
        if !background_path.is_empty() {
            args.push("-b".to_string());
            args.push(format!("data/{}", background_path));
            println!("Using specified background file: {}", background_path);
        }
    } else if let Some(background_file) = found_background_file {
        // Fallback to any found background file in the data directory
        args.push("-b".to_string());
        args.push(format!("data/{}", background_file));
        println!("Using background file found in data directory: {}", background_file);
    }

    // Add VCV matrix path if provided or found
    if let Some(vcv_path) = &params.vcv_matrix_path {
        if !vcv_path.is_empty() {
            println!("Custom VCV matrix file specified: data/{}", vcv_path);
            args.push("--vcv-matrix".to_string());
            args.push(format!("data/{}", vcv_path));
        }
    } else if let Some(vcv_file) = found_vcv_file {
        // Use the VCV matrix file found in the data directory
        println!("Found VCV matrix file in data directory: {}", vcv_file);
        args.push("--vcv-matrix".to_string());
        args.push(format!("data/{}", vcv_file));
    }

    args.push("-d".to_string());
    args.push("results".to_string());
    
    // Add optional parameters based on configuration
    args.push("-t".to_string());
    args.push(params.statistical_test.clone());
    
    let n_protein_str = params.n_protein_threshold.to_string();
    args.push("-m".to_string());
    args.push(n_protein_str);
    
    let score_str = params.odds_ratio_threshold.to_string();
    args.push("-r".to_string());
    args.push(score_str);
    
    // Add evidence code filter if provided
    if let Some(evidence_codes) = &params.evidence_codes {
        let trimmed = evidence_codes.trim();
        if !trimmed.is_empty() {
            args.push("-e".to_string());
            args.push(trimmed.to_string());
        }
    }
    
    if params.propagate_enabled {
        args.push("-p".to_string());
        // Add the algorithm value (classic, elim, weight) after -p
        if let Some(algo) = params.propagate_algorithm.as_ref() {
            args.push(algo.clone());
        } else {
            args.push("classic".to_string()); // fallback default
        }
        // Save the propagation setting to a marker file for other commands to use
        let marker_file = Path::new("data/.propagation_enabled");
        if let Err(e) = fs::write(marker_file, "1") {
            error!("Failed to create propagation marker file: {}", e);
            // Continue anyway, this is just a helper file
        }
    } else {
        // Remove the marker file if it exists
        let marker_file = Path::new("data/.propagation_enabled");
        if marker_file.exists() {
            if let Err(e) = fs::remove_file(marker_file) {
                error!("Failed to remove propagation marker file: {}", e);
                // Continue anyway, this is just a helper file
            }
        }
    }
    
    if params.fdr_enabled {
        if let Some(method) = &params.fdr_method {
            args.push("-c".to_string());
            args.push(method.clone());
        }
    }
    
    // Handle significance threshold independently
    if let Some(threshold) = params.fdr_threshold {
        // Only add if different from default (0.05)
        if (threshold - 0.05).abs() > 0.0001 {
            args.push("-a".to_string());
            args.push(threshold.to_string());
        }
    }
    
    if params.combine_enabled {
        if let Some(level) = &params.taxonomic_level {
            args.push("-g".to_string());
            args.push(level.clone());            
            args.push("-l".to_string());
            args.push(params.lineage_percentage.to_string());
            args.push("--permutations".to_string());
            args.push(params.pm_tolerance.to_string());
        } else {
            // If combination is enabled but no level is specified, log a warning
            warn!("Combination enabled but no taxonomic level specified");
        }
    }
    
    // Always add --save-plots to the command arguments
    args.push("--save-plots".to_string());
    
    // Add this before the Command execution
    match std::env::current_dir() {
        Ok(dir) => info!("Current working directory: {:?}", dir),
        Err(e) => error!("Failed to get current directory: {}", e),
    }
    
    // Create the full command string for display
    let command_str = format!("taxago {}", args.join(" "));
    
    // Print the command to terminal in a highly visible way
    println!("\n==================================================");
    println!("EXECUTING COMMAND:");
    println!("{}", command_str);
    println!("==================================================\n");
    
    info!("Executing command: taxago with args: {:?}", args);
    
    // Determine if this is a single analysis by checking the study population file
    // and the multiple_testing parameter from the frontend
    let is_multiple_testing = params.multiple_testing.unwrap_or(false);
    let is_single_analysis = !is_multiple_testing && !params.study_pop_path.contains(",") && !params.combine_enabled;
    
    // Log the determination for debugging
    debug!("Study population path: {}", params.study_pop_path);
    debug!("Contains comma: {}", params.study_pop_path.contains(","));
    debug!("Combine enabled: {}", params.combine_enabled);
    debug!("Multiple testing (from frontend): {:?}", params.multiple_testing);
    debug!("Is single analysis: {}", is_single_analysis);
    
    // Create a marker file to indicate if this is single testing
    if is_single_analysis {
        if let Err(e) = fs::write("results/single_testing_marker", "") {
            warn!("Failed to create single testing marker file: {}", e);
        } else {
            debug!("Created single testing marker file");
        }
    } else {
        // Remove the marker file if it exists
        if let Err(e) = fs::remove_file("results/single_testing_marker") {
            if e.kind() != std::io::ErrorKind::NotFound {
                warn!("Failed to remove single testing marker file: {}", e);
            }
        }
    }
    
    match Command::new("taxago")
        .args(&args)
        .output()
    {
        Ok(output) => {
            if output.status.success() {
                info!("Analysis completed successfully");
                
                // For single analysis (not combined), copy the results file automatically
                if !params.combine_enabled && !params.study_pop_path.contains(",") {
                    // Try to find the single result file and copy it
                    if let Ok(entries) = fs::read_dir("results/single_taxon_results") {
                        let files: Vec<_> = entries
                            .filter_map(Result::ok)
                            .filter(|entry| {
                                if let Ok(file_type) = entry.file_type() {
                                    file_type.is_file()
                                } else {
                                    false
                                }
                            })
                            .collect();
                        
                        // If there's exactly one file, copy it to results.txt
                        if files.len() == 1 {
                            let source_path = files[0].path();
                            let dest_path = "results/results.txt";
                            
                            // Read the source file and replace headers instead of just copying
                            match fs::read_to_string(&source_path) {
                                Ok(content) => {
                                    // Get the original header
                                    let original_header = content.lines().next().unwrap_or("");
                                    
                                    // Store the original filename and header in a metadata file for later use by the download manager
                                    if let Some(filename) = source_path.file_name().and_then(|f| f.to_str()) {
                                        let metadata_path = "results/original_file_info.txt";
                                        let metadata_content = format!("{}\n{}", filename, original_header);
                                        if let Err(e) = fs::write(metadata_path, metadata_content) {
                                            error!("Failed to write metadata file: {}", e);
                                            // Continue anyway, this is not critical
                                        }
                                    }
                                    
                                    // Replace the header
                                    let mut lines: Vec<&str> = content.lines().collect();
                                    if !lines.is_empty() {
                                        // Replace the first line (header) with the new header
                                        lines[0] = "GO ID\tName\tNamespace\tlog(Odds Ratio)\tStat. Sig.";
                                    }
                                    let modified_content = lines.join("\n");
                                    
                                    // Write the modified content to the destination file
                                    if let Err(e) = fs::write(dest_path, modified_content) {
                                        error!("Failed to write modified single result file: {}", e);
                                    } else {
                                        info!("Automatically copied and modified single result file to results.txt");
                                    }
                                },
                                Err(e) => {
                                    error!("Failed to read single result file: {}", e);
                                    // Fall back to simple copy if reading fails
                                    if let Err(e) = fs::copy(&source_path, dest_path) {
                                        error!("Failed to copy single result file: {}", e);
                                    }
                                }
                            }
                        }
                    }
                }
                // For combined analysis, check if there's a single combined file to copy
                else if params.combine_enabled {
                    // Try to find combined result files
                    if let Ok(entries) = fs::read_dir("results/combined_taxonomy_results") {
                        let files: Vec<_> = entries
                            .filter_map(Result::ok)
                            .filter(|entry| {
                                if let Ok(file_type) = entry.file_type() {
                                    file_type.is_file()
                                } else {
                                    false
                                }
                            })
                            .collect();
                        
                        // If there's exactly one file, copy it to results.txt
                        if files.len() == 1 {
                            let source_path = files[0].path();
                            let dest_path = "results/results.txt";
                            
                            // Read the source file and replace headers instead of just copying
                            match fs::read_to_string(&source_path) {
                                Ok(content) => {
                                    // Get the original header
                                    let original_header = content.lines().next().unwrap_or("");
                                    
                                    // Store the original filename and header in a metadata file for later use by the download manager
                                    if let Some(filename) = source_path.file_name().and_then(|f| f.to_str()) {
                                        let metadata_path = "results/original_file_info.txt";
                                        let metadata_content = format!("{}\n{}", filename, original_header);
                                        if let Err(e) = fs::write(metadata_path, metadata_content) {
                                            error!("Failed to write metadata file: {}", e);
                                            // Continue anyway, this is not critical
                                        }
                                    }
                                    
                                    // Replace the header
                                    let mut lines: Vec<&str> = content.lines().collect();
                                    if !lines.is_empty() {
                                        // Replace the first line (header) with the new header for combined results
                                        lines[0] = "GO ID\tName\tNamespace\tlog(Odds Ratio)\tStat. Sig.\tHetergnt\tSpecies %\tN w/\tN in tax";
                                    }
                                    let modified_content = lines.join("\n");
                                    
                                    // Write the modified content to the destination file
                                    if let Err(e) = fs::write(dest_path, modified_content) {
                                        error!("Failed to write modified combined result file: {}", e);
                                    } else {
                                        info!("Automatically copied and modified combined result file to results.txt");
                                    }
                                },
                                Err(e) => {
                                    error!("Failed to read combined result file: {}", e);
                                    // Fall back to simple copy if reading fails
                                    if let Err(e) = fs::copy(&source_path, dest_path) {
                                        error!("Failed to copy combined result file: {}", e);
                                    }
                                }
                            }
                        }
                    }
                }
                
                let stdout = String::from_utf8_lossy(&output.stdout).to_string();
                let stderr = String::from_utf8_lossy(&output.stderr).to_string();
                info!("Command executed successfully. Stdout: {}, Stderr: {}", stdout, stderr);
                
                // After the analysis is complete, check if we need to auto-copy single test results
                if !params.combine_enabled {
                    // Wait a moment for files to be written
                    tokio::time::sleep(tokio::time::Duration::from_secs(1)).await;
                    
                    // Auto-copy single test results if there's only one
                    if let Err(e) = auto_copy_single_test_results().await {
                        error!("Failed to auto-copy single test results: {}", e);
                        // Continue anyway, this is not critical
                    }
                }
                
                HttpResponse::Ok().json(json!({
                    "success": true,
                    "is_single_analysis": is_single_analysis,
                    "redirect": format!("/loading?odds_ratio_threshold={}&n_protein_threshold={}", 
                        params.odds_ratio_threshold, params.n_protein_threshold),
                    "result": if stdout.is_empty() { stderr } else { stdout },
                    "command": command_str
                }))
            } else {
                let stderr = String::from_utf8_lossy(&output.stderr).to_string();
                error!("Analysis failed: {}", stderr);
                
                // Print the error and command to terminal in a highly visible way
                println!("\n==================================================");
                println!("COMMAND FAILED:");
                println!("taxago {}", args.join(" "));
                println!("ERROR: {}", stderr);
                println!("==================================================\n");
                
                HttpResponse::InternalServerError().json(json!({
                    "success": false,
                    "error": stderr,
                    "command": command_str
                }))
            }
        },
        Err(e) => {
            error!("Failed to execute analysis command: {}", e);
            
            // Print the error and command to terminal in a highly visible way
            println!("\n==================================================");
            println!("COMMAND EXECUTION FAILED:");
            println!("taxago {}", args.join(" "));
            println!("ERROR: {}", e);
            println!("==================================================\n");
            
            HttpResponse::InternalServerError().json(json!({
                "success": false,
                "error": e.to_string(),
                "command": command_str
            }))
        }
    }
}

// Add this new endpoint handler
async fn serve_pdf() -> HttpResponse {
    info!("Attempting to serve PDF");
    let path = Path::new("results/ontology_graph.pdf");
    
    // Get current working directory
    match std::env::current_dir() {
        Ok(dir) => info!("Current working directory: {:?}", dir),
        Err(e) => error!("Failed to get current directory: {}", e),
    }
    
    info!("Looking for PDF at absolute path: {:?}", path.canonicalize().unwrap_or_default());
    
    if !path.exists() {
        error!("PDF file not found at {:?}", path);
        // List contents of results directory
        if let Ok(entries) = fs::read_dir("results") {
            info!("Contents of results directory:");
            for entry in entries {
                if let Ok(entry) = entry {
                    info!("  {:?}", entry.path());
                }
            }
        }
        return HttpResponse::NotFound().body("PDF file not found");
    }

    match fs::read(path) {
        Ok(contents) => {
            info!("Successfully read PDF file of size: {} bytes", contents.len());
            HttpResponse::Ok()
                .content_type("application/pdf")
                .append_header(("Content-Disposition", "inline"))
                .body(contents)
        },
        Err(e) => {
            error!("Failed to read PDF file: {}", e);
            HttpResponse::InternalServerError().body(format!("Failed to read PDF: {}", e))
        }
    }
}

// Add this function to list files in the combined_taxonomy_results directory
async fn get_combined_files() -> HttpResponse {
    let dir_path = "results/combined_taxonomy_results";
    
    // Create directory if it doesn't exist
    if !Path::new(dir_path).exists() {
        if let Err(e) = fs::create_dir_all(dir_path) {
            error!("Failed to create combined results directory: {}", e);
            return HttpResponse::InternalServerError().json(Vec::<String>::new());
        }
    }
    
    match fs::read_dir(dir_path) {
        Ok(entries) => {
            let files: Vec<String> = entries
                .filter_map(Result::ok)
                .filter(|entry| {
                    if let Ok(file_type) = entry.file_type() {
                        file_type.is_file()
                    } else {
                        false
                    }
                })
                .filter_map(|entry| {
                    entry.file_name().to_str().map(String::from)
                })
                .collect();
            
            HttpResponse::Ok().json(files)
        },
        Err(e) => {
            error!("Failed to read combined results directory: {}", e);
            HttpResponse::InternalServerError().json(Vec::<String>::new())
        }
    }
}

// Add this function to list files in the single_taxon_results directory
async fn get_single_files() -> HttpResponse {
    let dir_path = "results/single_taxon_results";
    
    // Create directory if it doesn't exist
    if !Path::new(dir_path).exists() {
        if let Err(e) = fs::create_dir_all(dir_path) {
            error!("Failed to create single results directory: {}", e);
            return HttpResponse::InternalServerError().json(Vec::<String>::new());
        }
    }
    
    match fs::read_dir(dir_path) {
        Ok(entries) => {
            let files: Vec<String> = entries
                .filter_map(Result::ok)
                .filter(|entry| {
                    if let Ok(file_type) = entry.file_type() {
                        file_type.is_file()
                    } else {
                        false
                    }
                })
                .filter_map(|entry| {
                    entry.file_name().to_str().map(String::from)
                })
                .collect();
            
            HttpResponse::Ok().json(files)
        },
        Err(e) => {
            error!("Failed to read single results directory: {}", e);
            HttpResponse::InternalServerError().json(Vec::<String>::new())
        }
    }
}

// Add this struct for deserializing the request
#[derive(Deserialize)]
struct CopyFileRequest {
    #[serde(rename = "type")]
    file_type: String,
    filename: String,
}

// Add this function to copy plots from single_taxon_results to results directory
async fn copy_plots(species_name: &str) -> Result<(), std::io::Error> {
    info!("Copying plots for species: {}", species_name);
    
    // Create plots directory if it doesn't exist
    let plots_dir = "results/plots";
    if let Err(e) = fs::create_dir_all(plots_dir) {
        error!("Failed to create plots directory: {}", e);
        return Err(e);
    }
    
    // Remove _GOEA_results.txt from species name if present
    let clean_species_name = species_name.replace("_GOEA_results.txt", "");
    
    // Use underscores for the new naming convention
    let species_name_with_underscores = clean_species_name.clone();

    info!("Looking for plots with species name: {}", species_name_with_underscores);
    
    // Define the namespaces and their standardized plot names
    let namespaces = [
        ("Biological_Process", "bp"),
        ("Molecular_Function", "mf"),
        ("Cellular_Component", "cc")
    ];
    
    // Define which plots to copy and their new names (new convention)
    let plot_types = [
        ("network_plot.html", "_network.html"),
        ("bubble_plot.html", "_bubble.html"),
        ("bar_plot.html", "_bar.html")
    ];
    
    // Copy plots for each namespace
    for (namespace, prefix) in namespaces.iter() {
        let source_dir = format!("results/single_taxon_results/plots/{}", namespace);
        
        // For each plot type
        for (source_suffix, dest_suffix) in plot_types.iter() {
            let source_path = format!("{}/{}_{}", source_dir, species_name_with_underscores, source_suffix);
            let dest_path = format!("{}/{}{}", plots_dir, prefix, dest_suffix);
            
            // Copy the file if it exists
            if Path::new(&source_path).exists() {
                info!("Copying {} to {}", source_path, dest_path);
                if let Err(e) = fs::copy(&source_path, &dest_path) {
                    error!("Failed to copy plot {}: {}", source_path, e);
                }
            } else {
                warn!("Plot not found: {}", source_path);
            }
        }
    }
    
    Ok(())
}

// Add this endpoint handler
#[post("/copy-results-file")]
async fn copy_results_file(data: web::Json<CopyFileRequest>) -> impl Responder {
    info!("Copying results file: {:?} - {}", data.file_type, data.filename);
    
    let source_path = if data.file_type == "combined_file" {
        format!("results/combined_taxonomy_results/{}", data.filename)
    } else {
        format!("results/single_taxon_results/{}", data.filename)
    };
    
    let dest_path = "results/results.txt";
    
    // Ensure the results directory exists
    if let Err(e) = fs::create_dir_all("results") {
        error!("Failed to create results directory: {}", e);
        return HttpResponse::InternalServerError().body("Failed to create results directory");
    }
    
    // Read the source file
    match fs::read_to_string(&source_path) {
        Ok(content) => {
            // Get the original header
            let original_header = content.lines().next().unwrap_or("");
            
            // Store the original filename and header in a metadata file for later use by the download manager
            let metadata_path = "results/original_file_info.txt";
            let metadata_content = format!("{}\n{}", data.filename, original_header);
            if let Err(e) = fs::write(metadata_path, metadata_content) {
                error!("Failed to write metadata file: {}", e);
                // Continue anyway, this is not critical
            }
            
            let modified_content = if data.file_type == "single_file" {
                // For single_taxon_results files, replace the header
                let mut lines: Vec<&str> = content.lines().collect();
                if !lines.is_empty() {
                    // Replace the first line (header) with the new header
                    lines[0] = "GO ID\tName\tNamespace\tlog(Odds Ratio)\tStat. Sig.";
                }
                lines.join("\n")
            } else if data.file_type == "combined_file" {
                // For combined_taxonomy_results files, replace the header
                let mut lines: Vec<&str> = content.lines().collect();
                if !lines.is_empty() {
                    // Replace the first line (header) with the new header for combined results
                    lines[0] = "GO ID\tName\tNamespace\tlog(Odds Ratio)\tStat. Sig.\tHetergnt\tSpecies %\tN w/\tN in tax";
                }
                lines.join("\n")
            } else {
                // For any other type, keep the original content
                content
            };
            
            // Write the modified content to the destination file
            match fs::write(dest_path, modified_content) {
                Ok(_) => {
                    info!("Successfully copied and modified {} to {}", source_path, dest_path);
                    
                    // Copy the appropriate plots based on the file type
                    if data.file_type == "single_file" {
                        if let Err(e) = copy_plots(&data.filename).await {
                            error!("Failed to copy single plots: {}", e);
                            // Continue anyway, this is not critical
                        }
                    } else if data.file_type == "combined_file" {
                        if let Err(e) = copy_combined_plots(&data.filename).await {
                            error!("Failed to copy combined plots: {}", e);
                            // Continue anyway, this is not critical
                        }
                    }
                    
                    HttpResponse::Ok().finish()
                },
                Err(e) => {
                    error!("Failed to write to {}: {}", dest_path, e);
                    HttpResponse::InternalServerError().body(format!("Failed to write file: {}", e))
                }
            }
        },
        Err(e) => {
            error!("Failed to read file from {}: {}", source_path, e);
            HttpResponse::InternalServerError().body(format!("Failed to read file: {}", e))
        }
    }
}

// Add this function to clear results
async fn clear_results() -> HttpResponse {
    info!("Clearing results directory");
    
    // Create a single cleanup command for the results directory
    let cleanup_script = "rm -rf results";

    // Execute cleanup as a single shell command
    let output = std::process::Command::new("sh")
        .arg("-c")
        .arg(cleanup_script)
        .output();

    match output {
        Ok(output) => {
            if output.status.success() {
                info!("✓ Successfully removed results directory");
            } else {
                error!("✗ Cleanup failed: {}", String::from_utf8_lossy(&output.stderr));
                // Fallback to Rust's fs::remove_dir_all
                if let Err(e) = fs::remove_dir_all("results") {
                    error!("Fallback cleanup also failed: {}", e);
                    return HttpResponse::InternalServerError().body("Failed to clear results directory");
                }
            }
        },
        Err(e) => {
            error!("✗ Failed to execute cleanup command: {}", e);
            // Fallback to Rust's fs::remove_dir_all
            if let Err(e) = fs::remove_dir_all("results") {
                error!("Fallback cleanup also failed: {}", e);
                return HttpResponse::InternalServerError().body("Failed to clear results directory");
            }
        }
    }

    // Verify final state
    let results_exists = Path::new("results").exists();
    info!("Results directory exists after cleanup: {}", results_exists);

    if !results_exists {
        info!("✓ Results cleanup successful!");
    } else {
        error!("✗ Results directory still exists!");
        return HttpResponse::InternalServerError().body("Failed to fully clear results directory");
    }

    HttpResponse::Ok().finish()
}

#[get("/check-testing-type")]
async fn check_testing_type() -> impl Responder {
    // First check if the marker file exists - this is the most reliable indicator
    let marker_exists = Path::new("results/single_testing_marker").exists();
    
    if marker_exists {
        debug!("Single testing marker file found - this is definitely single testing");
        return web::Json(json!({
            "is_single_testing": true
        }));
    }
    
    // If no marker file, check the directory structure
    let combined_dir = Path::new("results/combined_taxonomy_results");
    let has_combined_files = combined_dir.exists() && 
        fs::read_dir(combined_dir).map(|entries| entries.count() > 0).unwrap_or(false);
    
    let single_dir = Path::new("results/single_taxon_results");
    let single_files_count = if single_dir.exists() {
        let count = fs::read_dir(single_dir)
            .map(|entries| {
                let files: Vec<_> = entries
                    .filter_map(Result::ok)
                    .filter(|entry| {
                        if let Ok(file_type) = entry.file_type() {
                            file_type.is_file()
                        } else {
                            false
                        }
                    })
                    .collect();
                
                // Log the actual files found for debugging
                for (i, file) in files.iter().enumerate() {
                    if let Some(name) = file.file_name().to_str() {
                        debug!("Single taxon file {}: {}", i+1, name);
                    }
                }
                
                files.len()
            })
            .unwrap_or(0);
        
        debug!("Found {} files in single_taxon_results directory", count);
        count
    } else {
        debug!("single_taxon_results directory does not exist");
        0
    };
    
    debug!("Check testing type - Combined files: {}, Single files count: {}", 
           has_combined_files, single_files_count);
    
    // For multiple testing, we should have either:
    // 1. Files in the combined_taxonomy_results directory, or
    // 2. More than one file in the single_taxon_results directory
    let is_multiple_testing = has_combined_files || single_files_count > 1;
    
    // If there's exactly one single test result, automatically copy it and its plots
    if single_files_count == 1 && !has_combined_files {
        debug!("Found exactly one single test result, auto-copying");
        if let Err(e) = auto_copy_single_test_results().await {
            error!("Failed to auto-copy single test results: {}", e);
        }
    }
    
    web::Json(json!({
        "is_single_testing": !is_multiple_testing
    }))
}

#[derive(Deserialize)]
struct AncestorAnalysisRequest {
    go_terms: Vec<String>,
}

#[derive(Deserialize)]
struct SemanticSimilarityRequest {
    go_terms: Vec<String>,
    method: String,
}

#[post("/run-ancestor-analysis")]
async fn run_ancestor_analysis(data: web::Json<AncestorAnalysisRequest>) -> impl Responder {
    info!("Running ancestor analysis for GO terms: {:?}", data.go_terms);
    
    if data.go_terms.len() < 2 {
        return web::Json(json!({
            "success": false,
            "error": "At least two GO terms are required for ancestor analysis",
            "command": ""
        }));
    }
    
    // Delete the previous common_ontology_graph.pdf file if it exists
    let pdf_path = Path::new("results/common_ontology_graph.pdf");
    if pdf_path.exists() {
        info!("Deleting previous ancestor analysis PDF file");
        if let Err(e) = fs::remove_file(pdf_path) {
            error!("Failed to delete previous PDF file: {}", e);
            // Continue with the analysis even if deletion fails
        }
    }

    // Delete the previous common_ontology_graph.mmd file if it exists
    let mmd_path = Path::new("results/common_ontology_graph.mmd");
    if mmd_path.exists() {
        info!("Deleting previous ancestor analysis MMD file");
        if let Err(e) = fs::remove_file(mmd_path) {
            error!("Failed to delete previous MMD file: {}", e);
            // Continue with the analysis even if deletion fails
        }
    }
    
    // Check if any OBO file exists in the data directory
    let mut obo_param = String::new();
    let data_dir = Path::new("data");
    let mut obo_found = false;
    
    if data_dir.exists() && data_dir.is_dir() {
        // Look for .obo files in the data directory
        if let Ok(entries) = fs::read_dir(data_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_file() && 
                   path.extension().map_or(false, |ext| ext == "obo") {
                    if let Some(filename) = path.file_name() {
                        obo_param = format!("-o data/{}", filename.to_string_lossy());
                        obo_found = true;
                        info!("Using OBO file: {}", path.display());
                        break;
                    }
                }
            }
        }
    }
    
    // Build the command with the specified parameters
    let go_terms_arg = data.go_terms.join(",");
    let command = if obo_found {
        format!(
            "common-ancestors {} -t {} -d results/common_",
            obo_param, go_terms_arg
        )
    } else {
        format!(
            "common-ancestors -t {} -d results/common_",
            go_terms_arg
        )
    };
    
    // Print the command to the console for debugging
    println!("\n==================================================");
    println!("ANCESTOR ANALYSIS COMMAND:");
    println!("{}", command);
    println!("==================================================\n");
    
    // Create the results directory if it doesn't exist
    let results_dir = Path::new("results");
    if !results_dir.exists() {
        if let Err(e) = fs::create_dir_all(results_dir) {
            error!("Failed to create results directory: {}", e);
            return web::Json(json!({
                "success": false,
                "error": format!("Failed to create results directory: {}", e),
                "command": command
            }));
        }
    }
    
    // Execute the command
    let output = match Command::new("sh")
        .arg("-c")
        .arg(&command)
        .output() {
            Ok(output) => output,
            Err(e) => {
                error!("Failed to execute command: {}", e);
                
                // Print the error and command to terminal in a highly visible way
                println!("\n==================================================");
                println!("COMMAND EXECUTION FAILED:");
                println!("{}", command);
                println!("ERROR: {}", e);
                println!("==================================================\n");
                
                return web::Json(json!({
                    "success": false,
                    "error": e.to_string(),
                    "command": command
                }));
            }
        };
    
    // Check if the command was successful
    if output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout).to_string();
        let stderr = String::from_utf8_lossy(&output.stderr).to_string();
        info!("Command executed successfully. Stdout: {}, Stderr: {}", stdout, stderr);
        
        web::Json(json!({
            "success": true,
            "result": if stdout.is_empty() { stderr } else { stdout },
            "command": command
        }))
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr).to_string();
        error!("Analysis failed: {}", stderr);
        
        // Print the error and command to terminal in a highly visible way
        println!("\n==================================================");
        println!("COMMAND FAILED:");
        println!("{}", command);
        println!("ERROR: {}", stderr);
        println!("==================================================\n");
        
        web::Json(json!({
            "success": false,
            "error": stderr,
            "command": command
        }))
    }
}

// Add this function to serve the ancestor graph image
async fn serve_ancestor_graph() -> HttpResponse {
    let path = Path::new("results/common_ontology_graph.pdf");
    
    if path.exists() {
        match fs::read(path) {
            Ok(contents) => {
                HttpResponse::Ok()
                    .content_type("application/pdf")
                    .body(contents)
            }
            Err(e) => {
                error!("Failed to read ancestor graph PDF: {}", e);
                HttpResponse::InternalServerError().body("Failed to read ancestor graph PDF")
            }
        }
    } else {
        error!("Ancestor graph PDF not found");
        // Return a user-friendly message instead of a 404 error
        HttpResponse::Ok()
            .content_type("text/plain")
            .body("No common ancestors found between the provided GO terms")
    }
}

// Add this function to serve the ancestor graph mermaid file
async fn serve_ancestor_mermaid() -> HttpResponse {
    let path = Path::new("results/common_ontology_graph.mmd");
    
    if path.exists() {
        match fs::read_to_string(path) {
            Ok(contents) => {
                // Process the Mermaid content to remove asterisks from GO IDs
                // First, handle the common case of GO IDs with asterisks
                let processed_contents = contents
                    .replace("**GO:", "GO:")
                    .replace("GO:**", "GO:")
                    .replace("**GO", "GO")
                    .replace("GO**", "GO");
                
                // Handle more complex cases where asterisks might be in the middle or at the end
                let mut cleaned = String::with_capacity(processed_contents.len());
                let mut in_go_id = false;
                let mut go_id_buffer = String::new();
                let mut chars = processed_contents.chars().peekable();
                
                while let Some(c) = chars.next() {
                    if in_go_id {
                        // We're inside a GO ID
                        if c == '*' {
                            // Skip asterisks inside GO IDs
                            continue;
                        } else if c.is_ascii_digit() || c == ':' {
                            // Still in the GO ID
                            go_id_buffer.push(c);
                        } else {
                            // We're exiting the GO ID
                            in_go_id = false;
                            
                            // Check if there are trailing asterisks in the buffer
                            let clean_go_id = go_id_buffer.trim_end_matches('*');
                            cleaned.push_str(clean_go_id);
                            cleaned.push(c); // Add the current character
                            go_id_buffer.clear();
                        }
                    } else {
                        // Check if we're entering a GO ID
                        if c == 'G' && chars.peek() == Some(&'O') {
                            // Consume the 'O'
                            chars.next();
                            
                            // Check for the colon
                            if chars.peek() == Some(&':') {
                                in_go_id = true;
                                go_id_buffer.push_str("GO:");
                                chars.next(); // Consume the ':'
                            } else {
                                cleaned.push(c); // Just a 'G'
                                cleaned.push('O'); // Add the 'O' we peeked
                            }
                        } else {
                            cleaned.push(c);
                        }
                    }
                }
                
                // Handle case where the file ends with a GO ID
                if in_go_id && !go_id_buffer.is_empty() {
                    let clean_go_id = go_id_buffer.trim_end_matches('*');
                    cleaned.push_str(clean_go_id);
                }
                
                // Additional cleanup for any remaining patterns
                let final_cleaned = cleaned
                    .replace("GO:*", "GO:")
                    .replace("*GO:", "GO:")
                    .replace("GO:*", "GO:");
                
                HttpResponse::Ok()
                    .content_type("text/plain")
                    .body(final_cleaned)
            }
            Err(e) => {
                error!("Failed to read ancestor graph Mermaid file: {}", e);
                HttpResponse::InternalServerError().body("Failed to read ancestor graph Mermaid file")
            }
        }
    } else {
        error!("Ancestor graph Mermaid file not found");
        // Return a user-friendly message instead of a 404 error
        HttpResponse::Ok()
            .content_type("text/plain")
            .body("No common ancestors found between the provided GO terms")
    }
}

#[derive(Debug, Deserialize)]
struct DeleteFileInfoRequest {
    filename: String,
}

#[post("/delete-file-info")]
async fn delete_file_info(data: web::Json<DeleteFileInfoRequest>) -> impl Responder {
    info!("Deleting file info for: {}", data.filename);
    
    // Try to delete the file from both possible locations
    let data_path = format!("data/{}", data.filename);
    let background_path = format!("data/background_population/{}", data.filename);
    
    let mut success = false;
    let mut deleted_from_background = false;
    
    // Try to delete from data directory
    if Path::new(&data_path).exists() {
        if let Err(e) = fs::remove_file(&data_path) {
            error!("Failed to delete file from data directory: {}", e);
        } else {
            info!("Successfully deleted file from data directory: {}", data_path);
            success = true;
        }
    }
    
    // Try to delete from background_population directory
    if Path::new(&background_path).exists() {
        if let Err(e) = fs::remove_file(&background_path) {
            error!("Failed to delete file from background_population directory: {}", e);
        } else {
            info!("Successfully deleted file from background_population directory: {}", background_path);
            success = true;
            deleted_from_background = true;
        }
    }
    
    // If we deleted a file from the background_population directory,
    // check if the directory is now empty and delete it if it is
    if deleted_from_background {
        let background_dir = Path::new("data/background_population");
        if background_dir.exists() && background_dir.is_dir() {
            // Check if directory is empty
            if let Ok(entries) = fs::read_dir(background_dir) {
                let is_empty = entries.count() == 0;
                if is_empty {
                    // Directory is empty, delete it
                    if let Err(e) = fs::remove_dir(background_dir) {
                        error!("Failed to delete empty background_population directory: {}", e);
                    } else {
                        info!("Successfully deleted empty background_population directory");
                    }
                }
            }
        }
    }
    
    if success {
        HttpResponse::Ok().json(json!({
            "success": true,
            "message": "File deleted successfully"
        }))
    } else {
        HttpResponse::Ok().json(json!({
            "success": false,
            "error": "File not found or could not be deleted"
        }))
    }
}

#[post("/clear-data-folder")]
async fn clear_data_folder() -> HttpResponse {
    info!("Clearing data folder");
    
    // Create a single cleanup command for the data directory
    let cleanup_script = "rm -rf data";

    // Execute cleanup as a single shell command
    let output = std::process::Command::new("sh")
        .arg("-c")
        .arg(cleanup_script)
        .output();

    match output {
        Ok(output) => {
            if output.status.success() {
                info!("✓ Successfully removed data directory");
            } else {
                error!("✗ Cleanup failed: {}", String::from_utf8_lossy(&output.stderr));
                // Fallback to Rust's fs::remove_dir_all
                if let Err(e) = fs::remove_dir_all("data") {
                    error!("Fallback cleanup also failed: {}", e);
                    return HttpResponse::InternalServerError().body("Failed to clear data directory");
                }
            }
        },
        Err(e) => {
            error!("✗ Failed to execute cleanup command: {}", e);
            // Fallback to Rust's fs::remove_dir_all
            if let Err(e) = fs::remove_dir_all("data") {
                error!("Fallback cleanup also failed: {}", e);
                return HttpResponse::InternalServerError().body("Failed to clear data directory");
            }
        }
    }

    // Verify final state
    let data_exists = Path::new("data").exists();
    info!("Data directory exists after cleanup: {}", data_exists);

    if !data_exists {
        info!("✓ Data cleanup successful!");
        HttpResponse::Ok().finish()
    } else {
        error!("✗ Data directory still exists!");
        HttpResponse::InternalServerError().body("Failed to fully clear data directory")
    }
}

#[post("/upload-file")]
async fn upload_file(mut payload: Multipart) -> Result<HttpResponse> {
    info!("Handling file upload");
    
    // Create data directory if it doesn't exist
    if let Err(e) = fs::create_dir_all("data") {
        error!("Failed to create data directory: {}", e);
        return Ok(HttpResponse::InternalServerError().body("Failed to create data directory"));
    }
    
    let mut dropzone_type = String::new();
    let mut filename = String::new();
    let mut file_data = Vec::new();
    
    // Process the multipart form data
    while let Ok(Some(mut field)) = payload.try_next().await {
        let content_disposition = field.content_disposition();
        
        if let Some(name) = content_disposition.get_name() {
            if name == "file" {
                if let Some(fname) = content_disposition.get_filename() {
                    filename = fname.to_string();
                    info!("Processing file: {}", filename);
                    
                    // Read the file content
                    while let Some(chunk) = field.next().await {
                        let data = match chunk {
                            Ok(data) => data,
                            Err(e) => {
                                error!("Error while reading chunk: {}", e);
                                return Ok(HttpResponse::InternalServerError().body(format!("Error while reading chunk: {}", e)));
                            }
                        };
                        file_data.extend_from_slice(&data);
                    }
                }
            } else if name == "dropzone" {
                // Read the dropzone type
                while let Some(chunk) = field.next().await {
                    let data = match chunk {
                        Ok(data) => data,
                        Err(e) => {
                            error!("Error while reading dropzone type: {}", e);
                            return Ok(HttpResponse::InternalServerError().body(format!("Error while reading dropzone type: {}", e)));
                        }
                    };
                    if let Ok(s) = std::str::from_utf8(&data) {
                        dropzone_type = s.to_string();
                        info!("Dropzone type: {}", dropzone_type);
                    }
                }
            } else if name == "filename" {
                // Read the filename if provided separately
                while let Some(chunk) = field.next().await {
                    let data = match chunk {
                        Ok(data) => data,
                        Err(e) => {
                            error!("Error while reading filename: {}", e);
                            return Ok(HttpResponse::InternalServerError().body(format!("Error while reading filename: {}", e)));
                        }
                    };
                    if let Ok(s) = std::str::from_utf8(&data) {
                        if filename.is_empty() {
                            filename = s.to_string();
                            info!("Filename from field: {}", filename);
                        }
                    }
                }
            }
        }
    }
    
    if !filename.is_empty() && !file_data.is_empty() {
        // Determine the file path based on the dropzone type
        let filepath = if dropzone_type == "background-population" {
            // Create background_population directory only if we're uploading to it
            let background_dir = Path::new("data/background_population");
            if !background_dir.exists() {
                if let Err(e) = fs::create_dir_all(background_dir) {
                    error!("Failed to create background_population directory: {}", e);
                    return Ok(HttpResponse::InternalServerError().body(format!("Failed to create background_population directory: {}", e)));
                }
                info!("Created background_population directory");
            }
            
            // Check if this file also exists in the main data directory and remove it if it does
            let data_path = format!("data/{}", filename);
            if Path::new(&data_path).exists() {
                if let Err(e) = fs::remove_file(&data_path) {
                    warn!("Failed to remove duplicate file from data directory: {}", e);
                } else {
                    info!("Removed duplicate file from data directory: {}", data_path);
                }
            }
            
            format!("data/background_population/{}", filename)
        } else {
            format!("data/{}", filename)
        };
        
        info!("Saving file to: {}", filepath);
        
        // Create a file to save the uploaded content
        let mut file = match File::create(&filepath) {
            Ok(file) => file,
            Err(e) => {
                error!("Failed to create file: {}", e);
                return Ok(HttpResponse::InternalServerError().body(format!("Failed to create file: {}", e)));
            }
        };
        
        // Write the file content
        if let Err(e) = file.write_all(&file_data) {
            error!("Error while writing to file: {}", e);
            return Ok(HttpResponse::InternalServerError().body(format!("Error while writing to file: {}", e)));
        }
        
        info!("File saved successfully: {} ({} bytes)", filepath, file_data.len());
        
        // For OBO files, print additional information
        if filename.to_lowercase().ends_with(".obo") {
            if let Ok(metadata) = fs::metadata(&filepath) {
                info!("OBO file size: {} bytes", metadata.len());
                
                // Check if the file is empty or very small
                if metadata.len() < 1000 {
                    warn!("OBO file is suspiciously small (< 1KB). It may be incomplete or corrupted.");
                }
            }
            
            // Check if the file exists
            if Path::new(&filepath).exists() {
                info!("OBO file exists at path: {}", filepath);
            } else {
                error!("OBO file does not exist at path: {}", filepath);
            }
        }
    } else {
        if filename.is_empty() {
            error!("No filename provided in the upload");
        }
        if file_data.is_empty() {
            error!("No file data received in the upload");
        }
        return Ok(HttpResponse::BadRequest().body("Missing filename or file data"));
    }
    
    Ok(HttpResponse::Ok().body("File uploaded successfully"))
}

#[post("/run-semantic-similarity")]
async fn run_semantic_similarity(data: web::Json<SemanticSimilarityRequest>) -> impl Responder {
    info!("Running semantic similarity analysis for GO terms: {:?} with method: {}", data.go_terms, data.method);
    
    if data.go_terms.len() < 2 {
        return web::Json(json!({
            "success": false,
            "error": "At least two GO terms are required for semantic similarity analysis",
            "command": ""
        }));
    }
    
    // Check if any OBO file exists in the data directory
    let mut obo_param = String::new();
    let data_dir = Path::new("data");
    let mut obo_found = false;
    
    if data_dir.exists() && data_dir.is_dir() {
        // Look for .obo files in the data directory
        if let Ok(entries) = fs::read_dir(data_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_file() && 
                   path.extension().map_or(false, |ext| ext == "obo") {
                    if let Some(filename) = path.file_name() {
                        obo_param = format!("-o data/{}", filename.to_string_lossy());
                        obo_found = true;
                        info!("Using OBO file: {}", path.display());
                        break;
                    }
                }
            }
        }
    }
    
    // Check if background population files exist
    let mut background_param = String::new();
    let mut background_found = false;
    
    // First check if background_population directory exists and has background files
    let background_population_dir = Path::new("data/background_population");
    if background_population_dir.exists() && background_population_dir.is_dir() {
        // Check if the directory contains any background files
        if let Ok(entries) = fs::read_dir(background_population_dir) {
            let has_background_files = entries
                .filter_map(Result::ok)
                .any(|entry| {
                    if let Ok(metadata) = entry.metadata() {
                        if metadata.is_file() {
                            if let Some(file_name) = entry.path().file_name() {
                                let file_name_str = file_name.to_string_lossy();
                                return file_name_str.ends_with("_background.txt");
                            }
                        }
                    }
                    false
                });
            
            if has_background_files {
                background_param = format!("-b data/background_population");
                background_found = true;
                info!("Using background population files from background_population directory");
            }
        }
    }
    
    // If no background files found in background_population directory, check the data directory
    if !background_found && data_dir.exists() && data_dir.is_dir() {
        // Look for background files in the data directory
        if let Ok(entries) = fs::read_dir(data_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_file() {
                    if let Some(filename) = path.file_name() {
                        let filename_str = filename.to_string_lossy();
                        if filename_str.ends_with("_background.txt") {
                            background_param = format!("-b data");
                            background_found = true;
                            info!("Using background population files from data directory");
                            break;
                        }
                    }
                }
            }
        }
    }
    
    // Check if propagation was enabled
    let mut propagation_param = String::new();
    
    // Read the propagation setting from a file or environment variable
    let propagation_marker = Path::new("data/.propagation_enabled");
    if propagation_marker.exists() {
        propagation_param = "-p".to_string();
        info!("Propagation is enabled");
    }
    
    // Build the command with the specified parameters
    let go_terms_arg = data.go_terms.join(",");
    
    // Add method parameter if it's not the default (resnik)
    let method_param = if data.method != "resnik" {
        format!("-m {}", data.method)
    } else {
        String::new()
    };
    
    // Construct the full command
    let mut command_parts = Vec::new();
    command_parts.push("semantic-similarity".to_string());
    
    if obo_found {
        command_parts.push(obo_param);
    }
    
    command_parts.push(format!("-t {}", go_terms_arg));
    
    if background_found {
        command_parts.push(background_param);
    }
    
    command_parts.push("-d results".to_string());
    
    if !method_param.is_empty() {
        command_parts.push(method_param);
    }
    
    if !propagation_param.is_empty() {
        command_parts.push(propagation_param);
    } else {
        // Default propagation parameter if none specified
        command_parts.push("-p".to_string());
    }
    
    // Join all parts with spaces
    let command = command_parts.join(" ");
    
    // Print the command to the console for debugging
    println!("\n==================================================");
    println!("SEMANTIC SIMILARITY COMMAND:");
    println!("{}", command);
    println!("==================================================\n");
    
    // Create the results directory if it doesn't exist
    let results_dir = Path::new("results");
    if !results_dir.exists() {
        if let Err(e) = fs::create_dir_all(results_dir) {
            error!("Failed to create results directory: {}", e);
            return web::Json(json!({
                "success": false,
                "error": format!("Failed to create results directory: {}", e),
                "command": command
            }));
        }
    }
    
    // Execute the command
    let output = match Command::new("sh")
        .arg("-c")
        .arg(&command)
        .output() {
            Ok(output) => output,
            Err(e) => {
                error!("Failed to execute command: {}", e);
                return web::Json(json!({
                    "success": false,
                    "error": format!("Failed to execute command: {}", e),
                    "command": command
                }));
            }
        };
    
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        error!("Command execution failed: {}", stderr);
        return web::Json(json!({
            "success": false,
            "error": format!("Command execution failed: {}", stderr),
            "command": command
        }));
    }
    
    // Log the command output for debugging
    let stdout = String::from_utf8_lossy(&output.stdout);
    info!("Command output: {}", stdout);
    
    // Create a map to store similarity scores
    let mut similarity_scores = HashMap::new();
    
    // Read the similarity scores from the TSV file in the results folder
    // The file name might vary based on the method used, so we'll look for any TSV files
    let results_dir = Path::new("results");
    let mut tsv_file_path = None;
    
    if results_dir.exists() && results_dir.is_dir() {
        if let Ok(entries) = fs::read_dir(results_dir) {
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_file() && 
                   path.extension().map_or(false, |ext| ext == "tsv") {
                    // Check if the file name contains "similarity" or the method name
                    if let Some(filename) = path.file_name() {
                        let filename_str = filename.to_string_lossy().to_lowercase();
                        if filename_str.contains("similarity") || 
                           filename_str.contains(&data.method.to_lowercase()) {
                            let path_display = format!("{}", path.display());
                            tsv_file_path = Some(path);
                            info!("Found similarity results file: {}", path_display);
                            break;
                        }
                    }
                }
            }
        }
    }
    
    // If we found a TSV file, read the similarity scores from it
    if let Some(path) = tsv_file_path {
        match fs::read_to_string(path.clone()) {
            Ok(content) => {
                // Parse the TSV file
                // The format might vary, but typically it would have headers and then rows of data
                let lines: Vec<&str> = content.lines().collect();
                
                if lines.len() > 1 {  // At least one header line and one data line
                    // Get the headers (GO terms)
                    let headers: Vec<&str> = lines[0].split('\t').collect();
                    
                    // Process each data row
                    for i in 1..lines.len() {
                        let values: Vec<&str> = lines[i].split('\t').collect();
                        if values.len() > 1 {
                            let term1 = values[0].trim();
                            
                            // Each column (after the first) represents a similarity score with another term
                            for j in 1..values.len() {
                                if j < headers.len() {
                                    let term2 = headers[j].trim();
                                    if let Ok(score) = values[j].trim().parse::<f64>() {
                                        let key = format!("{}-{}", term1, term2);
                                        similarity_scores.insert(key, score);
                                    }
                                }
                            }
                        }
                    }
                }
                
                info!("Parsed {} similarity scores from file", similarity_scores.len());
            },
            Err(e) => {
                error!("Failed to read similarity results file: {}", e);
                // We'll fall back to dummy data below
            }
        }
    } else {
        warn!("No similarity results file found in the results directory");
    }
    
    // If no scores were found, create dummy data for testing
    if similarity_scores.is_empty() {
        warn!("No similarity scores found, using dummy data");
        for i in 0..data.go_terms.len() {
            for j in i..data.go_terms.len() {
                let term1 = &data.go_terms[i];
                let term2 = &data.go_terms[j];
                
                let score = if i == j { 1.0 } else { 0.5 }; // Dummy score
                let key = format!("{}-{}", term1, term2);
                similarity_scores.insert(key, score);
            }
        }
    }
    
    web::Json(json!({
        "success": true,
        "command": command,
        "similarityScores": similarity_scores
    }))
}

#[post("/clear-background-folder")]
async fn clear_background_folder() -> HttpResponse {
    info!("Deleting background_population folder");
    
    let background_dir = Path::new("data/background_population");
    
    // Delete the entire background_population directory if it exists
    if background_dir.exists() {
        if let Err(e) = fs::remove_dir_all(background_dir) {
            error!("Failed to remove background_population directory: {}", e);
            return HttpResponse::InternalServerError().json(json!({
                "success": false,
                "error": format!("Failed to remove background_population directory: {}", e)
            }));
        }
        info!("Successfully removed background_population directory");
    }
    
    // We don't recreate the directory - it will be created when needed
    
    HttpResponse::Ok().json(json!({
        "success": true,
        "message": "Background population folder deleted successfully"
    }))
}

#[get("/test-copy-plots/{species}")]
async fn test_copy_plots(species: web::Path<String>) -> HttpResponse {
    info!("Testing copy_plots with species: {}", species);
    
    match copy_plots(&species).await {
        Ok(_) => {
            HttpResponse::Ok().body(format!("Successfully copied plots for {}", species))
        },
        Err(e) => {
            HttpResponse::InternalServerError().body(format!("Failed to copy plots: {}", e))
        }
    }
}

async fn auto_copy_single_test_results() -> Result<(), std::io::Error> {
    info!("Checking for single test results to auto-copy");
    
    // Check if single_taxon_results directory exists
    let single_dir = Path::new("results/single_taxon_results");
    if !single_dir.exists() {
        info!("No single_taxon_results directory found");
        return Ok(());
    }
    
    // Get all files in the directory
    let entries = match fs::read_dir(single_dir) {
        Ok(entries) => entries,
        Err(e) => {
            error!("Failed to read single_taxon_results directory: {}", e);
            return Err(e);
        }
    };
    
    // Filter for actual files (not directories)
    let files: Vec<_> = entries
        .filter_map(Result::ok)
        .filter(|entry| {
            if let Ok(file_type) = entry.file_type() {
                file_type.is_file()
            } else {
                false
            }
        })
        .collect();
    
    // If there's exactly one file, copy it and its plots
    if files.len() == 1 {
        if let Some(file_name) = files[0].file_name().to_str() {
            info!("Found single test result: {}", file_name);
            
            // Copy the results file
            let source_path = format!("results/single_taxon_results/{}", file_name);
            let dest_path = "results/results.txt";
            
            // Read the source file
            match fs::read_to_string(&source_path) {
                Ok(content) => {
                    // Get the original header
                    let original_header = content.lines().next().unwrap_or("");
                    
                    // Store the original filename and header in a metadata file
                    let metadata_path = "results/original_file_info.txt";
                    let metadata_content = format!("{}\n{}", file_name, original_header);
                    if let Err(e) = fs::write(metadata_path, metadata_content) {
                        error!("Failed to write metadata file: {}", e);
                        // Continue anyway, this is not critical
                    }
                    
                    // For single_taxon_results files, replace the header
                    let mut lines: Vec<&str> = content.lines().collect();
                    if !lines.is_empty() {
                        // Replace the first line (header) with the new header
                        lines[0] = "GO ID\tName\tNamespace\tlog(Odds Ratio)\tStat. Sig.";
                    }
                    let modified_content = lines.join("\n");
                    
                    // Write the modified content to the destination file
                    if let Err(e) = fs::write(dest_path, modified_content) {
                        error!("Failed to write to {}: {}", dest_path, e);
                        return Err(e);
                    }
                    
                    info!("Successfully copied and modified {} to {}", source_path, dest_path);
                    
                    // Also copy the plots
                    if let Err(e) = copy_plots(file_name).await {
                        error!("Failed to copy plots: {}", e);
                        // Continue anyway, this is not critical
                    }
                    
                    return Ok(());
                },
                Err(e) => {
                    error!("Failed to read file from {}: {}", source_path, e);
                    return Err(e);
                }
            }
        }
    } else {
        info!("Found {} files in single_taxon_results directory, not auto-copying", files.len());
    }
    
    Ok(())
}

async fn copy_combined_plots(taxonomy_name: &str) -> Result<(), std::io::Error> {
    info!("Copying combined plots for taxonomy: {}", taxonomy_name);
    
    // Create plots directory if it doesn't exist
    let plots_dir = "results/plots";
    if let Err(e) = fs::create_dir_all(plots_dir) {
        error!("Failed to create plots directory: {}", e);
        return Err(e);
    }
    
    // Remove _GOEA_results.txt from taxonomy name if present
    let clean_taxonomy_name = taxonomy_name.replace("_GOEA_results.txt", "");
    
    // Use underscores for the new naming convention
    let taxonomy_name_with_underscores = clean_taxonomy_name.clone();
    
    // Define the namespaces and their standardized plot names
    let namespaces = [
        ("Biological_Process", "bp"),
        ("Molecular_Function", "mf"),
        ("Cellular_Component", "cc")
    ];
    
    // Define which plots to copy and their new names (new convention)
    let plot_types = [
        ("network_plot.html", "_network.html"),
        ("bubble_plot.html", "_bubble.html"),
        ("bar_plot.html", "_bar.html")
    ];
    
    // Copy plots for each namespace
    for (namespace, prefix) in namespaces.iter() {
        let source_dir = format!("results/combined_taxonomy_results/plots/{}", namespace);
        
        // For each plot type
        for (source_suffix, dest_suffix) in plot_types.iter() {
            let source_path = format!("{}/{}_{}", source_dir, taxonomy_name_with_underscores, source_suffix);
            let dest_path = format!("{}/{}{}", plots_dir, prefix, dest_suffix);
            
            // Copy the file if it exists
            if Path::new(&source_path).exists() {
                info!("Copying {} to {}", source_path, dest_path);
                if let Err(e) = fs::copy(&source_path, &dest_path) {
                    error!("Failed to copy plot {}: {}", source_path, e);
                }
            } else {
                warn!("Plot not found: {}", source_path);
            }
        }
    }
    
    Ok(())
}

#[actix_web::main]
async fn main() -> std::io::Result<()> {
    env_logger::init();
    info!("Starting Gene Ontology Analysis Tool");

    // Create data directory if it doesn't exist
    if let Err(e) = fs::create_dir_all("data") {
        error!("Failed to create data directory: {}", e);
    } else {
        info!("Data directory created or already exists");
    }

    // Set up cleanup handler before starting the server
    let current_dir = std::env::current_dir().expect("Failed to get current directory");
    println!("Current working directory: {}", current_dir.display());
    
    let results_path = current_dir.join("results");
    let data_path = current_dir.join("data");
    
    // Create a test file to verify if handler is called
    let test_file = current_dir.join("cleanup_test.txt");
    fs::write(&test_file, "Cleanup test file").expect("Failed to create test file");
    
    println!("Registered cleanup paths:");
    println!("Results path: {}", results_path.display());
    println!("Data path: {}", data_path.display());
    println!("Test file: {}", test_file.display());

    ctrlc::set_handler(move || {
        // Immediately disable further Ctrl+C handling to prevent interruption
        let _ = ctrlc::set_handler(move || {
            println!("Cleanup in progress, please wait...");
        });

        println!("\n==================================================");
        println!("CTRL+C received! Starting cleanup process...");
        println!("==================================================");

        // Write to test file first
        let _ = fs::write(&test_file, "Handler was called - cleanup in progress");
        
        // Create a single cleanup command for both directories
        let cleanup_script = format!(
            "rm -rf '{}' '{}' '{}'",
            results_path.display(),
            data_path.display(),
            test_file.display()
        );

        // Execute cleanup as a single shell command
        let output = std::process::Command::new("sh")
            .arg("-c")
            .arg(&cleanup_script)
            .output();

        match output {
            Ok(output) => {
                if output.status.success() {
                    println!("✓ Successfully removed all directories and files");
                } else {
                    println!("✗ Cleanup failed: {}", 
                        String::from_utf8_lossy(&output.stderr));
                    
                    // Fallback to individual removal
                    println!("Attempting individual removal...");
                    
                    // Remove results directory
                    if let Err(e) = fs::remove_dir_all(&results_path) {
                        println!("Failed to remove results directory: {}", e);
                    }
                    
                    // Remove data directory
                    if let Err(e) = fs::remove_dir_all(&data_path) {
                        println!("Failed to remove data directory: {}", e);
                    }
                    
                    // Remove test file
                    let _ = fs::remove_file(&test_file);
                }
            },
            Err(e) => {
                println!("✗ Failed to execute cleanup command: {}", e);
                // Fallback to individual removal
                println!("Attempting individual removal...");
                
                // Remove results directory
                if let Err(e) = fs::remove_dir_all(&results_path) {
                    println!("Failed to remove results directory: {}", e);
                }
                
                // Remove data directory
                if let Err(e) = fs::remove_dir_all(&data_path) {
                    println!("Failed to remove data directory: {}", e);
                }
                
                // Remove test file
                let _ = fs::remove_file(&test_file);
            }
        }

        // Verify final state
        let results_exists = results_path.exists();
        let data_exists = data_path.exists();
        let test_exists = test_file.exists();

        println!("\nFinal state verification:");
        println!("Results directory exists: {}", results_exists);
        println!("Data directory exists: {}", data_exists);
        println!("Test file exists: {}", test_exists);

        if !results_exists && !data_exists && !test_exists {
            println!("\n✓ All cleanup successful!");
        } else {
            println!("\n✗ Some items remain!");
        }

        println!("\n==================================================");
        println!("Cleanup process completed. Exiting...");
        println!("==================================================\n");

        // Force exit to ensure we terminate
        std::process::exit(0);
    }).expect("Error setting Ctrl-C handler");

    let app_state = web::Data::new(Arc::new(AppState {
    }));

    let server = HttpServer::new(move || {
        App::new()
            .app_data(app_state.clone())
            .route("/", web::get().to(index))
            .route("/loading", web::get().to(loading))
            .route("/results", web::get().to(results))
            .route("/get_results", web::get().to(get_results))
            .route("/download/{format}", web::get().to(download_results))
            .route("/plots/{name}", web::get().to(serve_plot))
            .route("/get_species", web::get().to(get_species))
            .route("/get_lineage_data", web::get().to(get_lineage_data))
            .route("/get-saved-categories", web::get().to(get_saved_categories))
            .route("/execute", web::post().to(execute_analysis))
            .route("/pdf/ontology_graph.pdf", web::get().to(serve_pdf))
            .route("/results/{namespace}", web::get().to(get_namespace_results))
            .route("/get-combined-files", web::get().to(get_combined_files))
            .route("/get-single-files", web::get().to(get_single_files))
            .route("/ancestor-graph", web::get().to(serve_ancestor_graph))
            .route("/ancestor-mermaid", web::get().to(serve_ancestor_mermaid))
            .service(copy_results_file)
            .service(update_categories)
            .service(actix_files::Files::new("/data", "data").show_files_listing())
            .service(actix_files::Files::new("/demo", "demo").show_files_listing())  // Add this line to serve demo files
            .route("/clear-results", web::post().to(clear_results))
            .service(check_testing_type)
            .service(run_ancestor_analysis)
            .service(run_semantic_similarity)
            .service(delete_file_info)
            .service(clear_data_folder)
            .service(upload_file)
            .service(clear_background_folder)
            .service(test_copy_plots)
    })
    .bind("127.0.0.1:8080")?
    .run();

    if webbrowser::open("http://127.0.0.1:8080").is_err() {
        println!("Failed to open browser automatically");
    }

    server.await
}
