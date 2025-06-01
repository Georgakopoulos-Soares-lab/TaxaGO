use serde::{Serialize, Deserialize};
use std::fs;
use std::path::Path;
use anyhow::{Result, Context};
use csv::WriterBuilder;
use std::io::Write;

#[derive(Debug, Serialize, Deserialize)]
pub struct GOResult {
    go_term: String,
    name: String,
    namespace: String,
    odds_ratio: f64,
    statistical_significance: f64,
}

pub struct DownloadManager {
    results_path: String,
}

impl DownloadManager {
    pub fn new(results_path: String) -> Self {
        Self { results_path }
    }

    pub fn generate_file_content(&self, format: &str) -> Result<(String, Vec<u8>)> {
        let (original_filename, original_header) = self.get_original_file_info()?;
        
        println!("Original filename: {}", original_filename);
        println!("Original header: {}", original_header);
        
        let output_filename = match format {
            "csv" => format!("{}.csv", original_filename),
            "tsv" => format!("{}.tsv", original_filename),
            "json" => format!("{}.json", original_filename),
            _ => format!("{}.{}", original_filename, format),
        };
        
        println!("Output filename: {}", output_filename);
        
        let results = self.read_results()?;
        
        match format {
            "csv" => {
                let mut wtr = WriterBuilder::new()
                    .has_headers(false) 
                    .from_writer(vec![]);
                
                let header_fields: Vec<&str> = original_header.split('\t').collect();
                wtr.write_record(&header_fields)?;
                
                for result in &results {
                        wtr.write_record(&[
                            &result.go_term,
                            &result.name,
                            &result.namespace,
                            &format!("{:.3}", result.odds_ratio),
                            &result.statistical_significance.to_string(),
                        ])?;
                }
                
                let content = wtr.into_inner()?;
                Ok((output_filename, content))
            }
            "tsv" => {
                let mut content = Vec::new();
                
                writeln!(content, "{}", original_header)?;
                
                for result in &results {
                    writeln!(content, "{}\t{}\t{}\t{:.3}\t{}",
                        result.go_term,
                        result.name,
                        result.namespace,
                        result.odds_ratio,
                        result.statistical_significance)?;
                }
            
                Ok((output_filename, content))
            }
            "json" => {
                let content = serde_json::to_vec_pretty(&results)?;
                Ok((output_filename, content))
            }
            _ => Err(anyhow::anyhow!("Unsupported format"))?
        }
    }

    fn get_original_file_info(&self) -> Result<(String, String)> {
        let results_dir = Path::new(&self.results_path).parent().unwrap_or(Path::new(""));
        
        let mut original_filename = String::from("results");
        let mut original_header;
        
        let content = fs::read_to_string(&self.results_path)
            .context("Failed to read results file")?;
        
        let lines: Vec<&str> = content.lines().collect();
        if !lines.is_empty() {
            original_header = lines[0].to_string();
        } else {
            return Err(anyhow::anyhow!("Results file is empty"));
        }
        
        let extract_taxonomy_name = |filename: &str| -> String {
            if let Some(suffix_pos) = filename.find("_GOEA_results.txt") {
                filename[0..suffix_pos].to_string()
            } else if let Some(name_without_ext) = filename.strip_suffix("_GOEA_results.txt") {
                name_without_ext.to_string()
            } else if let Some(name_without_ext) = filename.strip_suffix(".txt") {
                name_without_ext.to_string()
            } else {
                let clean_name = filename.trim_end_matches(|c: char| !c.is_alphanumeric() && c != '_' && c != '-' && c != '.');
                clean_name.to_string()
            }
        };
        
        let metadata_path = results_dir.join("original_file_info.txt");
        if metadata_path.exists() {
            if let Ok(metadata) = fs::read_to_string(&metadata_path) {
                let metadata_lines: Vec<&str> = metadata.lines().collect();
                
                if let Some(name) = metadata_lines.get(0) {
                    original_filename = extract_taxonomy_name(name);
                    println!("Original filename from metadata: {}", original_filename);
                }
                
                if let Some(header) = metadata_lines.get(1) {
                    original_header = header.to_string();
                    println!("Original header from metadata: {}", original_header);
                }
            }
        } else {
            let single_taxon_dir = results_dir.join("single_taxon_results");
            let combined_dir = results_dir.join("combined_taxonomy_results");
            
            if single_taxon_dir.exists() {
                if let Ok(entries) = fs::read_dir(single_taxon_dir) {
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
                    
                    if files.len() == 1 {
                        if let Some(filename) = files[0].file_name().to_str() {
                            original_filename = extract_taxonomy_name(filename);
                            println!("Using taxonomy name for download: {}", original_filename);
                        }
                    }
                }
            } else if combined_dir.exists() {
                if let Ok(entries) = fs::read_dir(combined_dir) {
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
                    
                    if files.len() == 1 {
                        if let Some(filename) = files[0].file_name().to_str() {
                            original_filename = extract_taxonomy_name(filename);
                            println!("Using taxonomy name for download: {}", original_filename);
                        }
                    }
                }
            }
        }
        
        Ok((original_filename, original_header))
    }

    fn read_results(&self) -> Result<Vec<GOResult>> {
        let content = fs::read_to_string(&self.results_path)
            .context("Failed to read results file")?;

        let mut results = Vec::new();
        let lines: Vec<&str> = content.lines().collect();
        
        if lines.is_empty() {
            return Err(anyhow::anyhow!("Results file is empty"));
        }    
            
        for line in lines.iter().skip(1) { 
            let fields: Vec<&str> = line.split('\t').collect();
            let odds_ratio: f64 = fields[3].parse().unwrap_or(0.0);
                results.push(GOResult {
                    go_term: fields[0].to_string(),
                    name: fields[1].to_string(),
                    namespace: fields[2].to_string(),
                    odds_ratio: (odds_ratio * 1000.0).round() / 1000.0,
                    statistical_significance: fields[4].parse().unwrap_or(0.0)});
        }

        results.sort_by(|a, b| {
            let odds_ratio_cmp = b.odds_ratio.partial_cmp(&a.odds_ratio).unwrap_or(std::cmp::Ordering::Equal);
            if odds_ratio_cmp == std::cmp::Ordering::Equal {
                a.statistical_significance.partial_cmp(&b.statistical_significance).unwrap_or(std::cmp::Ordering::Equal)
            } else {
                odds_ratio_cmp
            }
        });

        Ok(results)
    }
}
