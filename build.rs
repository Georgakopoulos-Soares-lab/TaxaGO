use std::env::var;
use std::fs;
use std::io::BufReader;
use std::path::Path;
use flate2::bufread::GzDecoder;
use tar::Archive;
use std::error::Error;
use std::collections::HashSet;

fn main() -> Result<(), Box<dyn Error>> {

    let cargo_home = var("CARGO_HOME")?;
    let source_file = Path::new("taxago_assets.tar.gz");
    let dest_path = Path::new(&cargo_home).join("taxago_assets");

    if !source_file.exists() {
        println!("cargo::error=taxago_assets.tar.gz is not found! Make sure you have downloaded it from the Zenodo")
    }

    let extract_marker = dest_path.join(".extracted");

    fs::create_dir_all(&dest_path)?;

    let should_extract = if extract_marker.exists() {
        let marker_modified = fs::metadata(&extract_marker)
            .and_then(|m| m.modified())?;

        let source_modified = fs::metadata(&source_file)
            .and_then(|m| m.modified())?;

        source_modified > marker_modified
    } else {
        true
    };

    if should_extract {
        extract_taxago_assets(&source_file, &dest_path)?;
        fs::write(&extract_marker, "")?;

    } 

    Ok(())
}

fn extract_taxago_assets(
    source: &Path, 
    destination: &Path
) -> Result<(), Box<dyn Error>> {
    
    let taxago_assets_tar_gz = fs::File::open(source)?;
    let taxago_assets_buffer = BufReader::with_capacity(256 * 1024 * 1024, taxago_assets_tar_gz);
    let taxago_assets_tar = GzDecoder::new(taxago_assets_buffer);
    let mut taxago_assets = Archive::new(taxago_assets_tar);

    let mut created_dirs = HashSet::new();
    for entry in taxago_assets.entries()? {
        let mut entry = entry?;
        let path = entry.path()?;
        
        if path.file_name()
            .and_then(|name| name.to_str())
            .map(|name| name.starts_with("._"))
            .unwrap_or(false) {
            continue;
        }

        let full_path = destination.join(&path);

        if let Some(parent) = full_path.parent() {
            if created_dirs.insert(parent.to_path_buf()) {
                fs::create_dir_all(parent)?;
            }
        }

        entry.unpack(&full_path)?;
    }

    Ok(())
}