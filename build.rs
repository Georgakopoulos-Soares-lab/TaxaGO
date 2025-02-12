use std::env;
use std::fs;
use std::path::Path;
use flate2::read::GzDecoder;
use tar::Archive;
use std::io::{self, Error, ErrorKind};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    
    let cargo_home = env::var("CARGO_HOME").unwrap_or_else(|_| {
        println!("cargo:warning=CARGO_HOME environment variable not set. Using default path.");
        String::from(".cargo")
    });

    let source_file = Path::new("taxago_assets.tar.gz");
    println!("cargo:rerun-if-changed=taxago_assets.tar.gz");

    if !source_file.exists() {
        return Err(format!("Assets file doesn't exist: {}", source_file.display()).into());
    }

    let dest_path = Path::new(&cargo_home).join("taxago_assets");
    let extract_marker = dest_path.join(".extracted");

    fs::create_dir_all(&dest_path)
        .map_err(|e| format!("Failed to create destination directory: {}", e))?;

    let should_extract = if extract_marker.exists() {
        let marker_modified = fs::metadata(&extract_marker)
            .and_then(|m| m.modified())
            .map_err(|e| format!("Failed to get marker modification time: {}", e))?;

        let source_modified = fs::metadata(&source_file)
            .and_then(|m| m.modified())
            .map_err(|e| format!("Failed to get source modification time: {}", e))?;

        source_modified > marker_modified
    } else {
        true
    };

    if should_extract {
        println!("cargo:warning=Extracting assets...");
        extract_tar_gz(&source_file, &dest_path)
            .map_err(|e| format!("Failed to extract assets: {}", e))?;
        
        fs::write(&extract_marker, "")
            .map_err(|e| format!("Failed to create extraction marker: {}", e))?;
        
        println!("cargo:warning=Assets extracted to: {}", dest_path.display());
    } 

    Ok(())
}

fn extract_tar_gz(source: &Path, target: &Path) -> io::Result<()> {
    let tar_gz = fs::File::open(source)?;
    let tar = GzDecoder::new(tar_gz);
    let mut archive = Archive::new(tar);

    for entry in archive.entries()? {
        let mut entry = entry?;
        let path = entry.path()
            .map_err(|e| Error::new(ErrorKind::InvalidData, e))?;

        let stripped_path = path.strip_prefix("taxago_assets/")
            .map_err(|_| Error::new(
                ErrorKind::InvalidData,
                format!("Failed to strip prefix from path: {}", path.display())
            ))?;

        let full_path = target.join(stripped_path);

        if let Some(parent) = full_path.parent() {
            fs::create_dir_all(parent)?;
        }

        entry.unpack(&full_path)?;
    }

    Ok(())
}