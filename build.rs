use std::env;
use std::fs;
use std::path::Path;
use flate2::read::GzDecoder;
use tar::Archive;
use std::io::{self, Error, ErrorKind};

fn main() {
    let cargo_home = env::var("CARGO_HOME").unwrap_or_else(|_| {
        println!("CARGO_HOME environment variable not set. Using default path.");
        String::from(".cargo")
    });

    let source_file = Path::new("taxago_assets.tar.gz");

    println!("cargo:rerun-if-changed=taxago_assets.tar.gz");

    if !source_file.exists() {
        eprintln!("Error: Assets file doesn't exist: {}", source_file.display());
        return;
    }

    let dest_path = Path::new(&cargo_home).join("taxago_assets");
    let extract_marker = dest_path.join(".extracted");

    fs::create_dir_all(&dest_path).expect("Failed to create destination directory");

    let should_extract = if extract_marker.exists() {
        let marker_modified = fs::metadata(&extract_marker)
            .and_then(|m| m.modified())
            .expect("Failed to get marker modification time");

        let source_modified = fs::metadata(&source_file)
            .and_then(|m| m.modified())
            .expect("Failed to get source modification time");
        source_modified > marker_modified
    } else {
        true
    };

    if should_extract {
        println!("Extracting assets...");
        match extract_tar_gz(&source_file, &dest_path) {
            Ok(_) => {
                fs::write(&extract_marker, "").expect("Failed to create extraction marker");
                println!("Assets extracted to: {}", dest_path.display());
            }
            Err(e) => eprintln!("Failed to extract assets: {}", e),
        }
    } else {
        println!("Assets are up to date");
    }
}

fn extract_tar_gz(source: &Path, target: &Path) -> io::Result<()> {
    let tar_gz = fs::File::open(source)?;
    let tar = GzDecoder::new(tar_gz);
    let mut archive = Archive::new(tar);

    for entry in archive.entries()? {
        let mut entry = entry?;
        let path = entry.path().map_err(|e| Error::new(ErrorKind::InvalidData, e))?;

        let stripped_path = match path.strip_prefix("taxago_assets/") {
            Ok(stripped) => stripped,
            Err(_) => {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    format!("Failed to strip prefix from path: {}", path.display()),
                ));
            }
        };

        let full_path = target.join(stripped_path);

        if let Some(parent) = full_path.parent() {
            fs::create_dir_all(parent)?;
        }

        entry.unpack(&full_path)?;
    }

    Ok(())
}