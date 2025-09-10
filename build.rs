use std::env::var;
use std::fs;
use std::io::{BufReader, copy};
use std::path::Path;
use flate2::bufread::GzDecoder;
use tar::Archive;
use std::error::Error;
use std::time::Duration;
use std::collections::HashSet;
use duct::cmd;

fn main() -> Result<(), Box<dyn Error>> {
    let cargo_home = var("CARGO_HOME")?;
    let source_file = Path::new("taxago_assets.tar.gz");
    let dest_path = Path::new(&cargo_home).join("taxago_assets");

    if !source_file.exists() {
        let url = "https://zenodo.org/records/14860780/files/taxago_assets.tar.gz?download=1";
        download_from_zenodo(url, source_file)?;
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
    
    run_mermaid_install()?;
    Ok(())
}

fn download_from_zenodo(url: &str, out_path: &Path) -> Result<(), Box<dyn Error>> {
    let client = reqwest::blocking::Client::builder()
        .user_agent("build.rs downloader")
        .timeout(Duration::from_secs(300))
        .redirect(reqwest::redirect::Policy::limited(10))
        .build()?;

    let mut resp = client.get(url).send()?.error_for_status()?;
    let mut file = fs::File::create(out_path)?;
    copy(&mut resp, &mut file)?;
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

fn run_mermaid_install() -> Result<(), Box<dyn Error>> {
    match std::env::consts::OS {
        "windows" => install_mermaid_windows(),
        "macos" => install_mermaid_unix(),
        "linux" => install_mermaid_unix(),
        os => {
            println!("cargo:warning=Unsupported OS: {}. Skipping mermaid installation.", os);
            Ok(())
        }
    }
}

fn install_mermaid_windows() -> Result<(), Box<dyn Error>> {
    if cmd!("where", "npm").run().is_err() {
        println!("cargo:warning=npm not found. Attempting to install Node.js via Chocolatey...");
        
        if cmd!("where", "choco").run().is_err() {
            println!("cargo:warning=Installing Chocolatey...");
            let install_result = cmd!("powershell", "-ExecutionPolicy", "Bypass", "-c", 
                "Set-ExecutionPolicy Bypass -Scope Process -Force; [System.Net.ServicePointManager]::SecurityProtocol = [System.Net.ServicePointManager]::SecurityProtocol -bor 3072; iex ((New-Object System.Net.WebClient).DownloadString('https://community.chocolatey.org/install.ps1'))")
                .run();
            
            if install_result.is_err() {
                println!("cargo:warning=Failed to install Chocolatey. Please install Node.js manually.");
                return Ok(());
            }
        }
        
        let nodejs_result = cmd!("choco", "install", "nodejs", "--version=20.18.0", "-y").run();
        if nodejs_result.is_err() {
            println!("cargo:warning=Failed to install Node.js via Chocolatey. Please install manually.");
            return Ok(());
        }
        
        cmd!("powershell", "-c", "refreshenv").run().ok();
    }

    if cmd!("where", "mmdc").run().is_err() {
        println!("cargo:warning=Installing Mermaid CLI...");
        let mermaid_result = cmd!("npm", "install", "-g", "@mermaid-js/mermaid-cli").run();
        if mermaid_result.is_err() {
            println!("cargo:warning=Failed to install Mermaid CLI. Please install manually with: npm install -g @mermaid-js/mermaid-cli");
        }
    }

    Ok(())
}

fn install_mermaid_unix() -> Result<(), Box<dyn Error>> {
    if cmd!("which", "npm").run().is_err() {
        println!("cargo:warning=npm not found. Attempting to install Node.js via nvm...");
        
        let nvm_install_cmd = r#"
            # Install nvm if not present
            if ! command -v nvm &> /dev/null; then
                curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
                export NVM_DIR="$HOME/.nvm"
                [ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"
            fi
            
            # Install and use Node.js
            nvm install 20
            nvm use 20
        "#;
        
        let install_result = cmd!("bash", "-c", nvm_install_cmd).run();
        if install_result.is_err() {
            println!("cargo:warning=Failed to install Node.js via nvm. Please install Node.js manually.");
            return Ok(());
        }
    }

    if cmd!("which", "mmdc").run().is_err() {
        println!("cargo:warning=Installing Mermaid CLI...");
        
        let mermaid_install_cmd = r#"
            # Source nvm if available
            export NVM_DIR="$HOME/.nvm"
            [ -s "$NVM_DIR/nvm.sh" ] && \. "$NVM_DIR/nvm.sh"
            
            # Install mermaid CLI
            npm install -g @mermaid-js/mermaid-cli
        "#;
        
        let mermaid_result = cmd!("bash", "-c", mermaid_install_cmd).run();
        if mermaid_result.is_err() {
            println!("cargo:warning=Failed to install Mermaid CLI. Please install manually with: npm install -g @mermaid-js/mermaid-cli");
        }
    }

    Ok(())
}