pub mod fasta;
pub mod nexus;
pub mod phylip;

use std::path::Path;
use crate::models::{Alignment, AlignmentFormat};

pub fn detect_format<P: AsRef<Path>>(path: P) -> Result<AlignmentFormat, String> {
    let path_ref = path.as_ref();
    let ext = path_ref
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();

    match ext.as_str() {
        "fa" | "fasta" | "fna" | "ffn" | "faa" => Ok(AlignmentFormat::Fasta),
        "phy" | "phylip" => Ok(AlignmentFormat::Phylip),
        "nex" | "nexus" => Ok(AlignmentFormat::Nexus),
        _ => {
            // Inspect first non-empty lines
            let file = std::fs::File::open(path_ref)
                .map_err(|e| format!("Failed to open file for format detection: {}", e))?;
            use std::io::BufRead;
            let reader = std::io::BufReader::new(file);

            for line in reader.lines().filter_map(|l| l.ok()) {
                let trimmed = line.trim();
                if trimmed.is_empty() {
                    continue;
                }
                if trimmed.starts_with('>') || trimmed.starts_with(';') {
                    return Ok(AlignmentFormat::Fasta);
                }
                if trimmed.to_ascii_uppercase().starts_with("#NEXUS") {
                    return Ok(AlignmentFormat::Nexus);
                }
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.len() == 2 && parts[0].parse::<usize>().is_ok() && parts[1].parse::<usize>().is_ok() {
                    return Ok(AlignmentFormat::Phylip);
                }
            }

            // Fallback default
            Ok(AlignmentFormat::Fasta)
        }
    }
}

pub fn parse_alignment<P: AsRef<Path>>(path: P) -> Result<Alignment, String> {
    let format = detect_format(&path)?;
    match format {
        AlignmentFormat::Fasta => fasta::parse_fasta(path),
        AlignmentFormat::Phylip => phylip::parse_phylip(path),
        AlignmentFormat::Nexus => nexus::parse_nexus(path),
    }
}

pub fn write_alignment<P: AsRef<Path>>(
    path: P,
    taxa: &[String],
    sequences: &[String],
    format: AlignmentFormat,
) -> Result<(), String> {
    match format {
        AlignmentFormat::Fasta => fasta::write_fasta(path, taxa, sequences, None),
        AlignmentFormat::Phylip => phylip::write_phylip(path, taxa, sequences, false),
        AlignmentFormat::Nexus => nexus::write_nexus(path, taxa, sequences, false),
    }
}
