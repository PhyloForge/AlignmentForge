pub mod fasta;
pub mod nexus;
pub mod phylip;

use std::path::Path;
use crate::models::{Alignment, AlignmentFormat};

/// Rejects an alignment whose rows do not form a rectangular matrix.
///
/// Every downstream step indexes columns across all samples, so a ragged matrix
/// silently produces wrong statistics rather than an error. `declared_length`
/// is checked only when a format states it in its header.
pub fn validate_alignment_shape(
    taxa: &[String],
    sequences: &[String],
    declared_length: Option<usize>,
    format_name: &str,
) -> Result<(), String> {
    if taxa.len() != sequences.len() {
        return Err(format!(
            "{format_name} file has {} taxon names but {} sequences",
            taxa.len(),
            sequences.len()
        ));
    }

    let Some(first_length) = sequences.first().map(String::len) else {
        return Ok(());
    };

    for (taxon, sequence) in taxa.iter().zip(sequences.iter()) {
        if sequence.len() != first_length {
            return Err(format!(
                "{format_name} alignment is not rectangular: '{taxon}' has {} sites but the first sequence has {first_length}",
                sequence.len()
            ));
        }
    }

    if let Some(declared) = declared_length {
        if declared > 0 && declared != first_length {
            return Err(format!(
                "{format_name} header declares {declared} sites but the sequences have {first_length}"
            ));
        }
    }

    Ok(())
}

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

/// Detects a format from the file name and, failing that, from the content.
/// The browser build has no filesystem, so detection must work from memory.
pub fn detect_format_from_text(file_name: &str, content: &str) -> AlignmentFormat {
    let ext = file_name
        .rsplit('.')
        .next()
        .map(|s| s.to_ascii_lowercase())
        .unwrap_or_default();

    match ext.as_str() {
        "fa" | "fasta" | "fna" | "ffn" | "faa" => return AlignmentFormat::Fasta,
        "phy" | "phylip" => return AlignmentFormat::Phylip,
        "nex" | "nexus" => return AlignmentFormat::Nexus,
        _ => {}
    }

    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('>') || trimmed.starts_with(';') {
            return AlignmentFormat::Fasta;
        }
        if trimmed.to_ascii_uppercase().starts_with("#NEXUS") {
            return AlignmentFormat::Nexus;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() == 2
            && parts[0].parse::<usize>().is_ok()
            && parts[1].parse::<usize>().is_ok()
        {
            return AlignmentFormat::Phylip;
        }
        break;
    }

    AlignmentFormat::Fasta
}

/// Parses an alignment held in memory. Used by the browser build and by tests.
pub fn parse_alignment_text(
    content: &str,
    id: &str,
    file_name: &str,
    file_path: &str,
) -> Result<Alignment, String> {
    match detect_format_from_text(file_name, content) {
        AlignmentFormat::Fasta => {
            fasta::parse_fasta_bytes(content.as_bytes(), id, file_name, file_path)
        }
        AlignmentFormat::Phylip => phylip::parse_phylip_str(content, id, file_name, file_path),
        AlignmentFormat::Nexus => nexus::parse_nexus_str(content, id, file_name, file_path),
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

#[cfg(test)]
mod tests {
    #[test]
    fn example_data_parses_cleanly() {
        let dir = std::path::Path::new("../public/example_data");
        if !dir.exists() {
            return;
        }
        let mut files = 0;
        let mut ok = 0;
        for entry in walkdir::WalkDir::new(dir).into_iter().filter_map(|e| e.ok()) {
            if !entry.file_type().is_file() {
                continue;
            }
            let p = entry.path();
            if p.extension().and_then(|s| s.to_str()) != Some("phy") {
                continue;
            }
            files += 1;
            match crate::parsers::parse_alignment(p) {
                Ok(a) => {
                    ok += 1;
                    assert!(a.num_taxa > 0, "{:?} parsed 0 taxa", p);
                    assert!(a.length > 0, "{:?} parsed 0 length", p);
                    for s in &a.sequences {
                        assert_eq!(s.len(), a.length, "{:?} ragged", p);
                    }
                }
                Err(e) => panic!("{:?} failed to parse: {}", p, e),
            }
        }
        eprintln!("EXAMPLE DATA: {}/{} parsed", ok, files);
        assert_eq!(files, 29);
    }

}
