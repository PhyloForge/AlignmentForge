use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use crate::models::{Alignment, AlignmentFormat};

pub fn parse_nexus<P: AsRef<Path>>(path: P) -> Result<Alignment, String> {
    let path_ref = path.as_ref();
    let content = fs::read_to_string(path_ref)
        .map_err(|e| format!("Failed to read NEXUS file: {}", e))?;

    let mut in_matrix = false;
    let mut taxa: Vec<String> = Vec::new();
    let mut sequences: Vec<String> = Vec::new();
    let mut taxon_index_map: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    for line in content.lines() {
        let trimmed = line.trim();

        if trimmed.is_empty() || trimmed.starts_with('[') {
            continue;
        }

        let upper = trimmed.to_ascii_uppercase();

        if upper.starts_with("MATRIX") {
            in_matrix = true;
            continue;
        }

        if in_matrix {
            if trimmed == ";" || upper.starts_with("END;") || upper == "END" {
                break;
            }

            let tokens: Vec<&str> = trimmed.split_whitespace().collect();
            if tokens.len() >= 2 {
                let name = tokens[0].trim_matches('\'').trim_matches('"').to_string();
                let seq = tokens[1..].concat();

                if let Some(&idx) = taxon_index_map.get(&name) {
                    sequences[idx].push_str(&seq);
                } else {
                    let idx = taxa.len();
                    taxon_index_map.insert(name.clone(), idx);
                    taxa.push(name);
                    sequences.push(seq);
                }
            }
        }
    }

    if taxa.is_empty() {
        return Err("No matrix data found in NEXUS file".to_string());
    }

    let file_name = path_ref
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());
    let id = path_ref
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());

    Ok(Alignment::new(
        id,
        file_name,
        path_ref.to_string_lossy().to_string(),
        AlignmentFormat::Nexus,
        taxa,
        sequences,
    ))
}

pub fn write_nexus<P: AsRef<Path>>(
    path: P,
    taxa: &[String],
    sequences: &[String],
    interleaved: bool,
) -> Result<(), String> {
    let mut file = File::create(path).map_err(|e| format!("Failed to create NEXUS: {}", e))?;

    let num_taxa = taxa.len();
    let length = sequences.first().map_or(0, |s| s.len());

    writeln!(file, "#NEXUS\n").map_err(|e| e.to_string())?;
    writeln!(file, "BEGIN DATA;").map_err(|e| e.to_string())?;
    writeln!(
        file,
        "  DIMENSIONS NTAX={} NCHAR={};",
        num_taxa, length
    )
    .map_err(|e| e.to_string())?;
    writeln!(
        file,
        "  FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=. {};",
        if interleaved { "INTERLEAVE" } else { "" }
    )
    .map_err(|e| e.to_string())?;
    writeln!(file, "  MATRIX").map_err(|e| e.to_string())?;

    let max_name_len = taxa.iter().map(|t| t.len()).max().unwrap_or(10);
    let pad = (max_name_len + 4).max(12);

    if !interleaved {
        for (name, seq) in taxa.iter().zip(sequences.iter()) {
            writeln!(file, "    {:<pad$} {}", name, seq, pad = pad).map_err(|e| e.to_string())?;
        }
    } else {
        let chunk_size = 60;
        let total_chunks = (length + chunk_size - 1) / chunk_size;

        for chunk_idx in 0..total_chunks {
            let start = chunk_idx * chunk_size;
            let end = (start + chunk_size).min(length);

            for (name, seq) in taxa.iter().zip(sequences.iter()) {
                let seq_chunk = if start < seq.len() {
                    let chunk_end = end.min(seq.len());
                    &seq[start..chunk_end]
                } else {
                    ""
                };
                writeln!(file, "    {:<pad$} {}", name, seq_chunk, pad = pad)
                    .map_err(|e| e.to_string())?;
            }
            if chunk_idx + 1 < total_chunks {
                writeln!(file).map_err(|e| e.to_string())?;
            }
        }
    }

    writeln!(file, "  ;").map_err(|e| e.to_string())?;
    writeln!(file, "END;").map_err(|e| e.to_string())?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nexus_roundtrip() {
        let temp_dir = std::env::temp_dir();
        let test_file = temp_dir.join("test_alignment.nex");

        let taxa = vec!["Taxon_A".to_string(), "Taxon_B".to_string()];
        let seqs = vec!["ATGCATGC".to_string(), "ATGCATGT".to_string()];

        write_nexus(&test_file, &taxa, &seqs, false).unwrap();
        let parsed = parse_nexus(&test_file).unwrap();

        assert_eq!(parsed.taxa, taxa);
        assert_eq!(parsed.sequences, seqs);
        assert_eq!(parsed.length, 8);
        assert_eq!(parsed.num_taxa, 2);

        let _ = std::fs::remove_file(test_file);
    }
}
