use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use crate::models::{Alignment, AlignmentFormat};

pub fn parse_phylip<P: AsRef<Path>>(path: P) -> Result<Alignment, String> {
    let path_ref = path.as_ref();
    let content = fs::read_to_string(path_ref)
        .map_err(|e| format!("Failed to read PHYLIP file: {}", e))?;

    let mut lines = content
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty());

    let header = lines.next().ok_or_else(|| "Empty PHYLIP file".to_string())?;
    let parts: Vec<&str> = header.split_whitespace().collect();
    if parts.len() < 2 {
        return Err(format!("Invalid PHYLIP header: '{}'", header));
    }

    let expected_taxa: usize = parts[0]
        .parse()
        .map_err(|_| format!("Invalid taxa count in PHYLIP header: {}", parts[0]))?;
    let _expected_length: usize = parts[1]
        .parse()
        .map_err(|_| format!("Invalid length in PHYLIP header: {}", parts[1]))?;

    let mut taxa = Vec::with_capacity(expected_taxa);
    let mut sequences = Vec::with_capacity(expected_taxa);

    let remaining_lines: Vec<&str> = lines.collect();

    // Check if sequential or interleaved
    if remaining_lines.len() == expected_taxa {
        // Sequential: 1 line per taxon
        for line in remaining_lines {
            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }
            if tokens.len() == 1 {
                // Strict fixed-10 width fallback
                if line.len() > 10 {
                    let (name, seq) = line.split_at(10);
                    taxa.push(name.trim().to_string());
                    sequences.push(seq.split_whitespace().collect::<String>());
                } else {
                    taxa.push(tokens[0].to_string());
                    sequences.push(String::new());
                }
            } else if tokens.len() == 2 {
                taxa.push(tokens[0].to_string());
                sequences.push(tokens[1].to_string());
            } else {
                taxa.push(tokens[0].to_string());
                let seq = tokens[1..].concat();
                sequences.push(seq);
            }
        }
    } else {
        // Multi-line sequential or interleaved
        let mut cur_taxon_idx = 0;

        for line in remaining_lines {
            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }

            if taxa.len() < expected_taxa {
                if tokens.len() >= 2 {
                    taxa.push(tokens[0].to_string());
                    sequences.push(tokens[1..].concat());
                } else if line.len() > 10 {
                    let (name, seq) = line.split_at(10);
                    taxa.push(name.trim().to_string());
                    sequences.push(seq.split_whitespace().collect::<String>());
                }
            } else {
                if cur_taxon_idx >= expected_taxa {
                    cur_taxon_idx = 0;
                }
                let seq_chunk = if tokens.len() >= 2 && tokens[0] == taxa[cur_taxon_idx] {
                    tokens[1..].concat()
                } else {
                    tokens.concat()
                };
                sequences[cur_taxon_idx].push_str(&seq_chunk);
                cur_taxon_idx += 1;
            }
        }
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
        AlignmentFormat::Phylip,
        taxa,
        sequences,
    ))
}

pub fn write_phylip<P: AsRef<Path>>(
    path: P,
    taxa: &[String],
    sequences: &[String],
    interleaved: bool,
) -> Result<(), String> {
    let mut file = File::create(path).map_err(|e| format!("Failed to create PHYLIP: {}", e))?;

    let num_taxa = taxa.len();
    let length = sequences.first().map_or(0, |s| s.len());

    // Relaxed header
    writeln!(file, "{} {}", num_taxa, length).map_err(|e| e.to_string())?;

    let max_name_len = taxa.iter().map(|t| t.len()).max().unwrap_or(10);
    let pad = (max_name_len + 4).max(12);

    if !interleaved {
        // Sequential relaxed PHYLIP
        for (name, seq) in taxa.iter().zip(sequences.iter()) {
            writeln!(file, "{:<pad$} {}", name, seq, pad = pad).map_err(|e| e.to_string())?;
        }
    } else {
        // Interleaved PHYLIP in 60bp chunks
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

                if chunk_idx == 0 {
                    writeln!(file, "{:<pad$} {}", name, seq_chunk, pad = pad)
                        .map_err(|e| e.to_string())?;
                } else {
                    writeln!(file, "{:<pad$} {}", "", seq_chunk, pad = pad)
                        .map_err(|e| e.to_string())?;
                }
            }
            if chunk_idx + 1 < total_chunks {
                writeln!(file).map_err(|e| e.to_string())?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phylip_roundtrip() {
        let temp_dir = std::env::temp_dir();
        let test_file = temp_dir.join("test_alignment.phy");

        let taxa = vec!["Taxon_A".to_string(), "Taxon_B".to_string()];
        let seqs = vec!["ATGCATGC".to_string(), "ATGCATGT".to_string()];

        write_phylip(&test_file, &taxa, &seqs, false).unwrap();
        let parsed = parse_phylip(&test_file).unwrap();

        assert_eq!(parsed.taxa, taxa);
        assert_eq!(parsed.sequences, seqs);
        assert_eq!(parsed.length, 8);
        assert_eq!(parsed.num_taxa, 2);

        let _ = std::fs::remove_file(test_file);
    }
}
