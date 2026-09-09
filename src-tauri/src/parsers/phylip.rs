use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use crate::models::{Alignment, AlignmentFormat};
use crate::parsers::validate_alignment_shape;

pub fn parse_phylip<P: AsRef<Path>>(path: P) -> Result<Alignment, String> {
    let path_ref = path.as_ref();
    let content = fs::read_to_string(path_ref)
        .map_err(|e| format!("Failed to read PHYLIP file: {}", e))?;
    let file_name = path_ref
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());
    let id = path_ref
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());
    parse_phylip_str(&content, &id, &file_name, &path_ref.to_string_lossy())
}

/// Parses PHYLIP from memory, for callers with no filesystem.
pub fn parse_phylip_str(
    content: &str,
    id: &str,
    file_name: &str,
    file_path: &str,
) -> Result<Alignment, String> {

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
    let expected_length: usize = parts[1]
        .parse()
        .map_err(|_| format!("Invalid length in PHYLIP header: {}", parts[1]))?;

    let mut taxa: Vec<String> = Vec::with_capacity(expected_taxa);
    let mut sequences: Vec<String> = Vec::with_capacity(expected_taxa);

    let remaining_lines: Vec<&str> = lines.collect();

    // Splits a record line into a taxon name and its sequence characters.
    // Strict PHYLIP puts the name in the first 10 columns with no separator, so
    // fall back to that only when the line holds a single token.
    let split_named_line = |line: &str| -> (String, String) {
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() >= 2 {
            return (tokens[0].to_string(), tokens[1..].concat());
        }
        if line.len() > 10 && line.is_char_boundary(10) {
            let (name, seq) = line.split_at(10);
            return (
                name.trim().to_string(),
                seq.split_whitespace().collect::<String>(),
            );
        }
        (line.trim().to_string(), String::new())
    };

    let sequence_chars = |line: &str| -> String { line.split_whitespace().collect() };

    // A continuation line holds only alignment states. Taxon names in practice
    // contain characters outside this alphabet, which keeps the two apart when
    // a file omits the blank line between interleaved blocks.
    let is_sequence_only = |token: &str| -> bool {
        !token.is_empty()
            && token.chars().all(|c| {
                matches!(
                    c.to_ascii_uppercase(),
                    'A' | 'C' | 'G' | 'T' | 'U' | 'N' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M'
                        | 'B' | 'D' | 'H' | 'V' | 'X' | '-' | '?' | '.' | '*'
                )
            })
    };

    if remaining_lines.len() == expected_taxa {
        // Exactly one line per taxon.
        for line in remaining_lines {
            let (name, seq) = split_named_line(line);
            taxa.push(name);
            sequences.push(seq);
        }
    } else {
        // Sequential records spread over several lines, or interleaved blocks.
        // Both layouts follow one rule: a line that carries only sequence
        // continues the first record that has not reached the declared length.
        // That fills records in order for sequential files and cycles through
        // them for interleaved files, without inventing a taxon from bases.
        let mut fallback_idx = 0;

        for line in remaining_lines {
            let tokens: Vec<&str> = line.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }

            let continues_a_record = tokens.len() == 1
                && is_sequence_only(tokens[0])
                && (!taxa.is_empty() || expected_length == 0);

            if continues_a_record {
                let target = if expected_length > 0 {
                    sequences.iter().position(|seq| seq.len() < expected_length)
                } else {
                    None
                };
                let target = target.unwrap_or_else(|| {
                    let idx = fallback_idx % sequences.len().max(1);
                    fallback_idx += 1;
                    idx
                });
                if let Some(sequence) = sequences.get_mut(target) {
                    sequence.push_str(&sequence_chars(line));
                    continue;
                }
            }

            if taxa.len() < expected_taxa {
                let (name, seq) = split_named_line(line);
                taxa.push(name);
                sequences.push(seq);
                continue;
            }

            // A named line in a later block identifies its own record.
            let named_target = tokens
                .first()
                .and_then(|name| taxa.iter().position(|taxon| taxon == name));
            let target = named_target.unwrap_or_else(|| {
                sequences
                    .iter()
                    .position(|seq| expected_length > 0 && seq.len() < expected_length)
                    .unwrap_or(0)
            });
            let seq_chunk = if named_target.is_some() {
                tokens[1..].concat()
            } else {
                tokens.concat()
            };
            if let Some(sequence) = sequences.get_mut(target) {
                sequence.push_str(&seq_chunk);
            }
        }
    }

    if taxa.is_empty() {
        return Err("No sequences found in PHYLIP file".to_string());
    }
    if taxa.len() != expected_taxa {
        return Err(format!(
            "PHYLIP header declares {} taxa but {} were read",
            expected_taxa,
            taxa.len()
        ));
    }

    validate_alignment_shape(&taxa, &sequences, Some(expected_length), "PHYLIP")?;

    Ok(Alignment::new(
        id.to_string(),
        file_name.to_string(),
        file_path.to_string(),
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

    #[test]
    fn sequential_multi_line_records_do_not_invent_taxa() {
        let temp_dir = std::env::temp_dir();
        let test_file = temp_dir.join("af_sequential_multiline.phy");
        // Two taxa, 20 sites each, written across two lines per record.
        std::fs::write(
            &test_file,
            "2 20\nTaxon_A    ACGTACGTAC\nGTACGTACGT\nTaxon_B    ACGTACGTAC\nGTACGTACGA\n",
        )
        .unwrap();

        let parsed = parse_phylip(&test_file).unwrap();
        assert_eq!(parsed.num_taxa, 2);
        assert_eq!(parsed.taxa, vec!["Taxon_A".to_string(), "Taxon_B".to_string()]);
        assert_eq!(parsed.sequences[0], "ACGTACGTACGTACGTACGT");
        assert_eq!(parsed.sequences[1], "ACGTACGTACGTACGTACGA");
        assert_eq!(parsed.length, 20);

        let _ = std::fs::remove_file(test_file);
    }

    #[test]
    fn interleaved_blocks_still_parse() {
        let temp_dir = std::env::temp_dir();
        let test_file = temp_dir.join("af_interleaved.phy");
        std::fs::write(
            &test_file,
            "2 20\nTaxon_A    ACGTACGTAC\nTaxon_B    ACGTACGTAC\n\nGTACGTACGT\nGTACGTACGA\n",
        )
        .unwrap();

        let parsed = parse_phylip(&test_file).unwrap();
        assert_eq!(parsed.num_taxa, 2);
        assert_eq!(parsed.sequences[0], "ACGTACGTACGTACGTACGT");
        assert_eq!(parsed.sequences[1], "ACGTACGTACGTACGTACGA");

        let _ = std::fs::remove_file(test_file);
    }

    #[test]
    fn ragged_alignment_is_rejected() {
        let temp_dir = std::env::temp_dir();
        let test_file = temp_dir.join("af_ragged.phy");
        std::fs::write(&test_file, "2 8\nTaxon_A    ACGTACGT\nTaxon_B    ACGTAC\n").unwrap();

        let error = parse_phylip(&test_file).unwrap_err();
        assert!(error.contains("not rectangular"), "unexpected error: {error}");

        let _ = std::fs::remove_file(test_file);
    }

    #[test]
    fn multibyte_taxon_name_does_not_panic() {
        let temp_dir = std::env::temp_dir();
        let test_file = temp_dir.join("af_multibyte.phy");
        // A single-token line whose 10th byte falls inside a multi-byte char.
        std::fs::write(&test_file, "1 4\nTaxonZo\u{00e9}\u{00e9}ACGT\n").unwrap();

        // Must return a result or an error, never abort the process.
        let _ = parse_phylip(&test_file);

        let _ = std::fs::remove_file(test_file);
    }
}
