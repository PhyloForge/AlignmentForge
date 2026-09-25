use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use crate::models::{Alignment, AlignmentFormat};
use crate::parsers::{normalize_sequence_symbols, validate_alignment_shape};

fn remove_comments(content: &str) -> Result<String, String> {
    let mut cleaned = String::with_capacity(content.len());
    let mut comment_depth = 0usize;
    let mut quote = None;

    for character in content.chars() {
        if comment_depth > 0 {
            match character {
                '[' => comment_depth += 1,
                ']' => comment_depth -= 1,
                _ => {}
            }
            continue;
        }
        if let Some(delimiter) = quote {
            cleaned.push(character);
            if character == delimiter {
                quote = None;
            }
            continue;
        }
        match character {
            '\'' | '"' => {
                quote = Some(character);
                cleaned.push(character);
            }
            '[' => comment_depth = 1,
            _ => cleaned.push(character),
        }
    }

    if comment_depth > 0 {
        return Err("NEXUS file contains an unterminated comment".to_string());
    }
    Ok(cleaned)
}

/// Finds `keyword` as a whole word in `upper`, at or after byte `from`.
fn find_word(upper: &str, keyword: &str, from: usize) -> Option<usize> {
    let bytes = upper.as_bytes();
    let is_word_byte =
        |index: usize| bytes.get(index).is_some_and(|byte| byte.is_ascii_alphanumeric() || *byte == b'_');
    let mut search_from = from;
    while let Some(offset) = upper.get(search_from..)?.find(keyword) {
        let start = search_from + offset;
        let end = start + keyword.len();
        if (start == 0 || !is_word_byte(start - 1)) && !is_word_byte(end) {
            return Some(start);
        }
        search_from = end;
    }
    None
}

/// Finds the MATRIX command. A command follows a `;`, so a title or a taxon
/// label that contains the word is not a match.
fn find_matrix_command(upper: &str) -> Option<usize> {
    let mut from = 0;
    while let Some(start) = find_word(upper, "MATRIX", from) {
        if upper[..start].trim_end().ends_with(';') {
            return Some(start);
        }
        from = start + "MATRIX".len();
    }
    None
}

/// Reads `KEY=value` from the header. Only a whole word followed by `=` is a
/// key, so a taxon label such as `Singapore_frog` cannot hide `GAP=`.
fn metadata_value(metadata: &str, key: &str) -> Option<String> {
    let upper = metadata.to_ascii_uppercase();
    let bytes = metadata.as_bytes();
    let mut from = 0;
    while let Some(key_start) = find_word(&upper, key, from) {
        from = key_start + key.len();
        let mut index = from;
        while index < bytes.len() && bytes[index].is_ascii_whitespace() {
            index += 1;
        }
        if bytes.get(index) != Some(&b'=') {
            continue;
        }
        index += 1;
        while index < bytes.len() && bytes[index].is_ascii_whitespace() {
            index += 1;
        }
        let end = metadata[index..]
            .find(|character: char| character.is_whitespace() || character == ';')
            .map(|offset| index + offset)
            .unwrap_or(metadata.len());
        return Some(metadata[index..end].trim_matches(['\'', '"']).to_string());
    }
    None
}

fn split_matrix_row(line: &str) -> Result<Vec<String>, String> {
    let mut tokens = Vec::new();
    let mut token = String::new();
    let mut characters = line.chars().peekable();

    while let Some(character) = characters.next() {
        if character.is_whitespace() {
            if !token.is_empty() {
                tokens.push(std::mem::take(&mut token));
            }
            continue;
        }
        if character == '\'' || character == '"' {
            let delimiter = character;
            let mut closed = false;
            while let Some(quoted) = characters.next() {
                if quoted == delimiter {
                    if characters.peek() == Some(&delimiter) {
                        token.push(quoted);
                        characters.next();
                    } else {
                        closed = true;
                        break;
                    }
                } else {
                    token.push(quoted);
                }
            }
            if !closed {
                return Err("NEXUS matrix contains an unterminated quoted taxon name".to_string());
            }
        } else {
            token.push(character);
        }
    }
    if !token.is_empty() {
        tokens.push(token);
    }
    Ok(tokens)
}

pub fn parse_nexus<P: AsRef<Path>>(path: P) -> Result<Alignment, String> {
    let path_ref = path.as_ref();
    let content = fs::read_to_string(path_ref)
        .map_err(|e| format!("Failed to read NEXUS file: {}", e))?;
    let file_name = path_ref
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());
    let id = path_ref
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());
    parse_nexus_str(&content, &id, &file_name, &path_ref.to_string_lossy())
}

/// Parses NEXUS from memory, for callers with no filesystem.
pub fn parse_nexus_str(
    content: &str,
    id: &str,
    file_name: &str,
    file_path: &str,
) -> Result<Alignment, String> {

    let cleaned = remove_comments(content)?;
    let upper = cleaned.to_ascii_uppercase();
    let matrix_start = find_matrix_command(&upper)
        .ok_or_else(|| "No matrix data found in NEXUS file".to_string())?;
    let metadata = &cleaned[..matrix_start];
    let declared_taxa = metadata_value(metadata, "NTAX")
        .map(|value| value.parse::<usize>().map_err(|_| "Invalid NEXUS NTAX value".to_string()))
        .transpose()?;
    let declared_length = metadata_value(metadata, "NCHAR")
        .map(|value| value.parse::<usize>().map_err(|_| "Invalid NEXUS NCHAR value".to_string()))
        .transpose()?;
    let datatype = metadata_value(metadata, "DATATYPE").unwrap_or_else(|| "DNA".to_string());
    if !matches!(datatype.to_ascii_uppercase().as_str(), "DNA" | "RNA" | "NUCLEOTIDE") {
        return Err(format!("Unsupported NEXUS DATATYPE={datatype}; only nucleotide data is supported"));
    }
    let gap = metadata_value(metadata, "GAP").and_then(|value| value.chars().next()).unwrap_or('-');
    let missing = metadata_value(metadata, "MISSING").and_then(|value| value.chars().next()).unwrap_or('?');
    let match_character = metadata_value(metadata, "MATCHCHAR").and_then(|value| value.chars().next());
    // INTERLEAVE may stand alone or carry YES or NO.
    let interleaved = match metadata_value(metadata, "INTERLEAVE") {
        Some(value) => !value.eq_ignore_ascii_case("NO"),
        None => find_word(&metadata.to_ascii_uppercase(), "INTERLEAVE", 0).is_some(),
    };

    let mut taxa: Vec<String> = Vec::new();
    let mut sequences: Vec<String> = Vec::new();
    let mut taxon_index_map: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut continuation_row = 0usize;
    let mut terminated = false;

    for line in cleaned[matrix_start + "MATRIX".len()..].lines() {
        let mut trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(position) = trimmed.find(';') {
            trimmed = trimmed[..position].trim();
            terminated = true;
        }
        if !trimmed.is_empty() {
            let tokens = split_matrix_row(trimmed)?;
            // Without INTERLEAVE, a row continues on the next lines until it
            // holds NCHAR sites, whatever the number of tokens on a line.
            let continues_last_row = !interleaved
                && match declared_length {
                    Some(length) => sequences.last().is_some_and(|sequence| {
                        let added: usize = tokens.iter().map(String::len).sum();
                        sequence.len() < length && sequence.len() + added <= length
                    }),
                    None => tokens.len() == 1 && !taxa.is_empty(),
                };
            if continues_last_row {
                if let Some(sequence) = sequences.last_mut() {
                    sequence.push_str(&tokens.concat());
                }
            } else if interleaved && tokens.len() == 1 && !taxa.is_empty() {
                sequences[continuation_row].push_str(&tokens[0]);
                continuation_row = (continuation_row + 1) % taxa.len();
            } else if let Some(name) = tokens.first() {
                let sequence = tokens[1..].concat();
                if let Some(&index) = taxon_index_map.get(name) {
                    sequences[index].push_str(&sequence);
                } else {
                    let index = taxa.len();
                    taxon_index_map.insert(name.clone(), index);
                    taxa.push(name.clone());
                    sequences.push(sequence);
                }
            } else {
                return Err("Invalid empty row in NEXUS matrix".to_string());
            }
        }
        if terminated {
            break;
        }
    }

    if taxa.is_empty() {
        return Err("No matrix data found in NEXUS file".to_string());
    }
    if !terminated {
        return Err("NEXUS matrix is missing its terminating semicolon".to_string());
    }
    if let Some(expected) = declared_taxa {
        if taxa.len() != expected {
            return Err(format!(
                "NEXUS declares {expected} taxa but {} were read",
                taxa.len()
            ));
        }
    }

    for sequence in &mut sequences {
        *sequence = sequence
            .chars()
            .map(|character| {
                if character == gap {
                    '-'
                } else if character == missing {
                    '?'
                } else {
                    character
                }
            })
            .collect();
    }
    if let Some(match_character) = match_character {
        let reference = sequences[0].clone();
        for (row, sequence) in sequences.iter_mut().enumerate().skip(1) {
            let mut resolved = String::with_capacity(sequence.len());
            for (column, character) in sequence.chars().enumerate() {
                if character == match_character {
                    let reference_character = reference.chars().nth(column).ok_or_else(|| {
                        format!("NEXUS match character in row {} exceeds the first sequence", row + 1)
                    })?;
                    resolved.push(reference_character);
                } else {
                    resolved.push(character);
                }
            }
            *sequence = resolved;
        }
    }

    normalize_sequence_symbols(&mut sequences);
    validate_alignment_shape(&taxa, &sequences, declared_length, "NEXUS")?;

    Ok(Alignment::new(
        id.to_string(),
        file_name.to_string(),
        file_path.to_string(),
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
    // No MATCHCHAR: the writer never writes match characters, and a declared
    // `.` makes readers replace any `.` in the data with the first taxon's base.
    writeln!(
        file,
        "  FORMAT DATATYPE=DNA GAP=- MISSING=?{};",
        if interleaved { " INTERLEAVE" } else { "" }
    )
    .map_err(|e| e.to_string())?;
    writeln!(file, "  MATRIX").map_err(|e| e.to_string())?;

    let quoted_taxa: Vec<String> = taxa
        .iter()
        .map(|name| format!("'{}'", name.replace('\'', "''")))
        .collect();
    let max_name_len = quoted_taxa.iter().map(|t| t.len()).max().unwrap_or(10);
    let pad = (max_name_len + 4).max(12);

    if !interleaved {
        for (name, seq) in quoted_taxa.iter().zip(sequences.iter()) {
            writeln!(file, "    {:<pad$} {}", name, seq, pad = pad).map_err(|e| e.to_string())?;
        }
    } else {
        let chunk_size = 60;
        let total_chunks = (length + chunk_size - 1) / chunk_size;

        for chunk_idx in 0..total_chunks {
            let start = chunk_idx * chunk_size;
            let end = (start + chunk_size).min(length);

            for (name, seq) in quoted_taxa.iter().zip(sequences.iter()) {
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

    #[test]
    fn written_nexus_does_not_declare_a_match_character() {
        let test_file = std::env::temp_dir().join("af_nexus_no_matchchar.nex");
        let taxa = vec!["a".to_string(), "b".to_string()];
        let seqs = vec!["ACGTACGT".to_string(), "ACGTTCGA".to_string()];

        for interleaved in [false, true] {
            write_nexus(&test_file, &taxa, &seqs, interleaved).unwrap();
            let written = std::fs::read_to_string(&test_file).unwrap();
            assert!(!written.contains("MATCHCHAR"), "{written}");
            assert_eq!(parse_nexus(&test_file).unwrap().sequences, seqs);
        }

        let _ = std::fs::remove_file(test_file);
    }

    #[test]
    fn reads_metadata_match_characters_and_quoted_names() {
        let content = "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=4;\nFORMAT DATATYPE=DNA GAP=~ MISSING=X MATCHCHAR=.;\nMATRIX\n'first sample' AT~X\n'second ''sample''' ....;\nEND;\n";
        let parsed = parse_nexus_str(content, "id", "test.nex", "test.nex").unwrap();

        assert_eq!(parsed.taxa, vec!["first sample", "second 'sample'"]);
        assert_eq!(parsed.sequences, vec!["AT-?", "AT-?"]);
    }

    #[test]
    fn parses_sequential_rows_that_span_several_lines() {
        let content = "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=12;\nFORMAT DATATYPE=DNA;\nMATRIX\na ACGT\nACGT\nAC GT\nb\nCCCC GGGG\nTTTT\n;\nEND;\n";
        let parsed = parse_nexus_str(content, "id", "test.nex", "test.nex").unwrap();

        assert_eq!(parsed.taxa, vec!["a", "b"]);
        assert_eq!(parsed.sequences, vec!["ACGTACGTACGT", "CCCCGGGGTTTT"]);
    }

    #[test]
    fn interleaved_rows_without_names_still_rotate() {
        let content = "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=2 NCHAR=8;\nFORMAT DATATYPE=DNA INTERLEAVE=YES;\nMATRIX\na ACGT\nb CCCC\nGGGG\nTTTT\n;\nEND;\n";
        let parsed = parse_nexus_str(content, "id", "test.nex", "test.nex").unwrap();

        assert_eq!(parsed.sequences, vec!["ACGTGGGG", "CCCCTTTT"]);
    }

    #[test]
    fn header_keys_and_matrix_match_whole_words_only() {
        let content = "#NEXUS\nBEGIN TAXA;\nDIMENSIONS NTAX=2;\nTAXLABELS Singapore_frog Matrix_toad;\nEND;\nBEGIN CHARACTERS;\nTITLE 'Frog matrix';\nDIMENSIONS NCHAR=4;\nFORMAT DATATYPE=DNA GAP=~ MISSING=?;\nMATRIX\nSingapore_frog AC~T\nMatrix_toad AC?T\n;\nEND;\n";
        let parsed = parse_nexus_str(content, "id", "test.nex", "test.nex").unwrap();

        assert_eq!(parsed.taxa, vec!["Singapore_frog", "Matrix_toad"]);
        assert_eq!(parsed.sequences, vec!["AC-T", "AC?T"]);
    }

    #[test]
    fn rejects_declared_dimension_mismatch_and_missing_terminator() {
        let mismatch = "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=1 NCHAR=40;\nMATRIX\na ATGC;\nEND;";
        assert!(parse_nexus_str(mismatch, "id", "test.nex", "test.nex").is_err());

        let unterminated = "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=1 NCHAR=4;\nMATRIX\na ATGC\nEND";
        let error = parse_nexus_str(unterminated, "id", "test.nex", "test.nex").unwrap_err();
        assert!(error.contains("terminating semicolon"));
    }
}
