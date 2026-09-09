use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use crate::models::{Alignment, AlignmentFormat};
use crate::parsers::validate_alignment_shape;

/// High-performance zero-copy byte streaming FASTA parser.
/// Reads the entire file into memory in a single syscall and parses headers and sequences
/// without intermediate line allocations.
/// Turns collected sequence bytes into text, naming the offending record rather
/// than silently yielding an empty sequence for one bad byte.
fn decode_sequence(bytes: Vec<u8>, taxon: &str) -> Result<String, String> {
    String::from_utf8(bytes).map_err(|error| {
        format!(
            "Sequence for '{}' contains a byte that is not valid text at offset {}",
            if taxon.is_empty() { "<unnamed record>" } else { taxon },
            error.utf8_error().valid_up_to()
        )
    })
}

pub fn parse_fasta<P: AsRef<Path>>(path: P) -> Result<Alignment, String> {
    let path_ref = path.as_ref();
    let content = fs::read(path_ref).map_err(|e| format!("Failed to read FASTA file: {}", e))?;
    let file_name = path_ref
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());
    let id = path_ref
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| "unknown".to_string());
    parse_fasta_bytes(&content, &id, &file_name, &path_ref.to_string_lossy())
}

/// Parses FASTA from memory. The file reading wrapper above calls this, and so
/// does the browser build, which has no filesystem.
pub fn parse_fasta_bytes(
    content: &[u8],
    id: &str,
    file_name: &str,
    file_path: &str,
) -> Result<Alignment, String> {

    let mut taxa = Vec::new();
    let mut sequences = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq_bytes = Vec::new();
    let mut header_for_error = String::new();

    let mut i = 0;
    let len = content.len();

    while i < len {
        // Find line end
        let mut line_end = i;
        while line_end < len && content[line_end] != b'\n' && content[line_end] != b'\r' {
            line_end += 1;
        }

        let line = &content[i..line_end];
        // Advance i past line break
        i = line_end;
        if i < len && content[i] == b'\r' {
            i += 1;
        }
        if i < len && content[i] == b'\n' {
            i += 1;
        }

        // Trim leading whitespace
        let mut start = 0;
        while start < line.len() && (line[start] == b' ' || line[start] == b'\t') {
            start += 1;
        }
        let trimmed = &line[start..];

        if trimmed.is_empty() {
            continue;
        }

        // ';' introduces a comment in the original FASTA definition. Treating it
        // as a record start invents a taxon, so skip the line instead.
        if trimmed[0] == b';' {
            continue;
        }

        if trimmed[0] == b'>' {
            if let Some(header) = current_header.take() {
                taxa.push(header);
                sequences.push(decode_sequence(
                    std::mem::take(&mut current_seq_bytes),
                    &header_for_error,
                )?);
            }
            let header_slice = &trimmed[1..];
            let header_str = String::from_utf8_lossy(header_slice).trim().to_string();
            header_for_error = header_str.clone();
            current_header = Some(header_str);
        } else {
            // Filter out internal whitespace and append directly
            for &b in trimmed {
                if b != b' ' && b != b'\t' {
                    current_seq_bytes.push(b);
                }
            }
        }
    }

    if let Some(header) = current_header {
        taxa.push(header);
        sequences.push(decode_sequence(current_seq_bytes, &header_for_error)?);
    }

    if taxa.is_empty() {
        return Err("No sequences found in FASTA file".to_string());
    }

    validate_alignment_shape(&taxa, &sequences, None, "FASTA")?;

    Ok(Alignment::new(
        id.to_string(),
        file_name.to_string(),
        file_path.to_string(),
        AlignmentFormat::Fasta,
        taxa,
        sequences,
    ))
}

pub fn write_fasta<P: AsRef<Path>>(
    path: P,
    taxa: &[String],
    sequences: &[String],
    line_width: Option<usize>,
) -> Result<(), String> {
    let mut file = File::create(path).map_err(|e| format!("Failed to create FASTA: {}", e))?;

    for (header, seq) in taxa.iter().zip(sequences.iter()) {
        writeln!(file, ">{}", header).map_err(|e| e.to_string())?;
        match line_width {
            Some(w) if w > 0 => {
                for chunk in seq.as_bytes().chunks(w) {
                    let chunk_str = std::str::from_utf8(chunk).unwrap_or("");
                    writeln!(file, "{}", chunk_str).map_err(|e| e.to_string())?;
                }
            }
            _ => {
                writeln!(file, "{}", seq).map_err(|e| e.to_string())?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fasta_roundtrip() {
        let temp_dir = std::env::temp_dir();
        let test_file = temp_dir.join("test_alignment.fasta");

        let taxa = vec!["Taxon_A".to_string(), "Taxon_B".to_string()];
        let seqs = vec!["ATGC-ATGC".to_string(), "ATGC-ATGT".to_string()];

        write_fasta(&test_file, &taxa, &seqs, Some(80)).unwrap();
        let parsed = parse_fasta(&test_file).unwrap();

        assert_eq!(parsed.taxa, taxa);
        assert_eq!(parsed.sequences, seqs);
        assert_eq!(parsed.length, 9);
        assert_eq!(parsed.num_taxa, 2);

        let _ = std::fs::remove_file(test_file);
    }

    #[test]
    fn semicolon_comment_lines_do_not_become_taxa() {
        let temp_dir = std::env::temp_dir();
        let test_file = temp_dir.join("af_comment.fa");
        std::fs::write(
            &test_file,
            "; a legacy FASTA comment\n>Taxon_A\nACGT\n>Taxon_B\nACGA\n",
        )
        .unwrap();

        let parsed = parse_fasta(&test_file).unwrap();
        assert_eq!(parsed.num_taxa, 2);
        assert_eq!(parsed.taxa, vec!["Taxon_A".to_string(), "Taxon_B".to_string()]);

        let _ = std::fs::remove_file(test_file);
    }

    #[test]
    fn invalid_utf8_reports_the_taxon_instead_of_emptying_it() {
        let temp_dir = std::env::temp_dir();
        let test_file = temp_dir.join("af_bad_utf8.fa");
        let mut bytes = b">Taxon_A\nACG".to_vec();
        bytes.push(0xFF);
        bytes.extend_from_slice(b"T\n");
        std::fs::write(&test_file, bytes).unwrap();

        let error = parse_fasta(&test_file).unwrap_err();
        assert!(error.contains("Taxon_A"), "unexpected error: {error}");

        let _ = std::fs::remove_file(test_file);
    }

    #[test]
    fn ragged_fasta_is_rejected() {
        let temp_dir = std::env::temp_dir();
        let test_file = temp_dir.join("af_ragged_fasta.fa");
        std::fs::write(&test_file, ">Taxon_A\nACGTACGT\n>Taxon_B\nACGT\n").unwrap();

        let error = parse_fasta(&test_file).unwrap_err();
        assert!(error.contains("not rectangular"), "unexpected error: {error}");

        let _ = std::fs::remove_file(test_file);
    }
}
