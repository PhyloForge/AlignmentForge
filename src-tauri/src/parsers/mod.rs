pub mod fasta;
pub mod nexus;
pub mod phylip;

use std::path::Path;
use std::collections::HashSet;
use crate::models::{Alignment, AlignmentFormat};

/// Converts symbols that have no base meaning in nucleotide data. `.` and `*`
/// become gaps, and `X` becomes `N`, so every filter reads them the same way.
/// NEXUS resolves its declared GAP, MISSING, and MATCHCHAR symbols first.
pub fn normalize_sequence_symbols(sequences: &mut [String]) {
    for sequence in sequences.iter_mut() {
        if !sequence.bytes().any(|byte| matches!(byte, b'.' | b'*' | b'X' | b'x')) {
            continue;
        }
        *sequence = sequence
            .chars()
            .map(|character| match character {
                '.' | '*' => '-',
                'X' => 'N',
                'x' => 'n',
                other => other,
            })
            .collect();
    }
}

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

    let mut seen_taxa = HashSet::with_capacity(taxa.len());
    for taxon in taxa {
        if taxon.trim().is_empty() {
            return Err(format!("{format_name} file contains an empty taxon name"));
        }
        if !seen_taxa.insert(taxon) {
            return Err(format!(
                "{format_name} file contains duplicate taxon name '{taxon}'"
            ));
        }
    }

    let Some(first_length) = sequences.first().map(String::len) else {
        return Ok(());
    };

    for (taxon, sequence) in taxa.iter().zip(sequences.iter()) {
        if sequence.is_empty() {
            return Err(format!(
                "{format_name} record '{taxon}' has an empty sequence"
            ));
        }
        for (position, byte) in sequence.bytes().enumerate() {
            let valid = matches!(
                byte.to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T' | b'U' | b'N' | b'R' | b'Y' | b'S' | b'W'
                    | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' | b'-' | b'?' | b'!'
            );
            if !byte.is_ascii() || !valid {
                let character = sequence
                    .get(position..)
                    .and_then(|rest| rest.chars().next())
                    .unwrap_or(char::REPLACEMENT_CHARACTER);
                return Err(format!(
                    "{format_name} record '{taxon}' contains unsupported character '{character}' at position {}; only nucleotide alignments are supported",
                    position + 1
                ));
            }
        }
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

fn format_from_extension(extension: &str) -> Option<AlignmentFormat> {
    match extension.to_ascii_lowercase().as_str() {
        "fa" | "fasta" | "fna" | "ffn" => Some(AlignmentFormat::Fasta),
        "phy" | "phylip" => Some(AlignmentFormat::Phylip),
        "nex" | "nexus" => Some(AlignmentFormat::Nexus),
        _ => None,
    }
}

/// Detects the format of a file with another extension (`.aln`, `.txt`) from
/// its first non-empty line. The desktop and browser builds both use this rule.
fn format_from_first_line(line: Option<&str>) -> AlignmentFormat {
    let Some(trimmed) = line.map(str::trim) else {
        return AlignmentFormat::Fasta;
    };
    if trimmed.to_ascii_uppercase().starts_with("#NEXUS") {
        return AlignmentFormat::Nexus;
    }
    let parts: Vec<&str> = trimmed.split_whitespace().collect();
    if parts.len() == 2 && parts[0].parse::<usize>().is_ok() && parts[1].parse::<usize>().is_ok() {
        return AlignmentFormat::Phylip;
    }
    AlignmentFormat::Fasta
}

pub fn detect_format<P: AsRef<Path>>(path: P) -> Result<AlignmentFormat, String> {
    let path_ref = path.as_ref();
    let extension = path_ref.extension().and_then(|s| s.to_str()).unwrap_or_default();
    if let Some(format) = format_from_extension(extension) {
        return Ok(format);
    }

    let file = std::fs::File::open(path_ref)
        .map_err(|e| format!("Failed to open file for format detection: {}", e))?;
    use std::io::BufRead;
    let first_line = std::io::BufReader::new(file)
        .lines()
        .map_while(Result::ok)
        .find(|line| !line.trim().is_empty());
    Ok(format_from_first_line(first_line.as_deref()))
}

/// Detects a format from the file name and, failing that, from the content.
/// The browser build has no filesystem, so detection must work from memory.
pub fn detect_format_from_text(file_name: &str, content: &str) -> AlignmentFormat {
    let extension = file_name.rsplit_once('.').map_or("", |(_, extension)| extension);
    if let Some(format) = format_from_extension(extension) {
        return format;
    }
    format_from_first_line(content.lines().find(|line| !line.trim().is_empty()))
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
    use crate::models::AlignmentFormat;

    #[test]
    fn dots_and_stars_become_gaps_and_x_becomes_n() {
        let content = ">a\nACGTACGT\n>b\nACGT.CGA\n>c\nACxT*CGX\n";
        let parsed = super::parse_alignment_text(content, "id", "test.fasta", "test.fasta").unwrap();

        assert_eq!(parsed.sequences, vec!["ACGTACGT", "ACGT-CGA", "ACnT-CGN"]);
    }

    #[test]
    fn unsupported_character_is_named_in_the_error() {
        let content = ">a\nMKLV\n>b\nMKLV\n";
        let error = super::parse_alignment_text(content, "id", "test.fasta", "test.fasta").unwrap_err();

        assert!(error.contains("'L' at position 3"), "{error}");
    }

    #[test]
    fn other_extensions_use_the_first_non_empty_line() {
        use super::detect_format_from_text;
        assert_eq!(detect_format_from_text("locus.txt", "\n\n2 8\na ACGT\n"), AlignmentFormat::Phylip);
        assert_eq!(detect_format_from_text("locus.aln", "#nexus\n"), AlignmentFormat::Nexus);
        // A later line that looks like a PHYLIP header does not change the result.
        assert_eq!(detect_format_from_text("notes.txt", "some text\n2 8\n"), AlignmentFormat::Fasta);
    }

    #[test]
    fn one_alignment_reads_the_same_from_every_format() {
        use crate::pipeline::engine::apply_recipe;
        use crate::pipeline::recipe::TrimmingRecipe;

        let fasta = ">Rana_one\nACGTACGT\nACGT\n>Rana_two\nACGTTCGT-CGA\n>Rana_three\nACG?ACGTACRT\n";
        let sequential_phylip =
            "3 12\nRana_one ACGTACGT\nACGT\nRana_two ACGTT CGT-C GA\nRana_three ACG?ACGTACRT\n";
        let interleaved_phylip =
            "3 12\nRana_one   ACGTAC\nRana_two   ACGTTC\nRana_three ACG?AC\n\nGTACGT\nGT-CGA\nGTACRT\n";
        // Match characters, a declared gap symbol, a comment, and a quoted name.
        let nexus = "#NEXUS\n[three frogs]\nBEGIN DATA;\nDIMENSIONS NTAX=3 NCHAR=12;\n\
            FORMAT DATATYPE=DNA GAP=~ MISSING=? MATCHCHAR=. INTERLEAVE;\nMATRIX\n\
            Rana_one ACGTAC\nRana_two ....T.\n'Rana_three' ...?..\n\n\
            Rana_one GTACGT\nRana_two ..~..A\n'Rana_three' ....R.\n;\nEND;\n";

        let read = |content: &str, file_name: &str| {
            super::parse_alignment_text(content, "locus", file_name, file_name).unwrap()
        };
        let expected = read(fasta, "locus.fasta");
        assert_eq!(expected.taxa, vec!["Rana_one", "Rana_two", "Rana_three"]);
        assert_eq!(expected.sequences, vec!["ACGTACGTACGT", "ACGTTCGT-CGA", "ACG?ACGTACRT"]);

        let recipe = TrimmingRecipe { trim_coverage: false, ..TrimmingRecipe::default() };
        let (expected_output, expected_diff) = apply_recipe(&expected, &recipe, 0);
        assert!(!expected_output.sequences.is_empty());

        for (content, file_name) in [
            (sequential_phylip, "sequential.phy"),
            (interleaved_phylip, "interleaved.phy"),
            (nexus, "locus.nex"),
        ] {
            let alignment = read(content, file_name);
            assert_eq!(alignment.taxa, expected.taxa, "{file_name}");
            assert_eq!(alignment.sequences, expected.sequences, "{file_name}");
            // The same matrix gives the same filter outcome, whatever its format.
            let (output, diff) = apply_recipe(&alignment, &recipe, 0);
            assert_eq!(output.sequences, expected_output.sequences, "{file_name}");
            assert_eq!(
                serde_json::to_value(&diff).unwrap(),
                serde_json::to_value(&expected_diff).unwrap(),
                "{file_name}"
            );
        }
    }

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
