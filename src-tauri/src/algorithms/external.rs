/// Port of PhyloProcessR `trimExternal.R`
/// Trims ragged 5' and 3' exterior columns from an alignment that are covered by fewer than
/// `min_taxa_percent` sequences.
/// If `codon_preserving` is true, snaps both boundaries to preserve complete triplet reading frames.
pub fn trim_external(
    sequences: &[String],
    min_taxa_percent: f64,
    codon_preserving: bool,
) -> (Vec<String>, Vec<usize>, (usize, usize)) {
    if sequences.is_empty() {
        return (Vec::new(), Vec::new(), (0, 0));
    }
    let length = sequences[0].len();
    let num_taxa = sequences.len();
    if length == 0 || num_taxa <= 2 {
        return (sequences.to_vec(), Vec::new(), (0, length));
    }

    let min_taxa_count = ((num_taxa as f64 * (min_taxa_percent / 100.0)).ceil() as usize).max(1);

    // Compute informative count per column
    let mut col_informative_counts = vec![0usize; length];
    for col in 0..length {
        let mut count = 0usize;
        for seq in sequences {
            if let Some(&b) = seq.as_bytes().get(col) {
                let u = b.to_ascii_uppercase();
                if u != b'-' && u != b'?' && u != b'N' {
                    count += 1;
                }
            }
        }
        col_informative_counts[col] = count;
    }

    // Find first column meeting threshold
    let mut start_col = None;
    for (col, &count) in col_informative_counts.iter().enumerate() {
        if count >= min_taxa_count {
            start_col = Some(col);
            break;
        }
    }

    // Find last column meeting threshold
    let mut end_col = None;
    for (col, &count) in col_informative_counts.iter().enumerate().rev() {
        if count >= min_taxa_count {
            end_col = Some(col + 1); // 1-past-the-end
            break;
        }
    }

    let (mut start, mut end) = match (start_col, end_col) {
        (Some(s), Some(e)) if s < e => (s, e),
        _ => return (vec![String::new(); num_taxa], (0..length).collect(), (0, 0)),
    };

    // Codon-preserving frame adjustment
    if codon_preserving {
        let rem = start % 3;
        if rem == 1 {
            start += 2;
        } else if rem == 2 {
            start += 1;
        }
        if start >= end {
            return (vec![String::new(); num_taxa], (0..length).collect(), (0, 0));
        }

        // Keep the retained interval divisible by three after shifting the start.
        end -= (end - start) % 3;
        if start >= end {
            return (vec![String::new(); num_taxa], (0..length).collect(), (0, 0));
        }
    }

    let mut trimmed_columns = Vec::new();
    for col in 0..start {
        trimmed_columns.push(col);
    }
    for col in end..length {
        trimmed_columns.push(col);
    }

    let trimmed_seqs: Vec<String> = sequences
        .iter()
        .map(|s| {
            if start < s.len() && end <= s.len() {
                s[start..end].to_string()
            } else {
                String::new()
            }
        })
        .collect();

    (trimmed_seqs, trimmed_columns, (start, end))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_external() {
        let seqs = vec![
            "--ATGCATGC--".to_string(), // 12 bp
            "--ATGCATGCAT".to_string(),
            "ATATGCATGC--".to_string(),
            "--ATGCATGC--".to_string(),
        ];
        // min 75% = 3 taxa required. Cols 0,1 only have 1 taxon. Cols 2..9 have 4 taxa. Cols 10,11 have 1 taxon.
        let (trimmed, dropped_cols, (s, e)) = trim_external(&seqs, 75.0, false);
        assert_eq!(s, 2);
        assert_eq!(e, 10);
        assert_eq!(trimmed[0], "ATGCATGC");
        assert_eq!(dropped_cols, vec![0, 1, 10, 11]);
    }

    #[test]
    fn test_trim_external_codon() {
        let seqs = vec![
            "-ATGCATGCATGC".to_string(),
            "-ATGCATGCATGC".to_string(),
            "-ATGCATGCATGC".to_string(),
        ];
        // start is 1, so 1 % 3 = 1 -> adjust forward by 2 to index 3
        let (trimmed, _, (s, e)) = trim_external(&seqs, 50.0, true);
        assert_eq!(s, 3);
        assert_eq!(e, 12);
        assert_eq!(trimmed[0], "GCATGCATG");
        assert_eq!(trimmed[0].len() % 3, 0);
    }
}
