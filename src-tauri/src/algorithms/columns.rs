/// Port of PhyloProcessR `trimAlignmentColumns.R`
/// Removes alignment columns where the proportion of gap / missing characters meets or exceeds
/// `max_gap_percent` (0.0 - 100.0).
pub fn trim_alignment_columns(
    sequences: &[String],
    max_gap_percent: f64,
    count_n_as_gap: bool,
) -> (Vec<String>, Vec<usize>) {
    if sequences.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let length = sequences[0].len();
    let num_taxa = sequences.len();
    if length == 0 || num_taxa <= 1 {
        return (sequences.to_vec(), Vec::new());
    }

    let mut kept_columns = Vec::new();
    let mut dropped_columns = Vec::new();

    for col in 0..length {
        let mut gap_count = 0usize;
        for seq in sequences {
            if let Some(&b) = seq.as_bytes().get(col) {
                let u = b.to_ascii_uppercase();
                if u == b'-' || u == b'?' || (count_n_as_gap && u == b'N') {
                    gap_count += 1;
                }
            }
        }

        let gap_percent = (gap_count as f64 / num_taxa as f64) * 100.0;
        if gap_percent < max_gap_percent {
            kept_columns.push(col);
        } else {
            dropped_columns.push(col);
        }
    }

    if dropped_columns.is_empty() {
        return (sequences.to_vec(), Vec::new());
    }

    let trimmed_seqs: Vec<String> = sequences
        .iter()
        .map(|s| {
            let bytes = s.as_bytes();
            kept_columns
                .iter()
                .filter_map(|&col| bytes.get(col).copied())
                .map(|b| b as char)
                .collect()
        })
        .collect();

    (trimmed_seqs, dropped_columns)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim_columns() {
        let seqs = vec![
            "A-T-G".to_string(),
            "A-T-G".to_string(),
            "A-T-G".to_string(),
            "A-T-G".to_string(),
        ];
        // Col 0: A (0% gaps) -> keep
        // Col 1: - (100% gaps) -> drop
        // Col 2: T (0% gaps) -> keep
        // Col 3: - (100% gaps) -> drop
        // Col 4: G (0% gaps) -> keep
        let (trimmed, dropped) = trim_alignment_columns(&seqs, 50.0, true);
        assert_eq!(trimmed[0], "ATG");
        assert_eq!(dropped, vec![1, 3]);
    }
}
