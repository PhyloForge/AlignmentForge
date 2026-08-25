use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum RelativeWidth {
    Alignment,
    Sample,
}

/// Port of PhyloProcessR `trimSampleCoverage.R`
/// Removes samples whose sequence coverage falls below minimum thresholds:
/// 1. `min_coverage_bp`: Minimum absolute non-gap base pairs.
/// 2. `min_coverage_percent`: Minimum non-gap bases relative to total alignment length or longest sample.
pub fn filter_sample_coverage(
    taxa: &[String],
    sequences: &[String],
    min_coverage_bp: usize,
    min_coverage_percent: f64,
    relative_width: RelativeWidth,
) -> (Vec<String>, Vec<String>, Vec<String>) {
    if sequences.len() <= 2 || sequences.is_empty() {
        return (taxa.to_vec(), sequences.to_vec(), Vec::new());
    }

    let alignment_len = sequences[0].len();
    if min_coverage_bp >= alignment_len && alignment_len > 0 {
        return (taxa.to_vec(), sequences.to_vec(), Vec::new());
    }

    // Count non-gap bases per sequence
    let base_counts: Vec<usize> = sequences
        .iter()
        .map(|s| {
            s.as_bytes()
                .iter()
                .filter(|&&b| b != b'-' && b != b'?' && b != b'N' && b != b'n')
                .count()
        })
        .collect();

    let max_sample_bp = base_counts.iter().copied().max().unwrap_or(1).max(1);

    let reference_width = match relative_width {
        RelativeWidth::Alignment => alignment_len.max(1),
        RelativeWidth::Sample => max_sample_bp,
    };

    let mut kept_taxa = Vec::new();
    let mut kept_seqs = Vec::new();
    let mut dropped_taxa = Vec::new();

    for i in 0..taxa.len() {
        let count = base_counts[i];
        let pct = (count as f64 / reference_width as f64) * 100.0;

        if count >= min_coverage_bp && pct >= min_coverage_percent {
            kept_taxa.push(taxa[i].clone());
            kept_seqs.push(sequences[i].clone());
        } else {
            dropped_taxa.push(taxa[i].clone());
        }
    }

    (kept_taxa, kept_seqs, dropped_taxa)
}

/// Port of PhyloProcessR `trimAlignmentRows.R`
/// Removes individual sequences where overall gap percentage meets or exceeds `max_gap_percent`.
pub fn filter_alignment_rows(
    taxa: &[String],
    sequences: &[String],
    max_gap_percent: f64,
    count_n_as_gap: bool,
) -> (Vec<String>, Vec<String>, Vec<String>) {
    if sequences.len() <= 2 || sequences.is_empty() {
        return (taxa.to_vec(), sequences.to_vec(), Vec::new());
    }

    let mut kept_taxa = Vec::new();
    let mut kept_seqs = Vec::new();
    let mut dropped_taxa = Vec::new();

    for (name, seq) in taxa.iter().zip(sequences.iter()) {
        let total = seq.len();
        if total == 0 {
            dropped_taxa.push(name.clone());
            continue;
        }

        let gaps = seq
            .as_bytes()
            .iter()
            .filter(|&&b| {
                let u = b.to_ascii_uppercase();
                u == b'-' || u == b'?' || (count_n_as_gap && u == b'N')
            })
            .count();

        let gap_pct = (gaps as f64 / total as f64) * 100.0;
        if gap_pct < max_gap_percent {
            kept_taxa.push(name.clone());
            kept_seqs.push(seq.clone());
        } else {
            dropped_taxa.push(name.clone());
        }
    }

    (kept_taxa, kept_seqs, dropped_taxa)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_coverage_trimming() {
        let taxa = vec![
            "Taxon_Good".to_string(),
            "Taxon_Short".to_string(),
            "Taxon_Empty".to_string(),
        ];
        let seqs = vec![
            "ATGCATGCATGC".to_string(), // 12 bp
            "AT----------".to_string(), // 2 bp
            "------------".to_string(), // 0 bp
        ];

        let (kept_t, kept_s, dropped) = filter_sample_coverage(
            &taxa,
            &seqs,
            6,
            50.0,
            RelativeWidth::Alignment,
        );
        assert_eq!(kept_t, vec!["Taxon_Good".to_string()]);
        assert_eq!(dropped, vec!["Taxon_Short".to_string(), "Taxon_Empty".to_string()]);
        assert_eq!(kept_s.len(), 1);
    }
}
