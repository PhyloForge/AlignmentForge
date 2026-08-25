use crate::algorithms::stats::{compute_majority_consensus, pairwise_distance_to_target};
use crate::models::MaskedSegment;

/// Port of PhyloProcessR `trimSampleSegments.R`
/// Masks localized divergent segments within individual sequences by sliding a window of `window_size`
/// base pairs across the alignment. For each window, if a sample's divergence from consensus exceeds
/// `threshold`, that localized segment is masked as gaps (`-`).
pub fn mask_divergent_segments(
    taxa: &[String],
    sequences: &[String],
    window_size: usize,
    threshold: f64,
) -> (Vec<String>, Vec<MaskedSegment>) {
    if sequences.len() <= 2 || sequences.is_empty() || window_size == 0 {
        return (sequences.to_vec(), Vec::new());
    }

    let length = sequences[0].len();
    if length == 0 {
        return (sequences.to_vec(), Vec::new());
    }

    let _num_seqs = sequences.len();
    let mut seq_chars: Vec<Vec<char>> = sequences.iter().map(|s| s.chars().collect()).collect();
    let mut masked_segments = Vec::new();

    let num_slices = (length + window_size - 1) / window_size;

    for slice_idx in 0..num_slices {
        let start = slice_idx * window_size;
        let end = (start + window_size).min(length);
        if start >= end {
            break;
        }

        // Extract slice sub-sequences
        let slice_seqs: Vec<String> = seq_chars
            .iter()
            .map(|chars| chars[start..end].iter().collect())
            .collect();

        // Majority consensus of this slice
        let slice_consensus = compute_majority_consensus(&slice_seqs, false);
        let dists = pairwise_distance_to_target(&slice_seqs, &slice_consensus);

        for (i, dist) in dists.iter().enumerate() {
            if *dist >= threshold {
                // Mask this slice with '-'
                for col in start..end {
                    seq_chars[i][col] = '-';
                }
                masked_segments.push(MaskedSegment {
                    taxon: taxa.get(i).cloned().unwrap_or_default(),
                    start,
                    end,
                });
            }
        }
    }

    let output_seqs = seq_chars.into_iter().map(|c| c.into_iter().collect()).collect();
    (output_seqs, masked_segments)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_segment_masking() {
        let taxa = vec![
            "Taxon_1".to_string(),
            "Taxon_2".to_string(),
            "Taxon_3".to_string(),
        ];
        // 20 bp alignment: first 10 bp identical, second 10 bp Taxon_3 is wildly divergent
        let seqs = vec![
            "AAAAAAAAAATTTTTTTTTT".to_string(),
            "AAAAAAAAAATTTTTTTTTT".to_string(),
            "AAAAAAAAAACCCCCCCCCC".to_string(),
        ];

        let (masked, records) = mask_divergent_segments(&taxa, &seqs, 10, 0.40);
        assert_eq!(masked[0], "AAAAAAAAAATTTTTTTTTT");
        assert_eq!(masked[1], "AAAAAAAAAATTTTTTTTTT");
        assert_eq!(masked[2], "AAAAAAAAAA----------"); // Second window masked
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].taxon, "Taxon_3");
        assert_eq!(records[0].start, 10);
        assert_eq!(records[0].end, 20);
    }
}
