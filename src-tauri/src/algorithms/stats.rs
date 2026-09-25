/// Fast alignment summary statistics, consensus calculation, and genetic distances.

pub fn calculate_gap_stats(sequences: &[String]) -> (usize, usize, f64) {
    if sequences.is_empty() {
        return (0, 0, 0.0);
    }
    let mut gap_count = 0usize;
    let mut total_chars = 0usize;

    for seq in sequences {
        for &b in seq.as_bytes() {
            total_chars += 1;
            match b {
                b'-' | b'?' | b'N' | b'n' | b'!' => gap_count += 1,
                _ => {}
            }
        }
    }

    let gap_percent = if total_chars > 0 {
        (gap_count as f64 / total_chars as f64) * 100.0
    } else {
        0.0
    };

    (gap_count, total_chars, gap_percent)
}

pub fn calculate_gc_percent(sequences: &[String]) -> f64 {
    let mut gc_count = 0usize;
    let mut non_gap_count = 0usize;

    for seq in sequences {
        for &b in seq.as_bytes() {
            match b.to_ascii_uppercase() {
                b'G' | b'C' => {
                    gc_count += 1;
                    non_gap_count += 1;
                }
                b'A' | b'T' | b'U' => {
                    non_gap_count += 1;
                }
                _ => {}
            }
        }
    }

    if non_gap_count > 0 {
        (gc_count as f64 / non_gap_count as f64) * 100.0
    } else {
        0.0
    }
}

/// Computes majority-rule consensus sequence from alignment (zero-allocation stack count)
pub fn compute_majority_consensus(sequences: &[String], remove_gaps: bool) -> String {
    if sequences.is_empty() {
        return String::new();
    }
    let length = sequences[0].len();
    let mut consensus = Vec::with_capacity(length);
    let mut char_counts = [0usize; 128];

    for col in 0..length {
        char_counts.fill(0);
        let mut total_informative = 0usize;

        for seq in sequences {
            if let Some(&b) = seq.as_bytes().get(col) {
                let upper = (b as char).to_ascii_uppercase() as usize;
                if upper < 128 && upper != b'-' as usize && upper != b'?' as usize && upper != b'N' as usize {
                    char_counts[upper] += 1;
                    total_informative += 1;
                }
            }
        }

        if total_informative == 0 {
            if !remove_gaps {
                consensus.push(b'-');
            }
        } else {
            // Find most frequent character
            let mut max_char = b'-';
            let mut max_cnt = 0usize;
            for (idx, &cnt) in char_counts.iter().enumerate() {
                if cnt > max_cnt {
                    max_cnt = cnt;
                    max_char = idx as u8;
                }
            }
            consensus.push(max_char);
        }
    }

    String::from_utf8(consensus).unwrap_or_default()
}

/// Computes pairwise genetic distance between each sequence in an alignment and a target sequence.
/// Matches PhyloProcessR `pairwiseDistanceTarget.R`:
/// Distance is calculated as the proportion of sites that differ between sample and reference,
/// counting only positions where BOTH sequences have non-gap, non-N, non-? data.
pub fn pairwise_distance_to_target(sequences: &[String], target: &str) -> Vec<f64> {
    let target_bytes = target.as_bytes();
    let length = target_bytes.len();

    sequences
        .iter()
        .map(|seq| {
            let seq_bytes = seq.as_bytes();
            let mut mismatches = 0usize;
            let mut overlap_informative = 0usize;

            for col in 0..length {
                let b_target = target_bytes[col];
                let b_seq = if col < seq_bytes.len() {
                    seq_bytes[col]
                } else {
                    b'-'
                };

                // Exclude positions where either is -, ?, N
                let is_invalid_target =
                    b_target == b'-' || b_target == b'?' || b_target == b'N' || b_target == b'n';
                let is_invalid_seq =
                    b_seq == b'-' || b_seq == b'?' || b_seq == b'N' || b_seq == b'n';

                if !is_invalid_target && !is_invalid_seq {
                    overlap_informative += 1;
                    if b_target.to_ascii_uppercase() != b_seq.to_ascii_uppercase() {
                        mismatches += 1;
                    }
                }
            }

            if overlap_informative > 0 {
                mismatches as f64 / overlap_informative as f64
            } else {
                0.0
            }
        })
        .collect()
}

/// Mean distance of the sequences from their majority consensus. The value is
/// approximate: with 40 or more sequences, it uses every (n / 20)th sequence,
/// which is 20 to 30 of them. The consensus uses every sequence. Only the
/// display and the statistics CSV use this value; no filter reads it.
pub fn compute_mean_divergence(sequences: &[String]) -> f64 {
    if sequences.len() <= 1 {
        return 0.0;
    }
    let consensus = compute_majority_consensus(sequences, false);

    // The sample keeps this display value cheap on large loci (numbers in HANDOFF.md).
    let dists = if sequences.len() > 20 {
        let step = sequences.len() / 20;
        let sample_seqs: Vec<String> = sequences.iter().step_by(step).cloned().collect();
        pairwise_distance_to_target(&sample_seqs, &consensus)
    } else {
        pairwise_distance_to_target(sequences, &consensus)
    };

    let sum: f64 = dists.iter().sum();
    sum / dists.len().max(1) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gap_stats() {
        let seqs = vec![
            "AAAA----".to_string(),
            "AAAA----".to_string(),
            "AAAAAAAA".to_string(),
        ];
        let (gap_count, total_chars, gap_pct) = calculate_gap_stats(&seqs);
        assert_eq!(gap_count, 8);
        assert_eq!(total_chars, 24);
        assert!((gap_pct - 33.333).abs() < 0.01);
    }

    #[test]
    fn consensus_uses_every_sequence() {
        let mut sequences = Vec::new();
        for index in 0..100 {
            if index < 80 && index % 2 == 0 {
                sequences.push("CCCC".to_string());
            } else {
                sequences.push("AAAA".to_string());
            }
        }

        assert_eq!(compute_majority_consensus(&sequences, false), "AAAA");
        sequences.reverse();
        assert_eq!(compute_majority_consensus(&sequences, false), "AAAA");
    }

    #[test]
    fn test_consensus_and_distance() {
        let seqs = vec![
            "AAAA".to_string(),
            "AAAT".to_string(),
            "AAAG".to_string(),
            "TTTT".to_string(),
        ];
        let consensus = compute_majority_consensus(&seqs, false);
        assert_eq!(consensus, "AAAT");

        let dists = pairwise_distance_to_target(&seqs, &consensus);
        assert_eq!(dists[0], 0.25);
        assert_eq!(dists[1], 0.0);
        assert_eq!(dists[2], 0.25);
        assert_eq!(dists[3], 0.75);
    }
}
