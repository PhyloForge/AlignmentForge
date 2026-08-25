use crate::models::MaskedSegment;

/// Port of TAPIR / Profile HMM alignment segment cleaner.
///
/// Builds a position-specific Profile HMM model with leave-one-out jackknife
/// emission probabilities (preventing an aberrant sample from inflating its own match profile).
///
/// For each sequence, calculates the per-residue posterior match confidence.
/// A sliding window average of the posterior is computed to prevent random matches
/// from breaking the contiguous bad segment. Windows whose average confidence falls
/// below `min_posterior` are masked.
pub fn clean_with_profile_hmm(
    taxa: &[String],
    sequences: &[String],
    min_posterior: f64,
    min_segment_length: usize,
    min_island_length: usize,
) -> (Vec<String>, Vec<MaskedSegment>) {
    if sequences.is_empty() || taxa.is_empty() {
        return (sequences.to_vec(), Vec::new());
    }

    let length = sequences[0].len();
    let num_seqs = sequences.len();

    if length == 0 || num_seqs < 3 {
        return (sequences.to_vec(), Vec::new());
    }

    let alpha = 0.5f64;
    let null_prob = 0.25f64;

    // 1. Precompute total column base counts (Length x 4: A, C, G, T)
    let mut col_counts: Vec<[f64; 4]> = vec![[0.0; 4]; length];
    let mut col_totals: Vec<f64> = vec![0.0; length];

    for col in 0..length {
        for seq in sequences {
            if let Some(&b) = seq.as_bytes().get(col) {
                let idx = match b.to_ascii_uppercase() {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' | b'U' => 3,
                    _ => 4, // gap or ambiguous
                };
                if idx < 4 {
                    col_counts[col][idx] += 1.0;
                    col_totals[col] += 1.0;
                }
            }
        }
    }

    let mut cleaned_sequences = Vec::with_capacity(num_seqs);
    let mut masked_segments = Vec::new();
    let half_w = min_segment_length / 2;

    for (taxon_idx, seq) in sequences.iter().enumerate() {
        let seq_bytes = seq.as_bytes();
        let mut chars: Vec<char> = seq.chars().collect();
        let mut confidences = vec![1.0f64; length];

        // 2. Score each sequence using Leave-One-Out (Jackknife) Profile HMM posterior
        for col in 0..length {
            let b = seq_bytes.get(col).copied().unwrap_or(b'-');
            let idx = match b.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' | b'U' => 3,
                _ => 4,
            };

            if idx < 4 {
                let count_without_i = (col_counts[col][idx] - 1.0).max(0.0);
                let total_without_i = (col_totals[col] - 1.0).max(0.0);

                let denom = total_without_i + alpha;
                let p_emit = (count_without_i + alpha * null_prob) / denom;

                let conf = p_emit / (p_emit + null_prob);
                confidences[col] = conf;
            } else {
                confidences[col] = 1.0; // Gaps are neutral
            }
        }

        // 3. Smooth confidences with sliding window (O(N) with cumsum)
        let mut cum_conf = vec![0.0f64; length + 1];
        for col in 0..length {
            cum_conf[col + 1] = cum_conf[col] + confidences[col];
        }

        let mut mask = vec![false; length];
        
        for col in 0..length {
            let start = col.saturating_sub(half_w);
            let end = (col + half_w).min(length.saturating_sub(1));
            
            let sum = cum_conf[end + 1] - cum_conf[start];
            let window_len = (end - start + 1) as f64;
            let avg_conf = sum / window_len;

            if avg_conf < min_posterior {
                mask[col] = true;
            }
        }

        
        // Bridge small unmasked islands
        if min_island_length > 0 {
            let mut start_false = None;
            for c in 0..=length {
                let is_false = c < length && !mask[c];
                if is_false {
                    if start_false.is_none() {
                        start_false = Some(c);
                    }
                } else {
                    if let Some(start) = start_false {
                        let island_len = c - start;
                        if island_len < min_island_length && island_len < length {
                            for i in start..c {
                                mask[i] = true;
                            }
                        }
                        start_false = None;
                    }
                }
            }
        }
        
        // 4. Apply mask and extract segments
        let mut seg_start: Option<usize> = None;
        for col in 0..length {
            if mask[col] {
                chars[col] = '-';
                if seg_start.is_none() {
                    seg_start = Some(col);
                }
            } else if let Some(start) = seg_start {
                masked_segments.push(MaskedSegment {
                    taxon: taxa.get(taxon_idx).cloned().unwrap_or_default(),
                    start,
                    end: col,
                });
                seg_start = None;
            }
        }

        if let Some(start) = seg_start {
            masked_segments.push(MaskedSegment {
                taxon: taxa.get(taxon_idx).cloned().unwrap_or_default(),
                start,
                end: length,
            });
        }

        cleaned_sequences.push(chars.into_iter().collect());
    }

    (cleaned_sequences, masked_segments)
}
