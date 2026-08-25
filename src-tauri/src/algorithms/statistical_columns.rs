use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum StatisticalColumnMethod {
    None,
    TrimalSimilarity,
    GblocksBlocks,
    Entropy,
}

impl Default for StatisticalColumnMethod {
    fn default() -> Self {
        Self::None
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum StatColGapTreatment {
    None,
    Half,
    All,
}

impl Default for StatColGapTreatment {
    fn default() -> Self {
        Self::Half
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum TrimalHeuristic {
    Custom,
    Gappyout,
    Strict,
    StrictPlus,
}

impl Default for TrimalHeuristic {
    fn default() -> Self {
        Self::Custom
    }
}

/// Computes pairwise similarity for a single column (0.0 to 1.0) among non-gap residues
pub fn compute_column_similarity_score(sequences: &[String], col: usize) -> f64 {
    let mut non_gap_chars = Vec::with_capacity(sequences.len());
    for s in sequences {
        if let Some(&b) = s.as_bytes().get(col) {
            let u = b.to_ascii_uppercase();
            if u != b'-' && u != b'?' && u != b'N' {
                non_gap_chars.push(u);
            }
        }
    }

    let n = non_gap_chars.len();
    if n <= 1 {
        return 0.0;
    }

    let total_pairs = (n * (n - 1)) / 2;
    let mut matches = 0usize;

    // Character frequency counting for fast O(N) pairwise similarity
    let mut counts = [0usize; 256];
    for &c in &non_gap_chars {
        counts[c as usize] += 1;
    }

    for &cnt in &counts {
        if cnt >= 2 {
            matches += (cnt * (cnt - 1)) / 2;
        }
    }

    matches as f64 / total_pairs as f64
}

/// Applies a sliding window smoothing of width W over a vector of raw column scores
pub fn apply_sliding_window(scores: &[f64], window_size: usize) -> Vec<f64> {
    let len = scores.len();
    if len == 0 || window_size <= 1 {
        return scores.to_vec();
    }

    let half = window_size / 2;
    let mut smoothed = Vec::with_capacity(len);

    for i in 0..len {
        let start = i.saturating_sub(half);
        let end = (i + half + 1).min(len);
        let window_slice = &scores[start..end];
        let sum: f64 = window_slice.iter().sum();
        smoothed.push(sum / window_slice.len() as f64);
    }

    smoothed
}

/// Computes column Shannon Entropy in bits (0.0 = completely identical, ~2.0 = uniform random)
pub fn compute_column_shannon_entropy(sequences: &[String], col: usize) -> f64 {
    let mut counts = [0usize; 256];
    let mut total = 0usize;

    for s in sequences {
        if let Some(&b) = s.as_bytes().get(col) {
            let u = b.to_ascii_uppercase();
            if u == b'A' || u == b'C' || u == b'G' || u == b'T' || u == b'U' {
                counts[u as usize] += 1;
                total += 1;
            }
        }
    }

    if total <= 1 {
        return 0.0;
    }

    let mut entropy = 0.0;
    for &cnt in &counts {
        if cnt > 0 {
            let p = cnt as f64 / total as f64;
            entropy -= p * p.log2();
        }
    }

    entropy
}

/// Computes column gap proportion (0.0 to 1.0)
pub fn compute_column_gap_fraction(sequences: &[String], col: usize) -> f64 {
    if sequences.is_empty() {
        return 1.0;
    }
    let mut gaps = 0usize;
    for s in sequences {
        if let Some(&b) = s.as_bytes().get(col) {
            let u = b.to_ascii_uppercase();
            if u == b'-' || u == b'?' || u == b'N' {
                gaps += 1;
            }
        }
    }
    gaps as f64 / sequences.len() as f64
}

/// Statistical Column Trimming Engine
/// Supports:
/// 1. trimAl Similarity & Sliding Window Consistency (Custom & Automated Heuristics)
/// 2. Gblocks Conserved Block Segmentation
/// 3. Shannon Information Entropy
pub fn trim_statistical_columns(
    sequences: &[String],
    method: StatisticalColumnMethod,
    similarity_threshold: f64,
    window_size: usize,
    heuristic: TrimalHeuristic,
    min_block_length: usize,
    _max_nonconserved: usize,
    gap_treatment: StatColGapTreatment,
    max_entropy_threshold: f64,
) -> (Vec<String>, Vec<usize>, HashMap<usize, String>) {
    if sequences.is_empty() || method == StatisticalColumnMethod::None {
        return (sequences.to_vec(), Vec::new(), HashMap::new());
    }

    let length = sequences[0].len();
    if length == 0 {
        return (sequences.to_vec(), Vec::new(), HashMap::new());
    }

    let mut dropped_cols = Vec::new();
    let mut reasons = HashMap::new();

    match method {
        StatisticalColumnMethod::TrimalSimilarity => {
            // Compute raw similarity scores
            let raw_scores: Vec<f64> = (0..length)
                .map(|col| compute_column_similarity_score(sequences, col))
                .collect();
            let smoothed_scores = apply_sliding_window(&raw_scores, window_size);

            let effective_threshold = match heuristic {
                TrimalHeuristic::Custom => similarity_threshold,
                TrimalHeuristic::Gappyout => {
                    // gappyout automatically focuses on columns above average alignment similarity
                    let avg: f64 = smoothed_scores.iter().sum::<f64>() / length.max(1) as f64;
                    (avg * 0.75).clamp(0.15, 0.60)
                }
                TrimalHeuristic::Strict => {
                    let avg: f64 = smoothed_scores.iter().sum::<f64>() / length.max(1) as f64;
                    (avg * 0.90).clamp(0.25, 0.70)
                }
                TrimalHeuristic::StrictPlus => {
                    let avg: f64 = smoothed_scores.iter().sum::<f64>() / length.max(1) as f64;
                    (avg * 1.10).clamp(0.35, 0.85)
                }
            };

            for col in 0..length {
                let score = smoothed_scores[col];
                let gap_frac = compute_column_gap_fraction(sequences, col);

                let is_dropped = if heuristic == TrimalHeuristic::Gappyout {
                    gap_frac > 0.60 || score < effective_threshold
                } else if heuristic == TrimalHeuristic::Strict || heuristic == TrimalHeuristic::StrictPlus {
                    gap_frac > 0.40 || score < effective_threshold
                } else {
                    score < effective_threshold
                };

                if is_dropped {
                    dropped_cols.push(col);
                    let label = match heuristic {
                        TrimalHeuristic::Custom => format!(
                            "trimAl Similarity ({:.2} < min {:.2}, win={})",
                            score, effective_threshold, window_size
                        ),
                        TrimalHeuristic::Gappyout => format!(
                            "trimAl Gappyout (Score {:.2}, Gap {:.0}%)",
                            score, gap_frac * 100.0
                        ),
                        TrimalHeuristic::Strict => format!(
                            "trimAl Strict (Score {:.2}, Gap {:.0}%)",
                            score, gap_frac * 100.0
                        ),
                        TrimalHeuristic::StrictPlus => format!(
                            "trimAl StrictPlus (Score {:.2}, Gap {:.0}%)",
                            score, gap_frac * 100.0
                        ),
                    };
                    reasons.insert(col, label);
                }
            }
        }

        StatisticalColumnMethod::GblocksBlocks => {
            // 1. Identify conserved status per column
            let raw_scores: Vec<f64> = (0..length)
                .map(|col| compute_column_similarity_score(sequences, col))
                .collect();
            let gap_fractions: Vec<f64> = (0..length)
                .map(|col| compute_column_gap_fraction(sequences, col))
                .collect();

            let mut is_conserved = vec![false; length];
            for col in 0..length {
                let sim = raw_scores[col];
                let gap = gap_fractions[col];

                let gap_ok = match gap_treatment {
                    StatColGapTreatment::None => gap == 0.0,
                    StatColGapTreatment::Half => gap <= 0.50,
                    StatColGapTreatment::All => true,
                };

                if sim >= 0.50 && gap_ok {
                    is_conserved[col] = true;
                }
            }

            // 2. Identify contiguous blocks of conserved positions (length >= min_block_length)
            let mut block_membership = vec![false; length];
            let mut col = 0;

            while col < length {
                if is_conserved[col] {
                    let block_start = col;
                    let mut block_end = col;
                    while block_end < length && is_conserved[block_end] {
                        block_end += 1;
                    }

                    let block_len = block_end - block_start;
                    if block_len >= min_block_length {
                        for b in block_start..block_end {
                            block_membership[b] = true;
                        }
                    }

                    col = block_end;
                } else {
                    col += 1;
                }
            }

            for c in 0..length {
                if !block_membership[c] {
                    dropped_cols.push(c);
                    let reason = if !is_conserved[c] {
                        format!("Gblocks Non-Conserved / Gap (sim {:.2})", raw_scores[c])
                    } else {
                        format!("Gblocks Fragment (< {} bp conserved block)", min_block_length)
                    };
                    reasons.insert(c, reason);
                }
            }
        }

        StatisticalColumnMethod::Entropy => {
            let raw_entropy: Vec<f64> = (0..length)
                .map(|col| compute_column_shannon_entropy(sequences, col))
                .collect();
            let smoothed_entropy = apply_sliding_window(&raw_entropy, window_size);

            for col in 0..length {
                let ent = smoothed_entropy[col];
                if ent > max_entropy_threshold {
                    dropped_cols.push(col);
                    reasons.insert(
                        col,
                        format!(
                            "High Entropy / Noise ({:.2} > max {:.2} bits, win={})",
                            ent, max_entropy_threshold, window_size
                        ),
                    );
                }
            }
        }

        StatisticalColumnMethod::None => {}
    }

    if dropped_cols.is_empty() {
        return (sequences.to_vec(), Vec::new(), HashMap::new());
    }

    let dropped_set: std::collections::HashSet<usize> = dropped_cols.iter().cloned().collect();
    let kept_cols: Vec<usize> = (0..length).filter(|c| !dropped_set.contains(c)).collect();

    let trimmed_seqs: Vec<String> = sequences
        .iter()
        .map(|s| {
            let bytes = s.as_bytes();
            kept_cols
                .iter()
                .filter_map(|&col| bytes.get(col).copied())
                .map(|b| b as char)
                .collect()
        })
        .collect();

    (trimmed_seqs, dropped_cols, reasons)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_column_similarity() {
        let seqs = vec![
            "AAAA".to_string(),
            "AAAA".to_string(),
            "AAAT".to_string(),
            "AATC".to_string(),
        ];
        // Col 0: 4 A's -> 1.0
        assert_eq!(compute_column_similarity_score(&seqs, 0), 1.0);
        // Col 3: A, A, T, C -> pairs: (A,A) match, (A,T) no, (A,C) no, (A,T) no, (A,C) no, (T,C) no -> 1/6
        assert!((compute_column_similarity_score(&seqs, 3) - (1.0 / 6.0)).abs() < 1e-5);
    }

    #[test]
    fn test_trimal_similarity_trimming() {
        let seqs = vec![
            "AAAAACCCCC".to_string(),
            "AAAAACCCCC".to_string(),
            "AAAAACCCCC".to_string(),
            "AAAATCCGCT".to_string(),
        ];
        let (trimmed, dropped, reasons) = trim_statistical_columns(
            &seqs,
            StatisticalColumnMethod::TrimalSimilarity,
            0.60,
            1,
            TrimalHeuristic::Custom,
            5,
            2,
            StatColGapTreatment::Half,
            1.5,
        );

        assert!(!dropped.is_empty());
        assert!(!reasons.is_empty());
        assert_eq!(trimmed[0].len(), seqs[0].len() - dropped.len());
    }

    #[test]
    fn test_gblocks_conserved_blocks() {
        let seqs = vec![
            "AAAAACGTACGAAAAA".to_string(),
            "AAAAATGCACTAAAAA".to_string(),
            "AAAAAGGTAACAAAAA".to_string(),
            "AAAAACCGTACAAAAA".to_string(),
        ];
        let (trimmed, dropped, _) = trim_statistical_columns(
            &seqs,
            StatisticalColumnMethod::GblocksBlocks,
            0.50,
            1,
            TrimalHeuristic::Custom,
            5,
            2,
            StatColGapTreatment::Half,
            1.5,
        );
        assert!(!dropped.is_empty());
        assert_eq!(trimmed[0], "AAAAAAAAAA");
    }

    #[test]
    fn test_shannon_entropy_trimming() {
        let seqs = vec![
            "AAAACGTN".to_string(),
            "AAAACGTN".to_string(),
            "AAAAACTG".to_string(),
            "AAAAGTCA".to_string(),
        ];
        let (trimmed, dropped, _) = trim_statistical_columns(
            &seqs,
            StatisticalColumnMethod::Entropy,
            0.50,
            1,
            TrimalHeuristic::Custom,
            5,
            2,
            StatColGapTreatment::Half,
            1.0,
        );
        // Col 0-3 (pure A) has entropy 0.0 <= 1.0 (kept)
        // High variable columns have entropy > 1.0 (dropped)
        assert!(dropped.contains(&4) || dropped.contains(&5) || dropped.contains(&6));
        assert!(!dropped.contains(&0));
    }
}
