use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum AmbiguityStrategy {
    Keep,
    MajorityBase,
    ConvertToGap,
    FixedStandard, // R->A, Y->T, S->G, W->A, K->T, M->A, B->T, D->T, H->T, V->A
}

pub fn replace_character(sequences: &[String], find: char, replace: char) -> Vec<String> {
    let find_upper = find.to_ascii_uppercase();
    let find_lower = find.to_ascii_lowercase();

    sequences
        .iter()
        .map(|s| {
            s.chars()
                .map(|c| {
                    if c == find_upper || c == find_lower {
                        replace
                    } else {
                        c
                    }
                })
                .collect()
        })
        .collect()
}

/// Removes alignment columns where every sample contains only a gap or an
/// explicit missing-data state. Unlike the general column-gap filter, this is
/// safe for single-sample alignments and never removes a column containing an
/// observed nucleotide or ambiguity code.
pub fn remove_gap_only_columns(sequences: &[String]) -> (Vec<String>, Vec<usize>) {
    if sequences.is_empty() {
        return (Vec::new(), Vec::new());
    }

    let length = sequences.iter().map(String::len).min().unwrap_or(0);
    let mut kept_columns = Vec::with_capacity(length);
    let mut dropped_columns = Vec::new();

    for column in 0..length {
        let has_observed_state = sequences.iter().any(|sequence| {
            sequence.as_bytes().get(column).is_some_and(|base| {
                !matches!(base.to_ascii_uppercase(), b'-' | b'N' | b'?')
            })
        });
        if has_observed_state {
            kept_columns.push(column);
        } else {
            dropped_columns.push(column);
        }
    }

    if dropped_columns.is_empty() {
        return (sequences.to_vec(), dropped_columns);
    }

    let cleaned = sequences
        .iter()
        .map(|sequence| {
            kept_columns
                .iter()
                .filter_map(|column| sequence.as_bytes().get(*column).copied())
                .map(char::from)
                .collect()
        })
        .collect();
    (cleaned, dropped_columns)
}

pub fn convert_ambiguous_consensus(
    sequences: &[String],
    strategy: AmbiguityStrategy,
) -> Vec<String> {
    match strategy {
        AmbiguityStrategy::Keep => sequences.to_vec(),
        AmbiguityStrategy::ConvertToGap => {
            let ambiguities = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
                               'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v'];
            sequences
                .iter()
                .map(|s| {
                    s.chars()
                        .map(|c| if ambiguities.contains(&c) { '-' } else { c })
                        .collect()
                })
                .collect()
        }
        AmbiguityStrategy::FixedStandard => {
            // Matching PhyloProcessR convertAmbiguousConsensus defaults
            sequences
                .iter()
                .map(|s| {
                    s.chars()
                        .map(|c| match c {
                            'R' | 'r' => 'A',
                            'Y' | 'y' => 'T',
                            'S' | 's' => 'G',
                            'W' | 'w' => 'A',
                            'K' | 'k' => 'T',
                            'M' | 'm' => 'A',
                            'B' | 'b' => 'T',
                            'D' | 'd' => 'T',
                            'H' | 'h' => 'T',
                            'V' | 'v' => 'A',
                            other => other,
                        })
                        .collect()
                })
                .collect()
        }
        AmbiguityStrategy::MajorityBase => {
            if sequences.is_empty() {
                return sequences.to_vec();
            }
            let length = sequences[0].len();
            let _num_seqs = sequences.len();

            // Find majority unambiguous base per column
            let mut majority_bases = Vec::with_capacity(length);
            for col in 0..length {
                let mut counts = [0usize; 4]; // A, C, G, T
                for seq in sequences {
                    let b = seq.as_bytes().get(col).copied().unwrap_or(b'-').to_ascii_uppercase();
                    match b {
                        b'A' => counts[0] += 1,
                        b'C' => counts[1] += 1,
                        b'G' => counts[2] += 1,
                        b'T' => counts[3] += 1,
                        _ => {}
                    }
                }
                let max_idx = counts
                    .iter()
                    .enumerate()
                    .max_by_key(|&(_, count)| *count)
                    .map(|(idx, _)| idx)
                    .unwrap_or(0);
                let maj_char = match max_idx {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    _ => 'T',
                };
                majority_bases.push(maj_char);
            }

            sequences
                .iter()
                .map(|s| {
                    s.char_indices()
                        .map(|(col, c)| match c.to_ascii_uppercase() {
                            'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' => {
                                majority_bases.get(col).copied().unwrap_or('A')
                            }
                            other => other,
                        })
                        .collect()
                })
                .collect()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_replace_n() {
        let seqs = vec!["ATGN-NATG".to_string(), "nnnnNNNN".to_string()];
        let res = replace_character(&seqs, 'N', '-');
        assert_eq!(res[0], "ATG---ATG");
        assert_eq!(res[1], "--------");
    }

    #[test]
    fn test_convert_ambiguities_fixed() {
        let seqs = vec!["ATGRYSWKM".to_string()];
        let res = convert_ambiguous_consensus(&seqs, AmbiguityStrategy::FixedStandard);
        assert_eq!(res[0], "ATGATGATA");
    }

    #[test]
    fn test_remove_gap_and_missing_only_columns() {
        let seqs = vec!["A-N?C".to_string(), "T-??G".to_string()];
        let (cleaned, dropped) = remove_gap_only_columns(&seqs);
        assert_eq!(cleaned, vec!["AC", "TG"]);
        assert_eq!(dropped, vec![1, 2, 3]);
    }
}
