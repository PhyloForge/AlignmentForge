use std::collections::HashMap;

use crate::algorithms::orf::reverse_complement_dna;
use crate::models::Alignment;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReferenceMatch {
    pub start: usize,
    pub end: usize,
    pub identity_percent: f64,
    pub coverage_percent: f64,
    pub is_reverse: bool,
}

#[derive(Clone, Copy, Default)]
struct LocalState {
    score: i32,
    start_target: usize,
    matches: usize,
    reference_bases: usize,
    target_bases: usize,
}

fn resolved_base(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(b'A'),
        b'C' => Some(b'C'),
        b'G' => Some(b'G'),
        b'T' | b'U' => Some(b'T'),
        _ => None,
    }
}

fn normalize_reference(reference: &str) -> Vec<u8> {
    reference
        .as_bytes()
        .iter()
        .filter_map(|base| resolved_base(*base))
        .collect()
}

fn consensus_with_column_map(sequences: &[String]) -> (Vec<u8>, Vec<usize>) {
    let length = sequences.iter().map(String::len).min().unwrap_or(0);
    let mut consensus = Vec::with_capacity(length);
    let mut column_map = Vec::with_capacity(length);

    for column in 0..length {
        let mut counts = [0usize; 4];
        for sequence in sequences {
            let Some(base) = sequence.as_bytes().get(column).and_then(|base| resolved_base(*base)) else {
                continue;
            };
            let index = match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                _ => 3,
            };
            counts[index] += 1;
        }
        let Some((index, count)) = counts.iter().enumerate().max_by_key(|(_, count)| *count) else {
            continue;
        };
        if *count == 0 {
            continue;
        }
        consensus.push([b'A', b'C', b'G', b'T'][index]);
        column_map.push(column);
    }

    (consensus, column_map)
}

fn encode_kmer(kmer: &[u8]) -> Option<u64> {
    let mut encoded = 0u64;
    for base in kmer {
        encoded <<= 2;
        encoded |= match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
    }
    Some(encoded)
}

fn dominant_anchor_diagonal(reference: &[u8], target: &[u8]) -> Option<isize> {
    for kmer_length in [13usize, 9, 7, 5] {
        if reference.len() < kmer_length || target.len() < kmer_length {
            continue;
        }
        let mut reference_kmers: HashMap<u64, Vec<usize>> = HashMap::new();
        for position in 0..=reference.len() - kmer_length {
            if let Some(kmer) = encode_kmer(&reference[position..position + kmer_length]) {
                reference_kmers.entry(kmer).or_default().push(position);
            }
        }
        let mut diagonals: HashMap<isize, usize> = HashMap::new();
        for target_position in 0..=target.len() - kmer_length {
            let Some(kmer) = encode_kmer(&target[target_position..target_position + kmer_length]) else {
                continue;
            };
            if let Some(reference_positions) = reference_kmers.get(&kmer) {
                for reference_position in reference_positions.iter().take(16) {
                    *diagonals
                        .entry(target_position as isize - *reference_position as isize)
                        .or_insert(0) += 1;
                }
            }
        }
        if let Some((diagonal, support)) = diagonals.into_iter().max_by_key(|(_, support)| *support) {
            if support >= 2 || kmer_length <= 7 {
                return Some(diagonal);
            }
        }
    }
    None
}

fn extend_state(mut state: LocalState, score_delta: i32, reference: bool, target: bool, matched: bool, start_target: usize) -> LocalState {
    if state.score <= 0 {
        state = LocalState {
            score: 0,
            start_target,
            ..LocalState::default()
        };
    }
    state.score += score_delta;
    state.reference_bases += usize::from(reference);
    state.target_bases += usize::from(target);
    state.matches += usize::from(matched);
    if state.score <= 0 {
        LocalState::default()
    } else {
        state
    }
}

fn banded_local_match(reference: &[u8], target: &[u8], diagonal: isize) -> Option<(usize, usize, f64, f64, i32)> {
    if reference.is_empty() || target.is_empty() {
        return None;
    }
    let band = 64isize;
    let target_bound = isize::try_from(target.len()).ok()?;
    let mut previous = vec![LocalState::default(); target.len() + 1];
    let mut current = vec![LocalState::default(); target.len() + 1];
    let mut previous_start = 1usize;
    let mut previous_end = 0usize;
    let mut best = LocalState::default();
    let mut best_end = 0usize;

    for reference_index in 1..=reference.len() {
        let reference_position = isize::try_from(reference_index).ok()?;
        let center = reference_position.saturating_add(diagonal);
        let signed_start = center.saturating_sub(band).max(1);
        let signed_end = center.saturating_add(band).min(target_bound);
        if signed_end < 1 || signed_start > signed_end {
            continue;
        }
        let start = usize::try_from(signed_start).ok()?;
        let end = usize::try_from(signed_end).ok()?;
        current[start.saturating_sub(1)..=end.saturating_add(1).min(target.len())]
            .fill(LocalState::default());

        for target_index in start..=end {
            let reference_base = reference[reference_index - 1];
            let target_base = target[target_index - 1];
            let diagonal_state = if target_index - 1 >= previous_start
                && target_index - 1 <= previous_end
            {
                previous[target_index - 1]
            } else {
                LocalState::default()
            };
            let up_state = if target_index >= previous_start && target_index <= previous_end {
                previous[target_index]
            } else {
                LocalState::default()
            };
            let left_state = current[target_index - 1];
            let mut selected = extend_state(
                diagonal_state,
                if reference_base == target_base { 3 } else { -2 },
                true,
                true,
                reference_base == target_base,
                target_index - 1,
            );
            let up = extend_state(up_state, -3, true, false, false, target_index - 1);
            if up.score > selected.score {
                selected = up;
            }
            let left = extend_state(left_state, -3, false, true, false, target_index - 1);
            if left.score > selected.score {
                selected = left;
            }
            current[target_index] = selected;
            if selected.score > best.score {
                best = selected;
                best_end = target_index;
            }
        }
        std::mem::swap(&mut previous, &mut current);
        previous_start = start;
        previous_end = end;
    }

    if best.score <= 0 || best_end <= best.start_target || best.reference_bases == 0 {
        return None;
    }
    let identity = best.matches as f64
        / best.reference_bases.max(best.target_bases).max(1) as f64
        * 100.0;
    let coverage = best.reference_bases as f64 / reference.len() as f64 * 100.0;
    Some((best.start_target, best_end, identity, coverage, best.score))
}

pub fn match_reference_to_alignment(sequences: &[String], reference: &str) -> Option<ReferenceMatch> {
    let (target, column_map) = consensus_with_column_map(sequences);
    let forward = normalize_reference(reference);
    if forward.len() < 12 || target.len() < 12 {
        return None;
    }
    let reverse = normalize_reference(&reverse_complement_dna(reference));
    let mut best: Option<(usize, usize, f64, f64, i32, bool)> = None;

    for (candidate, is_reverse) in [(&forward, false), (&reverse, true)] {
        let Some(diagonal) = dominant_anchor_diagonal(candidate, &target) else {
            continue;
        };
        let Some((start, end, identity, coverage, score)) =
            banded_local_match(candidate, &target, diagonal)
        else {
            continue;
        };
        if identity < 70.0 || coverage < 70.0 {
            continue;
        }
        if best.as_ref().is_none_or(|current| score > current.4) {
            best = Some((start, end, identity, coverage, score, is_reverse));
        }
    }

    let (target_start, target_end, identity, coverage, _, is_reverse) = best?;
    let start = *column_map.get(target_start)?;
    let end = column_map.get(target_end.saturating_sub(1))?.saturating_add(1);
    (end > start).then_some(ReferenceMatch {
        start,
        end,
        identity_percent: identity,
        coverage_percent: coverage,
        is_reverse,
    })
}

pub fn intron_alignment_from_reference(
    alignment: &Alignment,
    reference: &str,
) -> Option<(Alignment, ReferenceMatch)> {
    let reference_match = match_reference_to_alignment(&alignment.sequences, reference)?;
    let sequences = alignment
        .sequences
        .iter()
        .map(|sequence| {
            let prefix = sequence.get(..reference_match.start).unwrap_or_default();
            let suffix = sequence.get(reference_match.end..).unwrap_or_default();
            format!("{prefix}{suffix}")
        })
        .collect::<Vec<_>>();
    let intron = Alignment::new(
        format!("{}_intron", alignment.id),
        format!("{}_intron", alignment.file_name),
        alignment.file_path.clone(),
        alignment.format,
        alignment.taxa.clone(),
        sequences,
    );
    Some((intron, reference_match))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::models::AlignmentFormat;

    #[test]
    fn reference_anchors_exon_and_preserves_flanking_introns() {
        let exon = "ATGGCTGCTGCTGCTGCTGCTGCTGCTTAA";
        let alignment = Alignment::new(
            "locus_1".to_string(),
            "locus_1.fa".to_string(),
            "/tmp/locus_1.fa".to_string(),
            AlignmentFormat::Fasta,
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            vec![
                format!("CCCCCC{exon}GGGGGG"),
                format!("CCCCCC{exon}GGGGGG"),
                format!("CCCCCC{exon}GGGGGG"),
            ],
        );

        let (intron, matched) = intron_alignment_from_reference(&alignment, exon).unwrap();
        assert_eq!(matched.start, 6);
        assert_eq!(matched.end, 6 + exon.len());
        assert!(matched.identity_percent > 99.0);
        assert_eq!(intron.sequences[0], "CCCCCCGGGGGG");
    }

    #[test]
    fn banded_match_ignores_bands_completely_outside_the_target() {
        let reference = vec![b'A'; 256];
        let target = vec![b'A'; 128];

        let _ = banded_local_match(&reference, &target, -200);
        assert!(banded_local_match(&reference, &target, isize::MIN).is_none());
        assert!(banded_local_match(&reference, &target, isize::MAX).is_none());
    }

    #[test]
    fn reference_with_long_unmatched_prefix_does_not_overrun_target() {
        let exon = "ATGGCTGATCGTACCGGTTACGATGCTGACTAA";
        let reference = format!("{}{exon}", "C".repeat(256));
        let sequences = vec![exon.to_string(); 3];

        // Coverage is intentionally too low for a match, but the negative
        // anchor diagonal must be handled without panicking.
        assert!(match_reference_to_alignment(&sequences, &reference).is_none());
    }

    #[test]
    fn long_shifted_band_stays_within_allocated_columns() {
        let reference = vec![b'A'; 6_000];
        let target = vec![b'A'; 5_000];
        let result = banded_local_match(&reference, &target, -1_000);

        assert!(result.is_some());
    }
}
