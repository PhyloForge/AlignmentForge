/// Counts variable and parsimony-informative sites in one pass.
///
/// A variable site contains at least two distinct resolved nucleotide states. A
/// parsimony-informative site contains at least two resolved states that each
/// occur in at least two sequences. Gaps, `?`, and `N`/`n` are always missing;
/// IUPAC ambiguity codes are also missing when `exclude_ambiguities` is true.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteStatistics {
    pub variable_count: usize,
    pub variable_percent: f64,
    pub variable_mask: Vec<bool>,
    pub pis_count: usize,
    pub pis_percent: f64,
    pub pis_mask: Vec<bool>,
}

pub fn calculate_site_statistics(
    sequences: &[String],
    exclude_ambiguities: bool,
) -> SiteStatistics {
    let length = sequences.first().map_or(0, |sequence| sequence.len());
    let mut variable_mask = vec![false; length];
    let mut pis_mask = vec![false; length];
    let mut variable_count = 0usize;
    let mut pis_count = 0usize;

    // Stack array for character counts (zero heap allocation)
    let mut char_counts = [0usize; 128];

    for col in 0..length {
        char_counts.fill(0);
        let mut distinct_states = 0usize;
        let mut states_with_min_two = 0usize;

        for seq in sequences {
            if let Some(&b) = seq.as_bytes().get(col) {
                if let Some(state) = resolved_state(b, exclude_ambiguities) {
                    let state_index = state as usize;
                    if char_counts[state_index] == 0 {
                        distinct_states += 1;
                    }
                    char_counts[state_index] += 1;
                    if char_counts[state_index] == 2 {
                        states_with_min_two += 1;
                    }
                }
            }
        }

        if distinct_states >= 2 {
            variable_mask[col] = true;
            variable_count += 1;
        }
        if states_with_min_two >= 2 {
            pis_mask[col] = true;
            pis_count += 1;
        }
    }

    let variable_percent = if length > 0 {
        (variable_count as f64 / length as f64) * 100.0
    } else {
        0.0
    };
    let pis_percent = if length > 0 {
        (pis_count as f64 / length as f64) * 100.0
    } else {
        0.0
    };

    SiteStatistics {
        variable_count,
        variable_percent,
        variable_mask,
        pis_count,
        pis_percent,
        pis_mask,
    }
}

/// Port of PhyloProcessR `informativeSites.R`.
pub fn calculate_parsimony_informative_sites(
    sequences: &[String],
    exclude_ambiguities: bool,
) -> (usize, f64, Vec<bool>) {
    let stats = calculate_site_statistics(sequences, exclude_ambiguities);
    (stats.pis_count, stats.pis_percent, stats.pis_mask)
}

pub fn calculate_variable_sites(
    sequences: &[String],
    exclude_ambiguities: bool,
) -> (usize, f64, Vec<bool>) {
    let stats = calculate_site_statistics(sequences, exclude_ambiguities);
    (
        stats.variable_count,
        stats.variable_percent,
        stats.variable_mask,
    )
}

#[inline(always)]
fn resolved_state(character: u8, exclude_ambiguities: bool) -> Option<u8> {
    let upper = character.to_ascii_uppercase();
    match upper {
        b'A' | b'C' | b'G' | b'T' => Some(upper),
        b'U' => Some(b'T'),
        b'-' | b'?' | b'N' => None,
        ambiguity if exclude_ambiguities && is_iupac_ambiguity(ambiguity) => None,
        other if other < 128 => Some(other),
        _ => None,
    }
}

#[inline(always)]
fn is_iupac_ambiguity(c: u8) -> bool {
    matches!(
        c,
        b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V'
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parsimony_informative_sites() {
        let seqs = vec![
            "AAAA-".to_string(), // Col 0: AAAA (invariant -> not PIS)
            "AATT-".to_string(), // Col 1: AATT (A:2, T:2 -> PIS)
            "TTTT-".to_string(), // Col 2: TTTT (invariant -> not PIS)
            "TTCC-".to_string(), // Col 3: TTCC (T:2, C:2 -> PIS)
        ];

        let (count, pct, mask) = calculate_parsimony_informative_sites(&seqs, true);
        assert_eq!(count, 2);
        assert_eq!(mask, vec![true, true, false, false, false]);
        assert_eq!(pct, 40.0);
    }

    #[test]
    fn test_variable_and_informative_sites_ignore_missing_states() {
        let seqs = vec![
            "AAA-AA".to_string(),
            "AAANAA".to_string(),
            "AAA?NC".to_string(),
            "AACn?C".to_string(),
            "AACRAN".to_string(),
            "ATCY-?".to_string(),
        ];

        let stats = calculate_site_statistics(&seqs, true);
        assert_eq!(stats.variable_count, 3);
        assert_eq!(
            stats.variable_mask,
            vec![false, true, true, false, false, true]
        );
        assert_eq!(stats.pis_count, 2);
        assert_eq!(stats.pis_mask, vec![false, false, true, false, false, true]);
        assert_eq!(stats.variable_percent, 50.0);
        assert!((stats.pis_percent - (2.0 / 6.0 * 100.0)).abs() < 1e-10);
    }
}
