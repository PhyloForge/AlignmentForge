use serde::{Deserialize, Serialize};
use crate::models::MaskedSegment;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum GeneticCode {
    Standard,
    VertebrateMitochondrial,
    InvertebrateMitochondrial,
}

impl Default for GeneticCode {
    fn default() -> Self {
        GeneticCode::Standard
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum StopCodonAction {
    RemoveSample,
    MaskCodon,
    Keep,
}

impl Default for StopCodonAction {
    fn default() -> Self {
        StopCodonAction::RemoveSample
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum OrfSearchMode {
    ContinuousCds,
    BestSharedSegment,
    ReferenceGuided,
    ReferenceCandidateOrf,
}

impl Default for OrfSearchMode {
    fn default() -> Self {
        OrfSearchMode::ContinuousCds
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrfConfig {
    pub enable_orf: bool,
    pub auto_shift_frame: bool,
    pub auto_flip_reverse: bool,
    pub stop_codon_action: StopCodonAction,
    pub genetic_code: GeneticCode,
    pub search_mode: OrfSearchMode,
    pub min_shared_support_percent: f64,
    pub min_segment_aa: usize,
    pub min_coding_score: f64,
    pub exclude_uce: bool,
    pub fail_if_no_orf: bool,
    pub max_stop_codons_sample: usize,
    pub max_stop_codons_locus: usize,
    pub macse_trim_terminal: bool,
    pub macse_max_internal_sample: usize,
    pub macse_max_internal_locus: usize,
}

impl Default for OrfConfig {
    fn default() -> Self {
        Self {
            enable_orf: false,
            auto_shift_frame: true,
            auto_flip_reverse: true,
            stop_codon_action: StopCodonAction::RemoveSample,
            genetic_code: GeneticCode::Standard,
            search_mode: OrfSearchMode::ContinuousCds,
            min_shared_support_percent: 90.0,
            min_segment_aa: 35,
            min_coding_score: 40.0,
            exclude_uce: true,
            fail_if_no_orf: true, max_stop_codons_sample: 2, max_stop_codons_locus: 5, macse_trim_terminal: true, macse_max_internal_sample: 3, macse_max_internal_locus: 10,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReadingFrameResult {
    pub frame: i8, // +1, +2, +3, -1, -2, -3
    pub is_reverse: bool,
    pub offset: usize, // 0, 1, or 2 leading bp to trim
    pub clean_taxa_count: usize,
    pub stop_taxa_count: usize,
    pub score: i64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SharedOrfSegmentResult {
    pub frame: i8,
    pub is_reverse: bool,
    pub start: usize,
    pub end: usize,
    pub support_count: usize,
    pub length_codons: usize,
    pub informative_codons: usize,
    pub coding_score: f64,
    pub amino_acid_conservation: f64,
    pub frame_contrast: f64,
}

#[derive(Debug, Clone)]
pub struct OrfOptimizationResult {
    pub taxa: Vec<String>,
    pub sequences: Vec<String>,
    pub dropped_taxa: Vec<String>,
    pub trimmed_columns: Vec<usize>,
    pub masked_segments: Vec<MaskedSegment>,
    pub found_valid_orf: bool,
    pub is_reverse: bool,
    pub orf_evaluated: bool,
    pub candidate_found: bool,
    pub candidate_frame: Option<i8>,
    pub candidate_start: Option<usize>,
    pub candidate_end: Option<usize>,
    pub candidate_support_count: usize,
    pub candidate_support_percent: f64,
    pub candidate_length_aa: usize,
    pub coding_score: f64,
    pub amino_acid_conservation: f64,
    pub frame_contrast: f64,
    pub retained_samples: usize,
}

/// Checks whether a locus identifier indicates sequence that should not be
/// interpreted as exon-only CDS. This safeguard is intentionally conservative:
/// without explicit annotations, AlignmentForge can only use identifier hints.
pub fn should_skip_orf_locus(locus_id: &str, search_mode: OrfSearchMode) -> bool {
    let lower = locus_id.to_ascii_lowercase();
    let always_skip = lower.starts_with("uce-")
        || lower.starts_with("uce_")
        || lower.starts_with("uce")
        || lower.contains("noncoding")
        || lower.contains("non-coding")
        || lower.contains("intergenic")
        || lower.contains("intron");
    let continuous_cds_only = lower.contains("supercontig") || lower.contains("flanking");
    always_skip
        || (search_mode == OrfSearchMode::ContinuousCds && continuous_cds_only)
}

/// Computes reverse complement of a DNA sequence
pub fn reverse_complement_dna(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'a' => 't',
            'T' | 'U' => 'A',
            't' | 'u' => 'a',
            'C' => 'G',
            'c' => 'g',
            'G' => 'C',
            'g' => 'c',
            'R' => 'Y',
            'r' => 'y',
            'Y' => 'R',
            'y' => 'r',
            'S' => 'S',
            's' => 's',
            'W' => 'W',
            'w' => 'w',
            'K' => 'M',
            'k' => 'm',
            'M' => 'K',
            'm' => 'k',
            'B' => 'V',
            'b' => 'v',
            'D' => 'H',
            'd' => 'h',
            'H' => 'D',
            'h' => 'd',
            'V' => 'B',
            'v' => 'b',
            '-' => '-',
            '?' => '?',
            'N' => 'N',
            'n' => 'n',
            other => other,
        })
        .collect()
}

/// Translates a single 3-base codon into an Amino Acid or Stop ('*')
pub fn translate_codon(codon: &[u8], code: GeneticCode) -> char {
    if codon.len() < 3 {
        return '-';
    }
    let c1 = (codon[0] as char).to_ascii_uppercase();
    let c2 = (codon[1] as char).to_ascii_uppercase();
    let c3 = (codon[2] as char).to_ascii_uppercase();

    if c1 == '-' || c2 == '-' || c3 == '-' {
        return '-';
    }
    if c1 == '?' || c2 == '?' || c3 == '?' || c1 == 'N' || c2 == 'N' || c3 == 'N' {
        return 'X';
    }

    let trip = [c1 as u8, c2 as u8, c3 as u8];

    match code {
        GeneticCode::Standard => match &trip {
            b"TAA" | b"TAG" | b"TGA" => '*',
            b"TTT" | b"TTC" => 'F',
            b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
            b"ATT" | b"ATC" | b"ATA" => 'I',
            b"ATG" => 'M',
            b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
            b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
            b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
            b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
            b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
            b"TAT" | b"TAC" => 'Y',
            b"CAT" | b"CAC" => 'H',
            b"CAA" | b"CAG" => 'Q',
            b"AAT" | b"AAC" => 'N',
            b"AAA" | b"AAG" => 'K',
            b"GAT" | b"GAC" => 'D',
            b"GAA" | b"GAG" => 'E',
            b"TGT" | b"TGC" => 'C',
            b"TGG" => 'W',
            b"CGT" | b"CGC" | b"CGA" | b"CGG" | b"AGA" | b"AGG" => 'R',
            b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
            _ => 'X',
        },
        GeneticCode::VertebrateMitochondrial => match &trip {
            b"TAA" | b"TAG" | b"AGA" | b"AGG" => '*',
            b"TGA" => 'W',
            b"ATA" => 'M',
            b"TTT" | b"TTC" => 'F',
            b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
            b"ATT" | b"ATC" => 'I',
            b"ATG" => 'M',
            b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
            b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
            b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
            b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
            b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
            b"TAT" | b"TAC" => 'Y',
            b"CAT" | b"CAC" => 'H',
            b"CAA" | b"CAG" => 'Q',
            b"AAT" | b"AAC" => 'N',
            b"AAA" | b"AAG" => 'K',
            b"GAT" | b"GAC" => 'D',
            b"GAA" | b"GAG" => 'E',
            b"TGT" | b"TGC" => 'C',
            b"TGG" => 'W',
            b"CGT" | b"CGC" | b"CGA" | b"CGG" => 'R',
            b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
            _ => 'X',
        },
        GeneticCode::InvertebrateMitochondrial => match &trip {
            b"TAA" | b"TAG" => '*',
            b"TGA" => 'W',
            b"ATA" => 'M',
            b"AGA" | b"AGG" => 'S',
            b"TTT" | b"TTC" => 'F',
            b"TTA" | b"TTG" | b"CTT" | b"CTC" | b"CTA" | b"CTG" => 'L',
            b"ATT" | b"ATC" => 'I',
            b"ATG" => 'M',
            b"GTT" | b"GTC" | b"GTA" | b"GTG" => 'V',
            b"TCT" | b"TCC" | b"TCA" | b"TCG" | b"AGT" | b"AGC" => 'S',
            b"CCT" | b"CCC" | b"CCA" | b"CCG" => 'P',
            b"ACT" | b"ACC" | b"ACA" | b"ACG" => 'T',
            b"GCT" | b"GCC" | b"GCA" | b"GCG" => 'A',
            b"TAT" | b"TAC" => 'Y',
            b"CAT" | b"CAC" => 'H',
            b"CAA" | b"CAG" => 'Q',
            b"AAT" | b"AAC" => 'N',
            b"AAA" | b"AAG" => 'K',
            b"GAT" | b"GAC" => 'D',
            b"GAA" | b"GAG" => 'E',
            b"TGT" | b"TGC" => 'C',
            b"TGG" => 'W',
            b"CGT" | b"CGC" | b"CGA" | b"CGG" => 'R',
            b"GGT" | b"GGC" | b"GGA" | b"GGG" => 'G',
            _ => 'X',
        },
    }
}

/// Evaluates one explicit reading frame across the alignment. Positive frames
/// use the alignment as supplied; negative frames use its reverse complement.
pub fn evaluate_reading_frame(
    sequences: &[String],
    code: GeneticCode,
    frame: i8,
) -> ReadingFrameResult {
    let normalized_frame = if (-3..=-1).contains(&frame) || (1..=3).contains(&frame) {
        frame
    } else {
        1
    };
    let is_reverse = normalized_frame < 0;
    let offset = normalized_frame.unsigned_abs() as usize - 1;
    let reverse_sequences: Vec<String>;
    let oriented = if is_reverse {
        reverse_sequences = sequences
            .iter()
            .map(|sequence| reverse_complement_dna(sequence))
            .collect();
        &reverse_sequences
    } else {
        sequences
    };
    let mut clean_taxa_count = 0usize;
    let mut stop_taxa_count = 0usize;
    let mut total_clean_codons = 0usize;
    let mut start_codons = 0usize;
    let mut terminal_stops = 0usize;

    for sequence in oriented {
        if sequence.len() < offset + 3 {
            continue;
        }
        let codon_count = (sequence.len() - offset) / 3;
        let mut has_internal_stop = false;
        for codon_index in 0..codon_count {
            let start = offset + codon_index * 3;
            let codon = &sequence.as_bytes()[start..start + 3];
            let amino_acid = translate_codon(codon, code);
            if codon_index == 0 && codon.eq_ignore_ascii_case(b"ATG") {
                start_codons += 1;
            }
            if amino_acid == '*' && codon_index + 1 == codon_count {
                terminal_stops += 1;
            } else if amino_acid == '*' {
                has_internal_stop = true;
            }
        }
        if has_internal_stop {
            stop_taxa_count += 1;
        } else {
            clean_taxa_count += 1;
            total_clean_codons += codon_count;
        }
    }

    ReadingFrameResult {
        frame: normalized_frame,
        is_reverse,
        offset,
        clean_taxa_count,
        stop_taxa_count,
        score: (total_clean_codons as i64 * 100)
            + (start_codons as i64 * 250)
            + (terminal_stops as i64 * 150)
            - (stop_taxa_count as i64 * 80),
    }
}

/// Evaluates all 6 reading frames across the alignment and returns the best consensus frame.
pub fn find_optimal_reading_frame(sequences: &[String], code: GeneticCode) -> ReadingFrameResult {
    let mut best_result = evaluate_reading_frame(sequences, code, 1);
    for frame in [2, 3, -1, -2, -3] {
        let candidate = evaluate_reading_frame(sequences, code, frame);
        if candidate.score > best_result.score {
            best_result = candidate;
        }
    }
    best_result
}

fn protein_profile_conservation(
    sequences: &[String],
    offset: usize,
    code: GeneticCode,
) -> f64 {
    let length = sequences.iter().map(String::len).min().unwrap_or(0);
    if length < offset + 3 {
        return 0.0;
    }
    let codon_count = (length - offset) / 3;
    let mut conservation_sum = 0.0;
    let mut informative_columns = 0usize;

    for codon_index in 0..codon_count {
        let start = offset + codon_index * 3;
        let mut counts: HashMap<char, usize> = HashMap::new();
        for sequence in sequences {
            let amino_acid = translate_codon(&sequence.as_bytes()[start..start + 3], code);
            if !matches!(amino_acid, '*' | 'X' | '-') {
                *counts.entry(amino_acid).or_insert(0) += 1;
            }
        }
        let total: usize = counts.values().sum();
        if total >= 2 {
            conservation_sum += counts.values().copied().max().unwrap_or(0) as f64 / total as f64;
            informative_columns += 1;
        }
    }

    if informative_columns > 0 {
        conservation_sum / informative_columns as f64
    } else {
        0.0
    }
}

fn synonymous_change_fraction(sequences: &[String], code: GeneticCode) -> f64 {
    let length = sequences.iter().map(String::len).min().unwrap_or(0);
    let codon_count = length / 3;
    let mut substitutions = 0usize;
    let mut synonymous = 0usize;

    for codon_index in 0..codon_count {
        let start = codon_index * 3;
        let mut counts: HashMap<[u8; 3], usize> = HashMap::new();
        for sequence in sequences {
            let bytes = sequence.as_bytes();
            let codon = [
                bytes[start].to_ascii_uppercase(),
                bytes[start + 1].to_ascii_uppercase(),
                bytes[start + 2].to_ascii_uppercase(),
            ];
            if !matches!(translate_codon(&codon, code), '*' | 'X' | '-') {
                *counts.entry(codon).or_insert(0) += 1;
            }
        }
        let Some((&consensus, _)) = counts.iter().max_by_key(|(_, count)| *count) else {
            continue;
        };
        let consensus_amino_acid = translate_codon(&consensus, code);
        for (codon, count) in counts {
            if codon != consensus {
                substitutions += count;
                if translate_codon(&codon, code) == consensus_amino_acid {
                    synonymous += count;
                }
            }
        }
    }

    if substitutions >= 5 {
        synonymous as f64 / substitutions as f64
    } else {
        0.0
    }
}

fn codon_position_variability(sequences: &[String]) -> [f64; 3] {
    let length = sequences.iter().map(String::len).min().unwrap_or(0);
    let codon_count = length / 3;
    let mut sums = [0.0; 3];
    let mut columns = [0usize; 3];

    for codon_index in 0..codon_count {
        for position in 0..3 {
            let column = codon_index * 3 + position;
            let mut counts = [0usize; 4];
            for sequence in sequences {
                match sequence.as_bytes()[column].to_ascii_uppercase() {
                    b'A' => counts[0] += 1,
                    b'C' => counts[1] += 1,
                    b'G' => counts[2] += 1,
                    b'T' | b'U' => counts[3] += 1,
                    _ => {}
                }
            }
            let total: usize = counts.iter().sum();
            if total >= 2 {
                sums[position] +=
                    1.0 - counts.iter().copied().max().unwrap_or(0) as f64 / total as f64;
                columns[position] += 1;
            }
        }
    }

    for position in 0..3 {
        if columns[position] > 0 {
            sums[position] /= columns[position] as f64;
        }
    }
    sums
}

fn clamp_unit(value: f64) -> f64 {
    value.clamp(0.0, 1.0)
}

/// Scores evidence that a stop-free interval behaves like a shared protein
/// profile rather than a chance ORF in homologous non-coding DNA.
fn candidate_coding_evidence(
    oriented_sequences: &[String],
    start: usize,
    end: usize,
    code: GeneticCode,
) -> (f64, f64, f64) {
    if end <= start || oriented_sequences.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let segments: Vec<String> = oriented_sequences
        .iter()
        .filter_map(|sequence| sequence.get(start..end).map(str::to_string))
        .collect();
    if segments.len() < 2 {
        return (0.0, 0.0, 0.0);
    }

    let selected_conservation = protein_profile_conservation(&segments, 0, code);
    let reverse_segments: Vec<String> = segments
        .iter()
        .map(|segment| reverse_complement_dna(segment))
        .collect();
    let alternative_conservation = [
        protein_profile_conservation(&segments, 1, code),
        protein_profile_conservation(&segments, 2, code),
        protein_profile_conservation(&reverse_segments, 0, code),
        protein_profile_conservation(&reverse_segments, 1, code),
        protein_profile_conservation(&reverse_segments, 2, code),
    ]
    .into_iter()
    .fold(0.0_f64, f64::max);
    let frame_contrast = (selected_conservation - alternative_conservation).max(0.0);
    let synonymous_fraction = synonymous_change_fraction(&segments, code);
    let variability = codon_position_variability(&segments);
    let periodicity =
        (variability[2] - (variability[0] + variability[1]) / 2.0).max(0.0);

    let start_fraction = segments
        .iter()
        .filter(|segment| segment.as_bytes().get(0..3) == Some(b"ATG"))
        .count() as f64
        / segments.len() as f64;
    let boundary_stop_fraction = oriented_sequences
        .iter()
        .filter(|sequence| {
            sequence
                .as_bytes()
                .get(end..end.saturating_add(3))
                .is_some_and(|codon| translate_codon(codon, code) == '*')
        })
        .count() as f64
        / oriented_sequences.len() as f64;

    let score = clamp_unit((selected_conservation - 0.75) / 0.20) * 25.0
        + clamp_unit(frame_contrast / 0.08) * 30.0
        + clamp_unit((synonymous_fraction - 0.25) / 0.50) * 25.0
        + clamp_unit(periodicity / 0.10) * 10.0
        + start_fraction * 5.0
        + boundary_stop_fraction * 5.0;

    (
        score.clamp(0.0, 100.0),
        selected_conservation,
        frame_contrast,
    )
}

/// Finds one continuous, aligned stop-free segment supported by as many samples
/// as possible. Candidate intervals are the maximal stop-free runs observed in
/// individual samples. They are scored by sample support × length, with resolved
/// codon content used as a tie-breaker so missing-only regions are not preferred.
pub fn find_best_shared_orf_segment(
    sequences: &[String],
    code: GeneticCode,
    min_shared_support_percent: f64,
    min_segment_aa: usize,
    min_coding_score: f64,
) -> Option<SharedOrfSegmentResult> {
    if sequences.is_empty() {
        return None;
    }

    let min_support = (((sequences.len() as f64)
        * min_shared_support_percent.clamp(0.0, 100.0)
        / 100.0)
        .ceil() as usize)
        .max(1)
        .min(sequences.len());
    let minimum_codons = min_segment_aa.max(1);
    let reverse_sequences: Vec<String> = sequences
        .iter()
        .map(|sequence| reverse_complement_dna(sequence))
        .collect();
    let mut candidates: Vec<SharedOrfSegmentResult> = Vec::new();

    for is_reverse in [false, true] {
        let oriented = if is_reverse {
            &reverse_sequences
        } else {
            sequences
        };
        let alignment_length = oriented.iter().map(String::len).min().unwrap_or(0);

        for offset in 0..3 {
            if alignment_length < offset + 3 {
                continue;
            }
            let codon_count = (alignment_length - offset) / 3;
            if codon_count < minimum_codons {
                continue;
            }

            let mut stop_prefixes = Vec::with_capacity(oriented.len());
            let mut resolved_prefixes = Vec::with_capacity(oriented.len());
            let mut candidate_intervals: HashSet<(usize, usize)> = HashSet::new();

            for sequence in oriented {
                let bytes = sequence.as_bytes();
                let mut stop_prefix = vec![0usize; codon_count + 1];
                let mut resolved_prefix = vec![0usize; codon_count + 1];
                let mut run_start = 0usize;

                for codon_index in 0..codon_count {
                    let start = offset + codon_index * 3;
                    let amino_acid = translate_codon(&bytes[start..start + 3], code);
                    let is_stop = amino_acid == '*';
                    let is_resolved = !matches!(amino_acid, '*' | 'X' | '-');
                    stop_prefix[codon_index + 1] = stop_prefix[codon_index] + usize::from(is_stop);
                    resolved_prefix[codon_index + 1] =
                        resolved_prefix[codon_index] + usize::from(is_resolved);

                    if is_stop {
                        if codon_index.saturating_sub(run_start) >= minimum_codons {
                            candidate_intervals.insert((run_start, codon_index));
                        }
                        run_start = codon_index + 1;
                    }
                }

                if codon_count.saturating_sub(run_start) >= minimum_codons {
                    candidate_intervals.insert((run_start, codon_count));
                }
                stop_prefixes.push(stop_prefix);
                resolved_prefixes.push(resolved_prefix);
            }

            for (start_codon, end_codon) in candidate_intervals {
                let length_codons = end_codon - start_codon;
                let minimum_resolved = ((length_codons as f64) * 0.25).ceil() as usize;
                let mut support_count = 0usize;
                let mut informative_codons = 0usize;

                for (stop_prefix, resolved_prefix) in
                    stop_prefixes.iter().zip(resolved_prefixes.iter())
                {
                    let stops = stop_prefix[end_codon] - stop_prefix[start_codon];
                    let resolved = resolved_prefix[end_codon] - resolved_prefix[start_codon];
                    if stops == 0 && resolved >= minimum_resolved.max(1) {
                        support_count += 1;
                        informative_codons += resolved;
                    }
                }

                if support_count < min_support {
                    continue;
                }

                let candidate = SharedOrfSegmentResult {
                    frame: if is_reverse {
                        -(offset as i8 + 1)
                    } else {
                        offset as i8 + 1
                    },
                    is_reverse,
                    start: offset + start_codon * 3,
                    end: offset + end_codon * 3,
                    support_count,
                    length_codons,
                    informative_codons,
                    coding_score: 0.0,
                    amino_acid_conservation: 0.0,
                    frame_contrast: 0.0,
                };
                candidates.push(candidate);
                candidates.sort_by(|left, right| {
                    let left_score = left.support_count * left.length_codons;
                    let right_score = right.support_count * right.length_codons;
                    right_score
                        .cmp(&left_score)
                        .then_with(|| right.support_count.cmp(&left.support_count))
                        .then_with(|| right.length_codons.cmp(&left.length_codons))
                        .then_with(|| right.informative_codons.cmp(&left.informative_codons))
                });
                candidates.truncate(64);
            }
        }
    }

    for candidate in &mut candidates {
        let oriented = if candidate.is_reverse {
            &reverse_sequences
        } else {
            sequences
        };
        let (coding_score, amino_acid_conservation, frame_contrast) =
            candidate_coding_evidence(oriented, candidate.start, candidate.end, code);
        candidate.coding_score = coding_score;
        candidate.amino_acid_conservation = amino_acid_conservation;
        candidate.frame_contrast = frame_contrast;
    }

    if let Some(index) = candidates
        .iter()
        .position(|candidate| candidate.coding_score >= min_coding_score)
    {
        return Some(candidates.remove(index));
    }

    candidates.into_iter().max_by(|left, right| {
        left.coding_score
            .partial_cmp(&right.coding_score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                (left.support_count * left.length_codons)
                    .cmp(&(right.support_count * right.length_codons))
            })
    })
}

/// Executes full Open Reading Frame optimization and stop codon quality control
pub fn optimize_open_reading_frames(
    taxa: &[String],
    sequences: &[String],
    locus_id: &str,
    config: &OrfConfig,
) -> OrfOptimizationResult {
    optimize_open_reading_frames_guided(taxa, sequences, locus_id, config, None)
}

/// Runs ORF optimization with an optional frame supplied by a matched coding
/// reference. The reference fixes strand and phase; stop-codon and coding
/// evidence checks still validate the alignment itself.
pub fn optimize_open_reading_frames_guided(
    taxa: &[String],
    sequences: &[String],
    locus_id: &str,
    config: &OrfConfig,
    reference_frame: Option<i8>,
) -> OrfOptimizationResult {
    if !config.enable_orf || sequences.is_empty() || taxa.is_empty() {
        return OrfOptimizationResult {
            taxa: taxa.to_vec(),
            sequences: sequences.to_vec(),
            dropped_taxa: Vec::new(),
            trimmed_columns: Vec::new(),
            masked_segments: Vec::new(),
            found_valid_orf: true,
            is_reverse: false,
            orf_evaluated: false,
            candidate_found: false,
            candidate_frame: None,
            candidate_start: None,
            candidate_end: None,
            candidate_support_count: 0,
            candidate_support_percent: 0.0,
            candidate_length_aa: 0,
            coding_score: 0.0,
            amino_acid_conservation: 0.0,
            frame_contrast: 0.0,
            retained_samples: taxa.len(),
        };
    }

    // Check mode-aware non-coding exclusion.
    if config.exclude_uce && should_skip_orf_locus(locus_id, config.search_mode) {
        return OrfOptimizationResult {
            taxa: taxa.to_vec(),
            sequences: sequences.to_vec(),
            dropped_taxa: Vec::new(),
            trimmed_columns: Vec::new(),
            masked_segments: Vec::new(),
            found_valid_orf: true,
            is_reverse: false,
            orf_evaluated: false,
            candidate_found: false,
            candidate_frame: None,
            candidate_start: None,
            candidate_end: None,
            candidate_support_count: 0,
            candidate_support_percent: 0.0,
            candidate_length_aa: 0,
            coding_score: 0.0,
            amino_acid_conservation: 0.0,
            frame_contrast: 0.0,
            retained_samples: taxa.len(),
        };
    }

    let mut current_taxa = taxa.to_vec();
    let mut current_seqs = sequences.to_vec();
    let mut dropped_taxa = Vec::new();
    let mut trimmed_columns = Vec::new();
    let mut masked_segments = Vec::new();

    // 1. Detect either a continuous CDS frame or a candidate shared segment,
    // retaining the evidence needed to audit the decision in the catalog.
    let mut is_reverse = false;
    let mut trim_start = 0usize;
    let mut trim_end = current_seqs.first().map_or(0, String::len);
    let mut found_valid_orf = false;
    let mut candidate_found = false;
    let mut candidate_frame = None;
    let mut candidate_support_count = 0usize;
    let mut candidate_length_aa = 0usize;
    let mut coding_score = 0.0;
    let mut amino_acid_conservation = 0.0;
    let mut frame_contrast = 0.0;

    if let Some(frame_hint) = reference_frame {
        let frame = evaluate_reading_frame(&current_seqs, config.genetic_code, frame_hint);
        let old_len = current_seqs.first().map_or(0, String::len);
        let trailing_remainder = old_len.saturating_sub(frame.offset) % 3;
        is_reverse = frame.is_reverse;
        trim_start = frame.offset;
        trim_end = old_len.saturating_sub(trailing_remainder);
        candidate_found = frame.clean_taxa_count > 0;
        found_valid_orf = candidate_found;
        candidate_frame = Some(frame.frame);
        candidate_support_count = frame.clean_taxa_count;
        candidate_length_aa = trim_end.saturating_sub(trim_start) / 3;

        let reverse_sequences: Vec<String> = current_seqs
            .iter()
            .map(|sequence| reverse_complement_dna(sequence))
            .collect();
        let oriented = if is_reverse {
            &reverse_sequences
        } else {
            &current_seqs
        };
        (coding_score, amino_acid_conservation, frame_contrast) =
            candidate_coding_evidence(oriented, trim_start, trim_end, config.genetic_code);
    } else {
      match config.search_mode {
        OrfSearchMode::ContinuousCds | OrfSearchMode::ReferenceGuided => {
            let frame = find_optimal_reading_frame(&current_seqs, config.genetic_code);
            let old_len = current_seqs.first().map_or(0, String::len);
            let trailing_remainder = old_len.saturating_sub(frame.offset) % 3;
            is_reverse = frame.is_reverse;
            trim_start = frame.offset;
            trim_end = old_len.saturating_sub(trailing_remainder);
            candidate_found = frame.clean_taxa_count > 0;
            found_valid_orf = candidate_found;
            candidate_frame = Some(frame.frame);
            candidate_support_count = frame.clean_taxa_count;
            candidate_length_aa = trim_end.saturating_sub(trim_start) / 3;

            let reverse_sequences: Vec<String> = current_seqs
                .iter()
                .map(|sequence| reverse_complement_dna(sequence))
                .collect();
            let oriented = if is_reverse {
                &reverse_sequences
            } else {
                &current_seqs
            };
            (coding_score, amino_acid_conservation, frame_contrast) =
                candidate_coding_evidence(oriented, trim_start, trim_end, config.genetic_code);
        }
        OrfSearchMode::BestSharedSegment | OrfSearchMode::ReferenceCandidateOrf => {
            if let Some(segment) = find_best_shared_orf_segment(
                &current_seqs,
                config.genetic_code,
                config.min_shared_support_percent,
                config.min_segment_aa,
                config.min_coding_score,
            ) {
                is_reverse = segment.is_reverse;
                trim_start = segment.start;
                trim_end = segment.end;
                candidate_found = true;
                candidate_frame = Some(segment.frame);
                candidate_support_count = segment.support_count;
                candidate_length_aa = segment.length_codons;
                coding_score = segment.coding_score;
                amino_acid_conservation = segment.amino_acid_conservation;
                frame_contrast = segment.frame_contrast;
                found_valid_orf = coding_score >= config.min_coding_score;
            }
        }
      }
    }

    // A reverse-strand result is only usable when the caller permits the
    // alignment to be reoriented. Otherwise leave the locus unchanged and let
    // fail_if_no_orf report that no usable ORF was found.
    if is_reverse && !config.auto_flip_reverse {
        found_valid_orf = false;
    }

    // 2. Flip alignment if on reverse strand
    if is_reverse && config.auto_flip_reverse && found_valid_orf {
        current_seqs = current_seqs
            .into_iter()
            .map(|s| reverse_complement_dna(&s))
            .collect();
    }

    // 3. Trim leading offset and trailing partial codons to lock triplet boundaries
    let should_trim = reference_frame.is_some()
        || config.auto_shift_frame
        || matches!(
            config.search_mode,
            OrfSearchMode::BestSharedSegment | OrfSearchMode::ReferenceCandidateOrf
        );
    if should_trim && found_valid_orf && !current_seqs.is_empty() {
        let old_len = current_seqs[0].len();
        let start = trim_start.min(old_len);
        let end = trim_end.min(old_len).max(start);
        let new_len = end - start;

        // Record trimmed leading columns
        for col in 0..start {
            trimmed_columns.push(col);
        }
        // Record trimmed trailing columns
        for col in end..old_len {
            trimmed_columns.push(col);
        }

        current_seqs = current_seqs
            .into_iter()
            .map(|s| {
                if s.len() >= start + new_len {
                    s[start..start + new_len].to_string()
                } else {
                    String::new()
                }
            })
            .collect();
    }

    // 4. Codon QC: Terminal MACSE, Internal MACSE, and Stop Codon filtering
    let align_len = current_seqs.first().map_or(0, |s| s.len());
    let codon_count = align_len / 3;

    if found_valid_orf && codon_count > 0 {
        let mut kept_taxa = Vec::new();
        let mut kept_seqs = Vec::new();
        let mut total_locus_macse_cols = vec![false; align_len];
        let mut total_locus_stop_codons = 0;

        for (t_idx, taxon) in current_taxa.iter().enumerate() {
            let seq = &current_seqs[t_idx];
            let mut seq_bytes = seq.as_bytes().to_vec();
            
            // 4a. Trim Terminal MACSE Frameshifts (!)
            if config.macse_trim_terminal {
                // Trim 5'
                for c_idx in 0..codon_count {
                    let start = c_idx * 3;
                    if start + 3 > seq_bytes.len() { break; }
                    let codon = &seq_bytes[start..start + 3];
                    if codon.contains(&b'!') {
                        seq_bytes[start] = b'-';
                        seq_bytes[start + 1] = b'-';
                        seq_bytes[start + 2] = b'-';
                        masked_segments.push(MaskedSegment {
                            taxon: taxon.clone(),
                            start,
                            end: start + 3,
                        });
                    } else {
                        break;
                    }
                }
                // Trim 3'
                for c_idx in (0..codon_count).rev() {
                    let start = c_idx * 3;
                    if start + 3 > seq_bytes.len() { break; }
                    let codon = &seq_bytes[start..start + 3];
                    if codon.contains(&b'!') {
                        seq_bytes[start] = b'-';
                        seq_bytes[start + 1] = b'-';
                        seq_bytes[start + 2] = b'-';
                        masked_segments.push(MaskedSegment {
                            taxon: taxon.clone(),
                            start,
                            end: start + 3,
                        });
                    } else {
                        break;
                    }
                }
            }

            let mut sample_stop_count = 0;
            let mut sample_internal_macse = 0;
            let mut has_internal_stop = false;

            for c_idx in 0..codon_count {
                let start = c_idx * 3;
                if start + 3 > seq_bytes.len() {
                    break;
                }
                
                let codon = &seq_bytes[start..start + 3];
                
                // Count internal MACSE
                if codon.contains(&b'!') {
                    let macse_count = codon.iter().filter(|&&b| b == b'!').count();
                    sample_internal_macse += macse_count;
                    for i in 0..3 {
                        if codon[i] == b'!' {
                            total_locus_macse_cols[start + i] = true;
                        }
                    }
                }

                // Stop Codon checks
                if !codon.contains(&b'-') && !codon.contains(&b'N') && !codon.contains(&b'?') && !codon.contains(&b'!') {
                    let aa = translate_codon(codon, config.genetic_code);
                    if aa == '*' && c_idx < codon_count - 1 {
                        sample_stop_count += 1;
                        total_locus_stop_codons += 1;
                        has_internal_stop = true;

                        if config.stop_codon_action == StopCodonAction::MaskCodon {
                            seq_bytes[start] = b'-';
                            seq_bytes[start + 1] = b'-';
                            seq_bytes[start + 2] = b'-';
                            masked_segments.push(MaskedSegment {
                                taxon: taxon.clone(),
                                start,
                                end: start + 3,
                            });
                        }
                    }
                }
            }

            // Rejection checks
            let exceeds_macse = sample_internal_macse > config.macse_max_internal_sample;
            let exceeds_stop = sample_stop_count > config.max_stop_codons_sample;
            let legacy_remove = has_internal_stop && config.stop_codon_action == StopCodonAction::RemoveSample;

            if exceeds_macse || exceeds_stop || legacy_remove {
                dropped_taxa.push(taxon.clone());
            } else {
                kept_taxa.push(taxon.clone());
                kept_seqs.push(String::from_utf8(seq_bytes).unwrap());
            }
        }

        let locus_macse_cols = total_locus_macse_cols.iter().filter(|&&b| b).count();
        if locus_macse_cols > config.macse_max_internal_locus || total_locus_stop_codons > config.max_stop_codons_locus {
            // Drop entire locus
            return OrfOptimizationResult {
                taxa: vec![],
                sequences: vec![],
                dropped_taxa: current_taxa, // all dropped
                trimmed_columns: vec![],
                masked_segments: vec![],
                found_valid_orf: false, // failed QC
                is_reverse,
                orf_evaluated: true,
                candidate_found,
                candidate_frame,
                candidate_start: None,
                candidate_end: None,
                candidate_support_count,
                candidate_support_percent: 0.0,
                candidate_length_aa: 0,
                coding_score,
                amino_acid_conservation,
                frame_contrast,
                retained_samples: 0,
            };
        }

        current_taxa = kept_taxa;
        current_seqs = kept_seqs;
    }

    let retained_samples = current_taxa.len();
    OrfOptimizationResult {
        taxa: current_taxa,
        sequences: current_seqs,
        dropped_taxa,
        trimmed_columns,
        masked_segments,
        found_valid_orf,
        is_reverse,
        orf_evaluated: true,
        candidate_found,
        candidate_frame,
        candidate_start: candidate_found.then_some(trim_start),
        candidate_end: candidate_found.then_some(trim_end),
        candidate_support_count,
        candidate_support_percent: if taxa.is_empty() {
            0.0
        } else {
            candidate_support_count as f64 / taxa.len() as f64 * 100.0
        },
        candidate_length_aa,
        coding_score,
        amino_acid_conservation,
        frame_contrast,
        retained_samples,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement_dna("ATGC--N"), "N--GCAT");
    }

    #[test]
    fn test_detect_forward_offset_frame() {
        // Starts with 1 extra leading base 'G', then pure codons: ATG (Met), TTT (Phe), GGG (Gly), TAA (Stop)
        let seqs = vec![
            "GATGTTTGGGTAA".to_string(),
            "GATGTTTGGGTAA".to_string(),
            "GATGTTTGGGTAA".to_string(),
        ];
        let res = find_optimal_reading_frame(&seqs, GeneticCode::Standard);
        assert_eq!(res.frame, 2);
        assert_eq!(res.offset, 1);
        assert!(!res.is_reverse);
    }

    #[test]
    fn test_detect_reverse_complement_frame() {
        // Reverse complement of: ATGTTTGGGTAA (len 12) -> TTACCCAAACAT
        let seqs = vec![
            "TTACCCAAACAT".to_string(),
            "TTACCCAAACAT".to_string(),
            "TTACCCAAACAT".to_string(),
        ];
        let res = find_optimal_reading_frame(&seqs, GeneticCode::Standard);
        assert!(res.is_reverse);
        assert_eq!(res.offset, 0);
    }

    #[test]
    fn test_prune_sample_with_premature_stop() {
        let taxa = vec![
            "Taxon_1".to_string(),
            "Taxon_2".to_string(),
            "Taxon_Pseudogene".to_string(),
        ];
        // 9 bp: ATG AAA GGG (Taxon 1 & 2) vs ATG TAA GGG (Taxon_Pseudogene has premature TAA stop at codon 2)
        let seqs = vec![
            "ATGAAAGGG".to_string(),
            "ATGAAAGGG".to_string(),
            "ATGTAAGGG".to_string(),
        ];

        let config = OrfConfig {
            enable_orf: true,
            auto_shift_frame: true,
            auto_flip_reverse: true,
            stop_codon_action: StopCodonAction::RemoveSample,
            genetic_code: GeneticCode::Standard,
            search_mode: OrfSearchMode::ContinuousCds,
            min_shared_support_percent: 60.0,
            min_segment_aa: 3,
            min_coding_score: 0.0,
            exclude_uce: false,
            fail_if_no_orf: false, max_stop_codons_sample: 2, max_stop_codons_locus: 5, macse_trim_terminal: true, macse_max_internal_sample: 3, macse_max_internal_locus: 10,
            max_stop_codons_sample: 2,
            max_stop_codons_locus: 5,
            macse_trim_terminal: true,
            macse_max_internal_sample: 3,
            macse_max_internal_locus: 10,
        };

        let result = optimize_open_reading_frames(&taxa, &seqs, "exon_001", &config);

        assert_eq!(result.taxa.len(), 2);
        assert_eq!(result.dropped_taxa, vec!["Taxon_Pseudogene"]);
        assert_eq!(result.sequences[0], "ATGAAAGGG");
    }

    #[test]
    fn test_skip_non_coding_and_supercontig_loci() {
        let taxa = vec!["Taxon_1".to_string()];
        let seqs = vec!["GATGTTTGGGTAA".to_string()];

        let config = OrfConfig {
            enable_orf: true,
            auto_shift_frame: true,
            auto_flip_reverse: true,
            stop_codon_action: StopCodonAction::RemoveSample,
            genetic_code: GeneticCode::Standard,
            search_mode: OrfSearchMode::ContinuousCds,
            min_shared_support_percent: 60.0,
            min_segment_aa: 3,
            min_coding_score: 0.0,
            exclude_uce: true,
            fail_if_no_orf: false, max_stop_codons_sample: 2, max_stop_codons_locus: 5, macse_trim_terminal: true, macse_max_internal_sample: 3, macse_max_internal_locus: 10,
        };

        for locus_id in [
            "uce-1048",
            "gene_1_intron",
            "gene_2_supercontig",
            "gene_3_flanking",
            "gene_4_intergenic",
            "gene_5_non-coding",
        ] {
            let result = optimize_open_reading_frames(&taxa, &seqs, locus_id, &config);

            assert_eq!(result.taxa.len(), 1, "{locus_id}");
            assert_eq!(result.sequences[0], "GATGTTTGGGTAA", "{locus_id}");
            assert!(result.trimmed_columns.is_empty(), "{locus_id}");
        }
    }

    #[test]
    fn test_best_shared_mode_processes_supercontigs_but_still_skips_noncoding_loci() {
        assert!(should_skip_orf_locus(
            "uce-1048",
            OrfSearchMode::BestSharedSegment
        ));
        assert!(should_skip_orf_locus(
            "gene_intergenic",
            OrfSearchMode::BestSharedSegment
        ));
        assert!(!should_skip_orf_locus(
            "gene_supercontig",
            OrfSearchMode::BestSharedSegment
        ));
        assert!(should_skip_orf_locus(
            "gene_with_intron",
            OrfSearchMode::BestSharedSegment
        ));

        let taxa = vec!["Taxon_1".to_string(), "Taxon_2".to_string()];
        let seqs = vec![
            "ATGATGATGATGA".to_string(),
            "ATGATGATGATGA".to_string(),
        ];
        let config = OrfConfig {
            enable_orf: true,
            auto_shift_frame: true,
            auto_flip_reverse: true,
            stop_codon_action: StopCodonAction::RemoveSample,
            genetic_code: GeneticCode::Standard,
            search_mode: OrfSearchMode::BestSharedSegment,
            min_shared_support_percent: 100.0,
            min_segment_aa: 2,
            min_coding_score: 0.0,
            exclude_uce: true,
            fail_if_no_orf: true, max_stop_codons_sample: 2, max_stop_codons_locus: 5, macse_trim_terminal: true, macse_max_internal_sample: 3, macse_max_internal_locus: 10,
        };

        let result =
            optimize_open_reading_frames(&taxa, &seqs, "gene_supercontig", &config);
        assert!(result.found_valid_orf);
        assert_eq!(result.sequences[0].len(), 12);
        assert_eq!(result.sequences[0].len() % 3, 0);
        assert_eq!(result.trimmed_columns.len(), 1);

        let skipped = optimize_open_reading_frames(&taxa, &seqs, "uce-1048", &config);
        assert_eq!(skipped.sequences, seqs);
        assert!(skipped.trimmed_columns.is_empty());
    }

    #[test]
    fn test_best_shared_mode_honors_minimum_segment_length() {
        let seqs = vec!["ATGATGATG".to_string(), "ATGATGATG".to_string()];
        let result = find_best_shared_orf_segment(
            &seqs,
            GeneticCode::Standard,
            100.0,
            4,
            0.0,
        );
        assert!(result.is_none());
    }

    #[test]
    fn test_reference_frame_hint_controls_phase_and_is_still_validated() {
        let taxa = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let sequences = vec!["AAA".repeat(36); taxa.len()];
        let config = OrfConfig {
            enable_orf: true,
            auto_shift_frame: true,
            auto_flip_reverse: true,
            stop_codon_action: StopCodonAction::RemoveSample,
            genetic_code: GeneticCode::Standard,
            search_mode: OrfSearchMode::ContinuousCds,
            min_shared_support_percent: 90.0,
            min_segment_aa: 35,
            min_coding_score: 40.0,
            exclude_uce: false,
            fail_if_no_orf: true, max_stop_codons_sample: 2, max_stop_codons_locus: 5, macse_trim_terminal: true, macse_max_internal_sample: 3, macse_max_internal_locus: 10,
        };

        let de_novo = optimize_open_reading_frames(&taxa, &sequences, "exon", &config);
        assert_eq!(de_novo.candidate_frame, Some(1));

        let guided = optimize_open_reading_frames_guided(
            &taxa,
            &sequences,
            "exon",
            &config,
            Some(2),
        );
        assert_eq!(guided.candidate_frame, Some(2));
        assert_eq!(guided.candidate_support_count, 3);
        assert_eq!(guided.candidate_length_aa, 35);
        assert_eq!(guided.sequences[0].len(), 105);
        assert!(guided.found_valid_orf);
    }

    #[test]
    fn test_best_shared_mode_does_not_count_missing_codons_as_support() {
        let seqs = vec![
            "------------------".to_string(),
            "NNNNNNNNNNNNNNNNNN".to_string(),
        ];
        let result =
            find_best_shared_orf_segment(&seqs, GeneticCode::Standard, 50.0, 3, 0.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_coding_evidence_accepts_conserved_protein_with_synonymous_changes() {
        let seqs = ["GCT", "GCC", "GCA", "GCG"]
            .iter()
            .map(|codon| codon.repeat(90))
            .collect::<Vec<_>>();
        let result = find_best_shared_orf_segment(
            &seqs,
            GeneticCode::Standard,
            90.0,
            75,
            40.0,
        )
        .expect("expected a supported candidate");

        assert!(result.coding_score >= 40.0, "score={}", result.coding_score);
        assert!(result.amino_acid_conservation > 0.95);
    }

    #[test]
    fn test_stop_free_sequence_alone_is_not_sufficient_coding_evidence() {
        let taxa = (1..=4).map(|index| format!("Taxon_{index}")).collect::<Vec<_>>();
        let seqs = vec!["ATG".repeat(90); taxa.len()];
        let config = OrfConfig {
            enable_orf: true,
            auto_shift_frame: true,
            auto_flip_reverse: true,
            stop_codon_action: StopCodonAction::RemoveSample,
            genetic_code: GeneticCode::Standard,
            search_mode: OrfSearchMode::BestSharedSegment,
            min_shared_support_percent: 90.0,
            min_segment_aa: 75,
            min_coding_score: 40.0,
            exclude_uce: false,
            fail_if_no_orf: true, max_stop_codons_sample: 2, max_stop_codons_locus: 5, macse_trim_terminal: true, macse_max_internal_sample: 3, macse_max_internal_locus: 10,
        };

        let result = optimize_open_reading_frames(&taxa, &seqs, "exon_stop_free", &config);
        assert!(result.candidate_found);
        assert!(result.coding_score < 40.0, "score={}", result.coding_score);
        assert!(!result.found_valid_orf);
        assert_eq!(result.sequences, seqs, "rejected candidates must not trim the alignment");
    }

    #[test]
    fn test_homologous_noncoding_alignments_do_not_pass_coding_evidence_by_chance() {
        fn next_random(state: &mut u64) -> u64 {
            *state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            *state
        }

        let mut state = 0xA11C_EF0Au64;
        let mut accepted = 0usize;
        for _locus in 0..20 {
            let ancestor: Vec<u8> = (0..900)
                .map(|_| b"ACGT"[(next_random(&mut state) % 4) as usize])
                .collect();
            let sequences = (0..8)
                .map(|_| {
                    ancestor
                        .iter()
                        .map(|base| {
                            if next_random(&mut state) % 100 < 5 {
                                let replacement_index = (next_random(&mut state) % 4) as usize;
                                let mut replacement = b"ACGT"[replacement_index];
                                if replacement == *base {
                                    replacement = b"ACGT"[(replacement_index + 1) % 4];
                                }
                                replacement
                            } else {
                                *base
                            }
                        })
                        .map(char::from)
                        .collect::<String>()
                })
                .collect::<Vec<_>>();
            if find_best_shared_orf_segment(
                &sequences,
                GeneticCode::Standard,
                90.0,
                75,
                40.0,
            )
            .is_some_and(|candidate| candidate.coding_score >= 40.0)
            {
                accepted += 1;
            }
        }

        assert!(accepted <= 2, "accepted {accepted}/20 simulated non-coding loci");
    }
}
