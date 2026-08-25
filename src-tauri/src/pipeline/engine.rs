use std::collections::{HashMap, HashSet};
use crate::algorithms::assess::assess_alignment;
use crate::algorithms::columns::trim_alignment_columns;
use crate::algorithms::coverage::filter_sample_coverage;
use crate::algorithms::external::trim_external;
use crate::algorithms::hmm::clean_with_profile_hmm;
use crate::algorithms::informative::calculate_site_statistics;
use crate::algorithms::orf::{
    find_optimal_reading_frame, optimize_open_reading_frames_guided, should_skip_orf_locus,
    translate_codon, GeneticCode, OrfConfig, OrfSearchMode,
};
use crate::algorithms::reference::match_reference_to_alignment;
use crate::algorithms::sanitize::{
    convert_ambiguous_consensus, remove_gap_only_columns, replace_character,
};
use crate::algorithms::segments::mask_divergent_segments;
use crate::algorithms::similarity::filter_sample_similarity;
use crate::algorithms::statistical_columns::trim_statistical_columns;
use crate::algorithms::stats::calculate_gap_stats;
use crate::models::{Alignment, StopCodonPos, TrimmingDiff};
use crate::pipeline::recipe::TrimmingRecipe;

fn detect_stop_codons(
    taxa: &[String],
    sequences: &[String],
    genetic_code: GeneticCode,
) -> Vec<StopCodonPos> {
    let mut stops = Vec::new();
    for (taxon, sequence) in taxa.iter().zip(sequences.iter()) {
        let codon_count = sequence.len() / 3;
        for codon_index in 0..codon_count {
            let start = codon_index * 3;
            let codon = &sequence.as_bytes()[start..start + 3];
            if translate_codon(codon, genetic_code) == '*' {
                stops.push(StopCodonPos {
                    taxon: taxon.clone(),
                    start,
                    end: start + 3,
                    codon: String::from_utf8_lossy(codon).to_ascii_uppercase(),
                    is_terminal: codon_index + 1 == codon_count,
                });
            }
        }
    }
    stops
}

fn map_final_stops_to_raw(
    final_stops: &[StopCodonPos],
    raw_col_map: &[usize],
) -> Vec<StopCodonPos> {
    final_stops
        .iter()
        .filter_map(|stop| {
            let mapped = raw_col_map.get(stop.start..stop.end)?;
            let raw_start = *mapped.iter().min()?;
            let raw_end = mapped.iter().max()?.saturating_add(1);
            // Only annotate a raw triplet when the three retained columns are contiguous.
            if raw_end - raw_start != 3 {
                return None;
            }
            Some(StopCodonPos {
                taxon: stop.taxon.clone(),
                start: raw_start,
                end: raw_end,
                codon: stop.codon.clone(),
                is_terminal: stop.is_terminal,
            })
        })
        .collect()
}

fn detect_raw_stops_in_selected_frame(
    taxa: &[String],
    sequences: &[String],
    frame: i8,
    genetic_code: GeneticCode,
    raw_col_map: &[usize],
    selected_region_raw_map: &[usize],
) -> Vec<StopCodonPos> {
    if frame == 0 {
        return Vec::new();
    }
    let is_reverse = frame < 0;
    let frame_offset = frame.unsigned_abs().saturating_sub(1) as usize;
    let oriented_map: Vec<usize> = if is_reverse {
        raw_col_map.iter().rev().copied().collect()
    } else {
        raw_col_map.to_vec()
    };
    let oriented_sequences: Vec<String> = if is_reverse {
        sequences
            .iter()
            .map(|sequence| crate::algorithms::orf::reverse_complement_dna(sequence))
            .collect()
    } else {
        sequences.to_vec()
    };
    let selected_origin_column = if is_reverse {
        selected_region_raw_map.last()
    } else {
        selected_region_raw_map.first()
    };
    let selected_origin = selected_origin_column
        .and_then(|column| oriented_map.iter().position(|candidate| candidate == column))
        .unwrap_or(0);
    let offset = (selected_origin + frame_offset) % 3;
    let mut stops = Vec::new();

    for (taxon, sequence) in taxa.iter().zip(oriented_sequences.iter()) {
        if sequence.len() < offset + 3 {
            continue;
        }
        let codon_count = (sequence.len() - offset) / 3;
        for codon_index in 0..codon_count {
            let local_start = offset + codon_index * 3;
            let codon = &sequence.as_bytes()[local_start..local_start + 3];
            if translate_codon(codon, genetic_code) != '*' {
                continue;
            }
            let Some(mapped) = oriented_map.get(local_start..local_start + 3) else {
                continue;
            };
            let Some(raw_start) = mapped.iter().min().copied() else {
                continue;
            };
            let Some(raw_end_base) = mapped.iter().max().copied() else {
                continue;
            };
            let raw_end = raw_end_base + 1;
            if raw_end - raw_start == 3 {
                stops.push(StopCodonPos {
                    taxon: taxon.clone(),
                    start: raw_start,
                    end: raw_end,
                    codon: String::from_utf8_lossy(codon).to_ascii_uppercase(),
                    is_terminal: codon_index + 1 == codon_count,
                });
            }
        }
    }
    stops
}

/// Executes a non-destructive TrimmingRecipe against an Alignment.
/// Returns the transformed Alignment and the TrimmingDiff describing changes.
pub fn apply_recipe(alignment: &Alignment, recipe: &TrimmingRecipe, total_dataset_taxa: usize) -> (Alignment, TrimmingDiff) {
    let old_taxa_count = alignment.num_taxa;
    let old_length = alignment.length;

    let (_, _, old_gap_percent) = calculate_gap_stats(&alignment.sequences);
    let old_site_stats = calculate_site_statistics(&alignment.sequences, true);
    let old_variable = old_site_stats.variable_count;
    let old_pis = old_site_stats.pis_count;

    let mut current_taxa = alignment.taxa.clone();
    let mut current_seqs = alignment.sequences.clone();

    let mut raw_col_map: Vec<usize> = (0..old_length).collect();
    let mut column_reasons: HashMap<usize, String> = HashMap::new();
    let mut dropped_taxa_reasons: HashMap<String, String> = HashMap::new();
    let mut all_dropped_taxa = Vec::new();
    let mut all_masked_segments = Vec::new();
    let mut found_valid_orf = true;
    let mut orf_evaluated = false;
    let mut orf_candidate_found = false;
    let mut orf_frame = None;
    let mut orf_start = None;
    let mut orf_end = None;
    let mut orf_support_count = 0usize;
    let mut orf_support_percent = 0.0;
    let mut orf_retained_samples = 0usize;
    let mut orf_candidate_length_aa = 0usize;
    let mut orf_coding_score = 0.0;
    let mut orf_amino_acid_conservation = 0.0;
    let mut orf_frame_contrast = 0.0;
    let mut orf_reference_evaluated = false;
    let mut orf_reference_matched = false;
    let mut orf_reference_identity = 0.0;
    let mut orf_reference_coverage = 0.0;
    let mut orf_intron_length = 0usize;
    let mut reference_frame_hint = None;
    let mut raw_orf_stop_codons: Option<Vec<StopCodonPos>> = None;

    // Step 0: Dataset-wide sample occupancy filter. The catalog prepares this
    // derived list once so every locus, preview, and export drops the same taxa.
    if !recipe.excluded_taxa.is_empty() {
        let excluded: HashSet<&str> = recipe.excluded_taxa.iter().map(String::as_str).collect();
        let mut kept_taxa = Vec::with_capacity(current_taxa.len());
        let mut kept_seqs = Vec::with_capacity(current_seqs.len());

        for (taxon, sequence) in current_taxa.into_iter().zip(current_seqs.into_iter()) {
            if excluded.contains(taxon.as_str()) {
                dropped_taxa_reasons.insert(
                    taxon.clone(),
                    format!(
                        "Dataset-wide sample occupancy (< {:.0}% of loci)",
                        recipe.min_sample_locus_occupancy_percent
                    ),
                );
                all_dropped_taxa.push(taxon);
            } else {
                kept_taxa.push(taxon);
                kept_seqs.push(sequence);
            }
        }

        current_taxa = kept_taxa;
        current_seqs = kept_seqs;
    }

    // Step 1: Character Sanitation. Character conversions happen before sample
    // filtering, but no column is removed until all early sample filters finish.
    current_seqs = convert_ambiguous_consensus(&current_seqs, recipe.ambiguity_strategy);
    if recipe.replace_n_with_gap {
        current_seqs = replace_character(&current_seqs, 'N', '-');
        current_seqs = replace_character(&current_seqs, '?', '-');
    }

    // Step 2: Sample-level filtering. Remove samples before any column-level
    // operation so excluded sequences cannot influence gap, edge, or
    // statistical-column decisions.
    if recipe.trim_coverage && current_seqs.len() > 2 {
        let (cov_taxa, cov_seqs, dropped) = filter_sample_coverage(
            &current_taxa,
            &current_seqs,
            recipe.min_coverage_bp,
            recipe.min_coverage_percent,
            recipe.relative_width,
        );
        for t in &dropped {
            dropped_taxa_reasons.insert(
                t.clone(),
                format!("Low coverage (< {} bp / {:.0}%)", recipe.min_coverage_bp, recipe.min_coverage_percent),
            );
        }
        all_dropped_taxa.extend(dropped);
        current_taxa = cov_taxa;
        current_seqs = cov_seqs;
    }

    if recipe.trim_similarity && current_seqs.len() > 2 {
        let (kept_taxa, kept_seqs, dropped) =
            filter_sample_similarity(&current_taxa, &current_seqs, recipe.similarity_threshold);
        for t in &dropped {
            dropped_taxa_reasons.insert(
                t.clone(),
                format!("Divergent outlier (> {:.0}% divergence)", recipe.similarity_threshold * 100.0),
            );
        }
        all_dropped_taxa.extend(dropped);
        current_taxa = kept_taxa;
        current_seqs = kept_seqs;
    }

    // Step 3: Remove missing-only columns using only the surviving samples.
    if recipe.remove_gap_only_columns && !current_seqs.is_empty() {
        let (cleaned_seqs, dropped_columns) = remove_gap_only_columns(&current_seqs);
        if !dropped_columns.is_empty() {
            let dropped_set: HashSet<usize> = dropped_columns.iter().copied().collect();
            for local_column in &dropped_columns {
                if let Some(raw_column) = raw_col_map.get(*local_column) {
                    column_reasons.insert(
                        *raw_column,
                        "Gap-Only / Missing-Only Sanitation".to_string(),
                    );
                }
            }
            raw_col_map = raw_col_map
                .into_iter()
                .enumerate()
                .filter_map(|(local_column, raw_column)| {
                    (!dropped_set.contains(&local_column)).then_some(raw_column)
                })
                .collect();
            current_seqs = cleaned_seqs;
        }
    }

    // Preserve every original sample on the sanitized column coordinate system
    // so the raw overlay can show stops in trimmed regions and ORF-pruned taxa.
    let raw_orf_taxa = alignment.taxa.clone();
    let raw_orf_sequences: Vec<String> = alignment
        .sequences
        .iter()
        .map(|sequence| {
            raw_col_map
                .iter()
                .filter_map(|column| sequence.as_bytes().get(*column).copied())
                .map(char::from)
                .collect()
        })
        .collect();
    let raw_orf_col_map = raw_col_map.clone();

    let skip_orf = should_skip_orf_locus(&alignment.id, recipe.orf_search_mode);
    let is_coding_with_orf = recipe.enable_orf && !(recipe.exclude_uce && skip_orf);
    let uses_reference_mode = matches!(
        recipe.orf_search_mode,
        OrfSearchMode::ReferenceGuided | OrfSearchMode::ReferenceCandidateOrf
    );
    let allows_reference_fallback =
        recipe.orf_search_mode == OrfSearchMode::ReferenceCandidateOrf;

    // Keep the fully sanitized alignment so the hybrid mode can retry the
    // original region when a matched reference does not yield a usable ORF.
    let pre_reference_seqs = current_seqs.clone();
    let pre_reference_raw_col_map = raw_col_map.clone();
    let pre_reference_column_reasons = column_reasons.clone();
    let pre_reference_masked_segment_count = all_masked_segments.len();

    // Optional exact-name reference anchoring. The reference identifies the
    // exon span before coding-specific processing; excluded flanks are retained
    // independently during batch intron export.
    if is_coding_with_orf && uses_reference_mode && !current_seqs.is_empty() {
        orf_reference_evaluated = true;
        if let Some(reference) = recipe.orf_reference_sequences.get(&alignment.id) {
            if let Some(reference_match) = match_reference_to_alignment(&current_seqs, reference) {
                let current_length = current_seqs.first().map_or(0, String::len);
                let start = reference_match.start.min(current_length);
                let end = reference_match.end.min(current_length).max(start);
                if end > start {
                    orf_reference_matched = true;
                    orf_reference_identity = reference_match.identity_percent;
                    orf_reference_coverage = reference_match.coverage_percent;
                    orf_intron_length = current_length.saturating_sub(end - start);
                    let reference_frame =
                        find_optimal_reading_frame(&[reference.clone()], recipe.genetic_code).frame;
                    reference_frame_hint = Some(if reference_match.is_reverse {
                        -reference_frame
                    } else {
                        reference_frame
                    });
                    for local_column in 0..start {
                        if let Some(raw_column) = raw_col_map.get(local_column) {
                            column_reasons
                                .insert(*raw_column, "Reference-Anchored Intron".to_string());
                        }
                    }
                    for local_column in end..current_length {
                        if let Some(raw_column) = raw_col_map.get(local_column) {
                            column_reasons
                                .insert(*raw_column, "Reference-Anchored Intron".to_string());
                        }
                    }
                    raw_col_map = raw_col_map[start..end].to_vec();
                    current_seqs = current_seqs
                        .into_iter()
                        .map(|sequence| sequence.get(start..end).unwrap_or_default().to_string())
                        .collect();
                }
            }
        }
    }

    // Step 3a: Profile HMM Segment Cleaner (TAPIR-Style)
    if recipe.trim_hmm && current_seqs.len() > 2 {
        let (hmm_seqs, segments) = clean_with_profile_hmm(
            &current_taxa,
            &current_seqs,
            recipe.hmm_min_posterior,
            recipe.hmm_min_segment_length,
            recipe.hmm_min_island_length,
        );
        all_masked_segments.extend(segments);
        current_seqs = hmm_seqs;
    }

    // Step 3b: Sliding Window Segment Masking
    if recipe.trim_segments && current_seqs.len() > 2 {
        let (masked_seqs, segments) = mask_divergent_segments(
            &current_taxa,
            &current_seqs,
            recipe.segment_window_size,
            recipe.segment_threshold,
        );
        all_masked_segments.extend(segments);
        current_seqs = masked_seqs;
    }


    // Step 3c: Clean gap-only columns created by HMM or Segment masking
    if recipe.remove_gap_only_columns && !current_seqs.is_empty() {
        let (cleaned_seqs, dropped_columns) = remove_gap_only_columns(&current_seqs);
        if !dropped_columns.is_empty() {
            let dropped_set: std::collections::HashSet<usize> = dropped_columns.iter().copied().collect();
            for local_column in &dropped_columns {
                if let Some(raw_column) = raw_col_map.get(*local_column) {
                    column_reasons.insert(
                        *raw_column,
                        "Column masked completely by sequence trimming".to_string(),
                    );
                }
            }

            let mut new_raw_col_map = Vec::with_capacity(current_seqs[0].len());
            for idx in 0..current_seqs[0].len() {
                if !dropped_set.contains(&idx) {
                    new_raw_col_map.push(raw_col_map[idx]);
                }
            }
            raw_col_map = new_raw_col_map;
            current_seqs = cleaned_seqs;
        }
    }

    // Step 4: Open Reading Frame & Codon Optimization (Exons)
    let reference_blocks_orf = is_coding_with_orf
        && recipe.orf_search_mode == OrfSearchMode::ReferenceGuided
        && !orf_reference_matched;
    if reference_blocks_orf {
        found_valid_orf = false;
    }
    if recipe.enable_orf && !reference_blocks_orf && !current_seqs.is_empty() {
        let mut selected_region_raw_map = raw_col_map.clone();
        let orf_config = OrfConfig {
            enable_orf: recipe.enable_orf,
            auto_shift_frame: recipe.auto_shift_frame,
            auto_flip_reverse: recipe.auto_flip_reverse,
            stop_codon_action: recipe.stop_codon_action,
            genetic_code: recipe.genetic_code,
            search_mode: recipe.orf_search_mode,
            min_shared_support_percent: recipe.orf_min_shared_support_percent,
            min_segment_aa: recipe.orf_min_segment_aa,
            min_coding_score: recipe.orf_min_coding_score,
            exclude_uce: recipe.exclude_uce,
            max_stop_codons_sample: recipe.max_stop_codons_sample,
            max_stop_codons_locus: recipe.max_stop_codons_locus,
            macse_trim_terminal: recipe.macse_trim_terminal,
            macse_max_internal_sample: recipe.macse_max_internal_sample,
            macse_max_internal_locus: recipe.macse_max_internal_locus,
            fail_if_no_orf: recipe.fail_if_no_orf,
        };

        let mut orf_result = optimize_open_reading_frames_guided(
            &current_taxa,
            &current_seqs,
            &alignment.id,
            &orf_config,
            reference_frame_hint,
        );

        // A hybrid reference search must retry the unsliced alignment. This is
        // deliberately after reference-guided validation, not merely after a
        // failed sequence match, because a good match can still imply a frame
        // containing premature stops in the samples.
        let guided_has_internal_stop = detect_stop_codons(
            &orf_result.taxa,
            &orf_result.sequences,
            recipe.genetic_code,
        )
        .iter()
        .any(|stop| !stop.is_terminal);
        let guided_attempt_failed = !orf_result.found_valid_orf
            || orf_result.retained_samples == 0
            || guided_has_internal_stop;
        if allows_reference_fallback && orf_reference_matched && guided_attempt_failed {
            current_seqs = pre_reference_seqs.clone();
            raw_col_map = pre_reference_raw_col_map.clone();
            column_reasons = pre_reference_column_reasons.clone();
            all_masked_segments.truncate(pre_reference_masked_segment_count);

            if recipe.trim_hmm && current_seqs.len() > 2 {
                let (hmm_seqs, segments) = clean_with_profile_hmm(
                    &current_taxa,
                    &current_seqs,
                    recipe.hmm_min_posterior,
                    recipe.hmm_min_segment_length,
                    recipe.hmm_min_island_length,
                );
                all_masked_segments.extend(segments);
                current_seqs = hmm_seqs;
            }
            if recipe.trim_segments && current_seqs.len() > 2 {
                let (masked_seqs, segments) = mask_divergent_segments(
                    &current_taxa,
                    &current_seqs,
                    recipe.segment_window_size,
                    recipe.segment_threshold,
                );
                all_masked_segments.extend(segments);
                current_seqs = masked_seqs;
            }

            selected_region_raw_map = raw_col_map.clone();
            orf_result = optimize_open_reading_frames_guided(
                &current_taxa,
                &current_seqs,
                &alignment.id,
                &orf_config,
                None,
            );
        }
        found_valid_orf = orf_result.found_valid_orf;
        orf_evaluated = orf_result.orf_evaluated;
        orf_candidate_found = orf_result.candidate_found;
        orf_frame = orf_result.candidate_frame;
        orf_support_count = orf_result.candidate_support_count;
        orf_support_percent = orf_result.candidate_support_percent;
        orf_retained_samples = orf_result.retained_samples;
        orf_candidate_length_aa = orf_result.candidate_length_aa;
        orf_coding_score = orf_result.coding_score;
        orf_amino_acid_conservation = orf_result.amino_acid_conservation;
        orf_frame_contrast = orf_result.frame_contrast;

        let oriented_map: Vec<usize> = if orf_result.is_reverse {
            raw_col_map.iter().rev().copied().collect()
        } else {
            raw_col_map.clone()
        };
        if let (Some(start), Some(end)) =
            (orf_result.candidate_start, orf_result.candidate_end)
        {
            if let Some(mapped) = oriented_map.get(start..end) {
                orf_start = mapped.iter().min().copied();
                orf_end = mapped.iter().max().map(|column| column + 1);
            }
        }
        if let Some(frame) = orf_frame {
            raw_orf_stop_codons = Some(detect_raw_stops_in_selected_frame(
                &raw_orf_taxa,
                &raw_orf_sequences,
                frame,
                recipe.genetic_code,
                &raw_orf_col_map,
                &selected_region_raw_map,
            ));
        }

        if orf_result.is_reverse && recipe.auto_flip_reverse && orf_result.found_valid_orf {
            raw_col_map.reverse();
        }

        for t in &orf_result.dropped_taxa {
            dropped_taxa_reasons.insert(t.clone(), "Premature internal stop codon in exon".to_string());
        }
        all_dropped_taxa.extend(orf_result.dropped_taxa);
        all_masked_segments.extend(orf_result.masked_segments);

        let cur_len = raw_col_map.len();
        let mut kept_mask = vec![true; cur_len];
        for &local_col in &orf_result.trimmed_columns {
            if local_col < cur_len {
                kept_mask[local_col] = false;
                let raw_c = raw_col_map[local_col];
                column_reasons.insert(raw_c, "ORF Frame Shift / Codon Boundary".to_string());
            }
        }
        let mut new_raw_col_map = Vec::new();
        for (idx, &kept) in kept_mask.iter().enumerate() {
            if kept {
                new_raw_col_map.push(raw_col_map[idx]);
            }
        }
        raw_col_map = new_raw_col_map;
        current_taxa = orf_result.taxa;
        current_seqs = orf_result.sequences;
    }

    // Step 5: External Ragged Edge Trimming
    if recipe.trim_external && !current_seqs.is_empty() {
        let preserve_codon_frame = recipe.codon_preserving || is_coding_with_orf;
        let (ext_seqs, dropped_local_cols, _) = trim_external(
            &current_seqs,
            recipe.min_external_percent,
            preserve_codon_frame,
        );
        let cur_len = raw_col_map.len();
        let mut kept_mask = vec![true; cur_len];
        for &local_col in &dropped_local_cols {
            if local_col < cur_len {
                kept_mask[local_col] = false;
                let raw_c = raw_col_map[local_col];
                let reason = if local_col < cur_len / 2 {
                    format!("Ragged 5' End (< {:.0}% taxa coverage)", recipe.min_external_percent)
                } else {
                    format!("Ragged 3' End (< {:.0}% taxa coverage)", recipe.min_external_percent)
                };
                column_reasons.insert(raw_c, reason);
            }
        }
        let mut new_raw_col_map = Vec::new();
        for (idx, &kept) in kept_mask.iter().enumerate() {
            if kept {
                new_raw_col_map.push(raw_col_map[idx]);
            }
        }
        raw_col_map = new_raw_col_map;
        current_seqs = ext_seqs;
    }

    // Step 6: Column Gap Trimming (Bypassed on coding loci when ORF is enabled to preserve codon reading frames)
    if recipe.trim_columns && !is_coding_with_orf && !current_seqs.is_empty() {
        let (col_seqs, dropped_local_cols) =
            trim_alignment_columns(&current_seqs, recipe.min_column_gap_percent, recipe.count_n_as_gap);
        let cur_len = raw_col_map.len();
        let mut kept_mask = vec![true; cur_len];
        for &local_col in &dropped_local_cols {
            if local_col < cur_len {
                kept_mask[local_col] = false;
                let raw_c = raw_col_map[local_col];
                let mut gaps = 0usize;
                for s in &current_seqs {
                    if let Some(&b) = s.as_bytes().get(local_col) {
                        let u = b.to_ascii_uppercase();
                        if u == b'-' || u == b'?' || (recipe.count_n_as_gap && u == b'N') {
                            gaps += 1;
                        }
                    }
                }
                let gap_pct = if !current_seqs.is_empty() {
                    (gaps as f64 / current_seqs.len() as f64) * 100.0
                } else {
                    100.0
                };
                column_reasons.insert(
                    raw_c,
                    format!("High Gap Column ({:.1}% gaps > max {:.1}%)", gap_pct, recipe.min_column_gap_percent),
                );
            }
        }
        let mut new_raw_col_map = Vec::new();
        for (idx, &kept) in kept_mask.iter().enumerate() {
            if kept {
                new_raw_col_map.push(raw_col_map[idx]);
            }
        }
        raw_col_map = new_raw_col_map;
        current_seqs = col_seqs;
    }

    // Step 6b: Statistical Column Trimming (trimAl & Gblocks) (Bypassed on coding loci when ORF is enabled)
    if recipe.enable_statistical_columns && !is_coding_with_orf && !current_seqs.is_empty() {
        let (stat_seqs, dropped_local_cols, stat_reasons) = trim_statistical_columns(
            &current_seqs,
            recipe.stat_col_method,
            recipe.stat_col_similarity_threshold,
            recipe.stat_col_window_size,
            recipe.stat_col_heuristic,
            recipe.stat_col_min_block_length,
            recipe.stat_col_max_nonconserved,
            recipe.stat_col_gap_treatment,
            recipe.stat_col_entropy_threshold,
        );
        let cur_len = raw_col_map.len();
        let mut kept_mask = vec![true; cur_len];
        for &local_col in &dropped_local_cols {
            if local_col < cur_len {
                kept_mask[local_col] = false;
                let raw_c = raw_col_map[local_col];
                let reason = stat_reasons.get(&local_col).cloned().unwrap_or_else(|| {
                    "Statistical Column Quality Threshold".to_string()
                });
                column_reasons.insert(raw_c, reason);
            }
        }
        let mut new_raw_col_map = Vec::new();
        for (idx, &kept) in kept_mask.iter().enumerate() {
            if kept {
                new_raw_col_map.push(raw_col_map[idx]);
            }
        }
        raw_col_map = new_raw_col_map;
        current_seqs = stat_seqs;
    }

    let kept_set: HashSet<usize> = raw_col_map.iter().copied().collect();
    let mut all_trimmed_columns: Vec<usize> = (0..old_length).filter(|c| !kept_set.contains(c)).collect();
    all_trimmed_columns.sort();

    let new_taxa_count = current_taxa.len();
    let new_length = current_seqs.first().map_or(0, |s| s.len());

    let (_, _, new_gap_percent) = calculate_gap_stats(&current_seqs);
    let new_site_stats = calculate_site_statistics(&current_seqs, true);
    let new_variable = new_site_stats.variable_count;
    let new_variable_percent = new_site_stats.variable_percent;
    let new_pis = new_site_stats.pis_count;
    let new_pis_percent = new_site_stats.pis_percent;
    let final_stop_codons = detect_stop_codons(&current_taxa, &current_seqs, recipe.genetic_code);
    let stop_codons = raw_orf_stop_codons
        .unwrap_or_else(|| map_final_stops_to_raw(&final_stop_codons, &raw_col_map));

    // Step 8: Locus Assessment Gating
    let mut fail_reasons = Vec::new();

    if new_taxa_count == 0 {
        fail_reasons.push("0 surviving taxa (all samples pruned)".to_string());
    }

    if is_coding_with_orf {
        let has_internal_stop = final_stop_codons.iter().any(|stop| !stop.is_terminal);

    }

    if recipe.assess_alignment {
        let (_, assess_reasons) = assess_alignment(
            new_taxa_count,
            if total_dataset_taxa > 0 { total_dataset_taxa } else { alignment.num_taxa },
            new_length,
            new_gap_percent,
            recipe.min_taxa,
            recipe.min_taxa_occupancy_percent,
            recipe.min_length,
            recipe.max_gap_percent,
            new_pis,
            new_pis_percent,
            recipe.min_pis_count,
            recipe.min_pis_percent,
            new_variable,
            new_variable_percent,
            recipe.min_variable_count,
            recipe.min_variable_percent,
        );
        for r in assess_reasons {
            if !fail_reasons.contains(&r) {
                fail_reasons.push(r);
            }
        }
    }

    let pass = fail_reasons.is_empty();

    // Final Sanitization: Convert surviving MACSE frameshifts to standard missing data
    crate::algorithms::codon_qc::convert_macse_to_n(&mut current_seqs);

    let transformed_alignment = Alignment::new(
        alignment.id.clone(),
        alignment.file_name.clone(),
        alignment.file_path.clone(),
        alignment.format,
        current_taxa,
        current_seqs,
    );

    let diff = TrimmingDiff {
        id: alignment.id.clone(),
        old_taxa_count,
        new_taxa_count,
        dropped_taxa: all_dropped_taxa,
        old_length,
        new_length,
        trimmed_columns: all_trimmed_columns,
        masked_segments: all_masked_segments,
        column_reasons,
        dropped_taxa_reasons,
        stop_codons,
        final_stop_codons,
        old_gap_percent,
        new_gap_percent,
        old_variable,
        new_variable,
        old_pis,
        new_pis,
        found_valid_orf,
        orf_evaluated,
        orf_candidate_found,
        orf_frame,
        orf_start,
        orf_end,
        orf_support_count,
        orf_support_percent,
        orf_retained_samples,
        orf_candidate_length_aa,
        orf_coding_score,
        orf_amino_acid_conservation,
        orf_frame_contrast,
        orf_reference_evaluated,
        orf_reference_matched,
        orf_reference_identity,
        orf_reference_coverage,
        orf_intron_length,
        pass,
        fail_reasons,
    };

    (transformed_alignment, diff)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithms::orf::StopCodonAction;
    use crate::models::AlignmentFormat;

    #[test]
    fn test_apply_recipe_pipeline() {
        let taxa = vec![
            "Taxon_1".to_string(),
            "Taxon_2".to_string(),
            "Taxon_3".to_string(),
            "Taxon_4".to_string(),
        ];
        let seqs = vec![
            "--ATGCATGCATGC--".to_string(),
            "--ATGCATGCATGC--".to_string(),
            "--ATGCATGCATGC--".to_string(),
            "--ATGCATGCATGC--".to_string(),
        ];

        let align = Alignment::new(
            "locus_1".to_string(),
            "locus_1.fa".to_string(),
            "/dummy/locus_1.fa".to_string(),
            AlignmentFormat::Fasta,
            taxa,
            seqs,
        );

        let mut recipe = TrimmingRecipe::default();
        recipe.min_length = 10;
        let (transformed, diff) = apply_recipe(&align, &recipe, 0);

        assert_eq!(transformed.length, 12);
        assert_eq!(diff.new_length, 12);
        assert!(diff.pass);
        assert_eq!(diff.dropped_taxa.len(), 0);
    }

    #[test]
    fn test_sample_coverage_runs_before_column_filters() {
        let alignment = Alignment::new(
            "sample_first".to_string(),
            "sample_first.fa".to_string(),
            "/dummy/sample_first.fa".to_string(),
            AlignmentFormat::Fasta,
            vec![
                "LowCoverage".to_string(),
                "Taxon_1".to_string(),
                "Taxon_2".to_string(),
                "Taxon_3".to_string(),
            ],
            vec![
                "A---".to_string(),
                "-CGT".to_string(),
                "-CGT".to_string(),
                "-CGT".to_string(),
            ],
        );

        let mut recipe = TrimmingRecipe::default();
        recipe.trim_coverage = true;
        recipe.min_coverage_bp = 2;
        recipe.min_coverage_percent = 50.0;
        recipe.trim_similarity = false;
        recipe.remove_gap_only_columns = false;
        recipe.trim_external = false;
        recipe.trim_columns = true;
        recipe.min_column_gap_percent = 90.0;
        recipe.enable_statistical_columns = false;
        recipe.enable_orf = false;
        recipe.assess_alignment = false;

        let (transformed, diff) = apply_recipe(&alignment, &recipe, 4);

        assert_eq!(transformed.taxa, vec!["Taxon_1", "Taxon_2", "Taxon_3"]);
        assert_eq!(transformed.sequences, vec!["CGT", "CGT", "CGT"]);
        assert_eq!(diff.trimmed_columns, vec![0]);
        assert_eq!(
            diff.dropped_taxa_reasons.get("LowCoverage").map(String::as_str),
            Some("Low coverage (< 2 bp / 50%)")
        );
    }

    #[test]
    fn test_post_orf_edge_trim_preserves_frame_and_remaps_stops() {
        let alignment = Alignment::new(
            "exon_frame_test".to_string(),
            "exon_frame_test.fa".to_string(),
            "/dummy/exon_frame_test.fa".to_string(),
            AlignmentFormat::Fasta,
            vec![
                "Taxon_1".to_string(),
                "Taxon_2".to_string(),
                "Taxon_3".to_string(),
                "Taxon_4".to_string(),
            ],
            vec![
                "--GATGAAAGGGTAA".to_string(),
                "--GATGAAAGGGTAA".to_string(),
                "--GATGAAAGGGTAA".to_string(),
                "AAGATGAAAGGGTAA".to_string(),
            ],
        );

        let mut recipe = TrimmingRecipe::default();
        recipe.trim_similarity = false;
        recipe.trim_hmm = false;
        recipe.enable_orf = true;
        recipe.exclude_uce = false;
        recipe.stop_codon_action = StopCodonAction::RemoveSample;
        recipe.trim_external = true;
        recipe.min_external_percent = 75.0;
        recipe.codon_preserving = false; // ORF must enforce this downstream.
        recipe.trim_columns = false;
        recipe.trim_coverage = false;
        recipe.assess_alignment = false;

        let (transformed, diff) = apply_recipe(&alignment, &recipe, 4);

        assert_eq!(transformed.length, 12);
        assert_eq!(transformed.length % 3, 0);
        assert!(transformed.sequences.iter().all(|sequence| sequence.starts_with("ATG")));
        assert!(diff
            .final_stop_codons
            .iter()
            .all(|stop| stop.is_terminal && stop.start == 9));
        assert!(diff
            .stop_codons
            .iter()
            .all(|stop| stop.is_terminal && stop.start == 12));
        assert!(diff.pass);
    }

    #[test]
    fn test_pruned_orf_sample_does_not_fail_surviving_locus() {
        let alignment = Alignment::new(
            "exon_orf_pruning".to_string(),
            "exon_orf_pruning.fa".to_string(),
            "/dummy/exon_orf_pruning.fa".to_string(),
            AlignmentFormat::Fasta,
            vec![
                "Taxon_1".to_string(),
                "Taxon_2".to_string(),
                "Taxon_3".to_string(),
                "Taxon_Stop".to_string(),
            ],
            vec![
                "ATGAAAGGG".to_string(),
                "ATGAAAGGG".to_string(),
                "ATGAAAGGG".to_string(),
                "ATGTAAGGG".to_string(),
            ],
        );

        let mut recipe = TrimmingRecipe::default();
        recipe.trim_similarity = false;
        recipe.trim_hmm = false;
        recipe.enable_orf = true;
        recipe.exclude_uce = true;
        recipe.fail_if_no_orf = true;
        recipe.stop_codon_action = StopCodonAction::RemoveSample;
        recipe.trim_external = false;
        recipe.trim_columns = false;
        recipe.trim_coverage = false;
        recipe.min_taxa = 3;
        recipe.min_taxa_occupancy_percent = 0.0;
        recipe.min_length = 0;

        let (transformed, diff) = apply_recipe(&alignment, &recipe, 4);

        assert_eq!(transformed.num_taxa, 3);
        assert_eq!(diff.dropped_taxa, vec!["Taxon_Stop".to_string()]);
        assert!(diff
            .stop_codons
            .iter()
            .any(|stop| stop.taxon == "Taxon_Stop" && stop.start == 3));
        assert!(diff
            .final_stop_codons
            .iter()
            .all(|stop| stop.taxon != "Taxon_Stop"));
        assert!(diff.pass);
        assert!(!diff
            .fail_reasons
            .iter()
            .any(|reason| reason.starts_with("ORF check failed")));
    }

    #[test]
    fn test_reference_guided_orf_separates_flanking_intron() {
        let exon = "ATGGCTGCTGCTGCTGCTGCTGCTGCTTAA";
        let alignment = Alignment::new(
            "exon_reference_test".to_string(),
            "exon_reference_test.fa".to_string(),
            "/dummy/exon_reference_test.fa".to_string(),
            AlignmentFormat::Fasta,
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            vec![
                format!("CCCCCC{exon}GGGGGG"),
                format!("CCCCCC{exon}GGGGGG"),
                format!("CCCCCC{exon}GGGGGG"),
            ],
        );
        let mut recipe = TrimmingRecipe::default();
        recipe.enable_orf = true;
        recipe.exclude_uce = true;
        recipe.orf_use_references = true;
        recipe.orf_search_mode = OrfSearchMode::ReferenceGuided;
        recipe
            .orf_reference_sequences
            .insert(alignment.id.clone(), exon.to_string());
        recipe.trim_similarity = false;
        recipe.trim_hmm = false;
        recipe.trim_external = false;
        recipe.trim_coverage = false;
        recipe.assess_alignment = false;

        let (transformed, diff) = apply_recipe(&alignment, &recipe, 3);

        assert!(diff.orf_reference_evaluated);
        assert!(diff.orf_reference_matched);
        assert_eq!(diff.orf_intron_length, 12);
        assert_eq!(transformed.length, exon.len());
        assert!(transformed.sequences.iter().all(|sequence| sequence == exon));
    }

    #[test]
    fn test_reference_mode_requires_an_exact_name_match_but_skips_uces() {
        let sequences = vec!["ATGGCTGCTTAA".to_string(); 3];
        let coding = Alignment::new(
            "exon_missing_reference".to_string(),
            "exon_missing_reference.fa".to_string(),
            "/dummy/exon_missing_reference.fa".to_string(),
            AlignmentFormat::Fasta,
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            sequences.clone(),
        );
        let uce = Alignment::new(
            "uce-100".to_string(),
            "uce-100.fa".to_string(),
            "/dummy/uce-100.fa".to_string(),
            AlignmentFormat::Fasta,
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            sequences,
        );
        let mut recipe = TrimmingRecipe::default();
        recipe.enable_orf = true;
        recipe.orf_use_references = true;
        recipe.orf_search_mode = OrfSearchMode::ReferenceGuided;
        recipe.trim_similarity = false;
        recipe.trim_hmm = false;
        recipe.trim_external = false;
        recipe.trim_coverage = false;
        recipe.assess_alignment = false;

        let (_, coding_diff) = apply_recipe(&coding, &recipe, 3);
        assert!(!coding_diff.pass);
        assert!(coding_diff.orf_reference_evaluated);
        assert!(coding_diff
            .fail_reasons
            .iter()
            .any(|reason| reason.starts_with("Reference-guided exon failed")));

        let (_, uce_diff) = apply_recipe(&uce, &recipe, 3);
        assert!(uce_diff.pass);
        assert!(!uce_diff.orf_reference_evaluated);
        assert!(!uce_diff.orf_evaluated);
    }

    #[test]
    fn test_reference_candidate_mode_falls_back_when_name_does_not_match() {
        let alignment = Alignment::new(
            "exon_hybrid_missing_reference".to_string(),
            "exon_hybrid_missing_reference.fa".to_string(),
            "/dummy/exon_hybrid_missing_reference.fa".to_string(),
            AlignmentFormat::Fasta,
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            vec!["GCT".repeat(45), "GCC".repeat(45), "GCA".repeat(45)],
        );
        let mut recipe = TrimmingRecipe::default();
        recipe.enable_orf = true;
        recipe.orf_search_mode = OrfSearchMode::ReferenceCandidateOrf;
        recipe.orf_use_references = true;
        recipe.orf_min_coding_score = 0.0;
        recipe.trim_similarity = false;
        recipe.trim_hmm = false;
        recipe.trim_external = false;
        recipe.trim_coverage = false;
        recipe.assess_alignment = false;

        let (_, diff) = apply_recipe(&alignment, &recipe, 3);

        assert!(diff.orf_reference_evaluated);
        assert!(!diff.orf_reference_matched);
        assert!(diff.orf_evaluated);
        assert!(diff.orf_candidate_found);
        assert!(diff.found_valid_orf);
        assert!(diff.pass);
    }

    #[test]
    fn test_reference_candidate_mode_retries_unsliced_alignment_after_guided_failure() {
        let reference = format!("ATG{}TAA", "GCT".repeat(65));
        let sample = format!("ATG{}TAA{}TAA", "GCT".repeat(19), "GCT".repeat(45));
        let alignment = Alignment::new(
            "exon_hybrid_guided_failure".to_string(),
            "exon_hybrid_guided_failure.fa".to_string(),
            "/dummy/exon_hybrid_guided_failure.fa".to_string(),
            AlignmentFormat::Fasta,
            vec!["A".to_string(), "B".to_string(), "C".to_string()],
            vec![sample.clone(), sample.clone(), sample],
        );
        let mut recipe = TrimmingRecipe::default();
        recipe.enable_orf = true;
        recipe.orf_search_mode = OrfSearchMode::ReferenceCandidateOrf;
        recipe.orf_use_references = true;
        recipe.orf_min_coding_score = 0.0;
        recipe
            .orf_reference_sequences
            .insert(alignment.id.clone(), reference);
        recipe.trim_similarity = false;
        recipe.trim_hmm = false;
        recipe.trim_external = false;
        recipe.trim_coverage = false;
        recipe.assess_alignment = false;

        let mut strict_recipe = recipe.clone();
        strict_recipe.orf_search_mode = OrfSearchMode::ReferenceGuided;
        let (_, strict_diff) = apply_recipe(&alignment, &strict_recipe, 3);
        assert!(strict_diff.orf_reference_matched);
        assert!(!strict_diff.found_valid_orf);
        assert!(!strict_diff.pass);

        let (transformed, diff) = apply_recipe(&alignment, &recipe, 3);

        assert!(diff.orf_reference_evaluated);
        assert!(diff.orf_reference_matched);
        assert!(diff.orf_evaluated);
        assert!(diff.orf_candidate_found);
        assert!(diff.found_valid_orf);
        assert!(diff.pass);
        assert!(transformed.length >= 35 * 3);
    }
}
