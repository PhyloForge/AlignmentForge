use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use rayon::prelude::*;
use walkdir::WalkDir;

use crate::algorithms::assess::assess_alignment;
use crate::algorithms::informative::calculate_site_statistics;
use crate::algorithms::stats::{calculate_gap_stats, calculate_gc_percent, compute_mean_divergence};
use crate::models::{Alignment, AlignmentSummary, DatasetOverview, TaxonOccupancy};
use crate::parsers::parse_alignment;
use crate::pipeline::engine::apply_recipe;
use crate::pipeline::recipe::TrimmingRecipe;

/// Builds the runtime recipe for a dataset-wide sample occupancy threshold.
/// A taxon below the threshold is excluded from every processed alignment.
pub fn recipe_with_dataset_sample_filter(
    recipe: &TrimmingRecipe,
    alignments: &[Alignment],
) -> TrimmingRecipe {
    let mut runtime_recipe = recipe.clone();
    runtime_recipe.excluded_taxa.clear();

    let threshold = recipe.min_sample_locus_occupancy_percent;
    if threshold <= 0.0 || alignments.is_empty() {
        return runtime_recipe;
    }

    let mut presence_counts: HashMap<String, usize> = HashMap::new();
    for alignment in alignments {
        let taxa_in_locus: HashSet<&str> = alignment.taxa.iter().map(String::as_str).collect();
        for taxon in taxa_in_locus {
            *presence_counts.entry(taxon.to_string()).or_insert(0) += 1;
        }
    }

    let total_loci = alignments.len() as f64;
    runtime_recipe.excluded_taxa = presence_counts
        .into_iter()
        .filter_map(|(taxon, count)| {
            let occupancy_percent = (count as f64 / total_loci) * 100.0;
            (occupancy_percent < threshold).then_some(taxon)
        })
        .collect();
    runtime_recipe.excluded_taxa.sort();
    runtime_recipe
}

/// Scans a directory of alignments in parallel and returns summaries for all files in a single pass.
pub fn scan_alignment_directory<P: AsRef<Path>, F>(
    dir: P,
    progress_callback: Option<F>,
) -> Result<(Vec<AlignmentSummary>, DatasetOverview, Vec<TaxonOccupancy>, Vec<Alignment>), String>
where
    F: Fn(usize, usize, &str) + Send + Sync + 'static,
{
    let dir_path = dir.as_ref();
    if !dir_path.exists() || !dir_path.is_dir() {
        return Err(format!("Directory does not exist: {:?}", dir_path));
    }

    // Collect all candidate alignment files (skipping hidden and build dirs)
    let file_paths: Vec<PathBuf> = WalkDir::new(dir_path)
        .max_depth(4)
        .into_iter()
        .filter_entry(|e| {
            let name = e.file_name().to_string_lossy();
            !name.starts_with('.') && name != "node_modules" && name != "target" && name != "__MACOSX"
        })
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .map(|e| e.path().to_path_buf())
        .filter(|p| {
            if let Some(ext) = p.extension().and_then(|s| s.to_str()) {
                let ext_lower = ext.to_ascii_lowercase();
                matches!(
                    ext_lower.as_str(),
                    "fa" | "fasta" | "fna" | "faa" | "phy" | "phylip" | "nex" | "nexus" | "aln" | "txt"
                )
            } else {
                false
            }
        })
        .collect();

    let total_files = file_paths.len();
    if total_files == 0 {
        return Err(
            "No supported alignment files (.fa, .fasta, .phy, .nex) found in directory"
                .to_string(),
        );
    }

    // Notify initial count
    if let Some(ref cb) = progress_callback {
        cb(0, total_files, "Starting index...");
    }

    let default_recipe = TrimmingRecipe::default();
    let counter = Arc::new(AtomicUsize::new(0));
    let progress_ref = progress_callback.map(Arc::new);

    // Single-pass ultra-fast parallel extraction: returns (Summary, Vec<(Taxon, non_gap_bp, gap_pct)>, Alignment)
    let processed: Vec<(AlignmentSummary, Vec<(String, usize, f64)>, Alignment)> = file_paths
        .par_iter()
        .filter_map(|path| {
            let parsed = parse_alignment(path).ok().map(|align| {
                let (summary, taxa_info) = fast_index_alignment(&align, &default_recipe, 0);
                (summary, taxa_info, align)
            });

            // Update atomic progress counter (throttled to every 100 files to minimize IPC overhead)
            let cur = counter.fetch_add(1, Ordering::Relaxed) + 1;
            if let Some(ref cb) = progress_ref {
                if cur % 100 == 0 || cur == total_files {
                    let file_name = path
                        .file_name()
                        .map(|s| s.to_string_lossy().to_string())
                        .unwrap_or_default();
                    cb(cur, total_files, &file_name);
                }
            }

            parsed
        })
        .collect();

    if processed.is_empty() {
        return Err("Failed to parse any alignment files in the selected folder".to_string());
    }

    let total_loci = processed.len();
    let mut summaries = Vec::with_capacity(total_loci);
    let mut alignments = Vec::with_capacity(total_loci);
    let mut taxa_stats_map: HashMap<String, (usize, usize, f64)> = HashMap::new();

    for (summary, taxa_info, align) in processed {
        summaries.push(summary);
        alignments.push(align);
        for (taxon, bp, gap_pct) in taxa_info {
            let entry = taxa_stats_map.entry(taxon).or_insert((0, 0, 0.0));
            entry.0 += 1;
            entry.1 += bp;
            entry.2 += gap_pct;
        }
    }

    let total_unique_taxa = taxa_stats_map.len();

    // Re-evaluate quality & occupancy gating with known total dataset taxa
    if total_unique_taxa > 0
        && default_recipe.assess_alignment
        && default_recipe.min_taxa_occupancy_percent > 0.0
    {
        for summary in &mut summaries {
            let (pass, fail_reasons) = assess_alignment(
                summary.num_taxa,
                total_unique_taxa,
                summary.length,
                summary.gap_percent,
                default_recipe.min_taxa,
                default_recipe.min_taxa_occupancy_percent,
                default_recipe.min_length,
                default_recipe.max_gap_percent,
                summary.pis_count,
                summary.pis_percent,
                default_recipe.min_pis_count,
                default_recipe.min_pis_percent,
                summary.variable_count,
                summary.variable_percent,
                default_recipe.min_variable_count,
                default_recipe.min_variable_percent,
            );
            summary.pass = pass;
            summary.fail_reasons = fail_reasons;
        }
    }

    let mut overview = compute_dataset_overview(&summaries);
    overview.total_unique_taxa = total_unique_taxa;

    let mut occupancy: Vec<TaxonOccupancy> = taxa_stats_map
        .into_iter()
        .map(|(taxon_name, (count, total_bp, sum_gap))| TaxonOccupancy {
            taxon_name,
            present_loci_count: count,
            present_loci_percent: (count as f64 / total_loci as f64) * 100.0,
            mean_gap_percent: sum_gap / count.max(1) as f64,
            total_bp,
        })
        .collect();

    occupancy.sort_by(|a, b| b.present_loci_count.cmp(&a.present_loci_count));

    Ok((summaries, overview, occupancy, alignments))
}

/// Ultra-fast single-pass alignment metadata extractor for directory indexing.
/// Processes raw sequences in a single linear pass over byte slices without extra memory allocations.
pub fn fast_index_alignment(
    align: &Alignment,
    recipe: &TrimmingRecipe,
    total_dataset_taxa: usize,
) -> (AlignmentSummary, Vec<(String, usize, f64)>) {
    let num_taxa = align.num_taxa;
    let length = align.length;

    let mut total_bp = 0usize;
    let mut gap_count = 0usize;
    let mut gc_count = 0usize;
    let mut total_valid_bases = 0usize;

    let mut taxa_stats = Vec::with_capacity(num_taxa);

    for (taxon, seq) in align.taxa.iter().zip(align.sequences.iter()) {
        let seq_len = seq.len();
        let mut sample_gaps = 0usize;

        for &b in seq.as_bytes() {
            match b {
                b'-' | b'?' => {
                    sample_gaps += 1;
                    gap_count += 1;
                }
                b'N' | b'n' => {
                    if recipe.count_n_as_gap {
                        sample_gaps += 1;
                        gap_count += 1;
                    }
                }
                b'G' | b'g' | b'C' | b'c' => {
                    gc_count += 1;
                    total_valid_bases += 1;
                    total_bp += 1;
                }
                b'A' | b'a' | b'T' | b't' | b'U' | b'u' => {
                    total_valid_bases += 1;
                    total_bp += 1;
                }
                _ => {
                    total_valid_bases += 1;
                    total_bp += 1;
                }
            }
        }

        let non_gap_bp = seq_len.saturating_sub(sample_gaps);
        let sample_gap_pct = if seq_len > 0 {
            (sample_gaps as f64 / seq_len as f64) * 100.0
        } else {
            100.0
        };
        taxa_stats.push((taxon.clone(), non_gap_bp, sample_gap_pct));
    }

    let total_matrix_cells = num_taxa * length;
    let gap_percent = if total_matrix_cells > 0 {
        (gap_count as f64 / total_matrix_cells as f64) * 100.0
    } else {
        0.0
    };

    let gc_percent = if total_valid_bases > 0 {
        (gc_count as f64 / total_valid_bases as f64) * 100.0
    } else {
        0.0
    };

    let site_stats = calculate_site_statistics(&align.sequences, true);
    let variable_count = site_stats.variable_count;
    let variable_percent = site_stats.variable_percent;
    let pis_count = site_stats.pis_count;
    let pis_percent = site_stats.pis_percent;

    let mean_divergence = compute_mean_divergence(&align.sequences);

    let (pass, fail_reasons) = if recipe.assess_alignment {
        assess_alignment(
            num_taxa,
            total_dataset_taxa,
            length,
            gap_percent,
            recipe.min_taxa,
            recipe.min_taxa_occupancy_percent,
            recipe.min_length,
            recipe.max_gap_percent,
            pis_count,
            pis_percent,
            recipe.min_pis_count,
            recipe.min_pis_percent,
            variable_count,
            variable_percent,
            recipe.min_variable_count,
            recipe.min_variable_percent,
        )
    } else {
        (true, Vec::new())
    };

    let summary = AlignmentSummary {
        id: align.id.clone(),
        file_name: align.file_name.clone(),
        file_path: align.file_path.clone(),
        format: align.format,
        num_taxa,
        length,
        total_basepairs: total_bp,
        gap_count,
        gap_percent,
        variable_count,
        variable_percent,
        pis_count,
        pis_percent,
        mean_divergence,
        gc_percent,
        pass,
        fail_reasons,
        orf_valid: true,
        orf_evaluated: false,
        orf_candidate_found: false,
        orf_frame: None,
        orf_start: None,
        orf_end: None,
        orf_support_count: 0,
        orf_support_percent: 0.0,
        orf_retained_samples: 0,
        orf_candidate_length_aa: 0,
        orf_coding_score: 0.0,
        orf_amino_acid_conservation: 0.0,
        orf_frame_contrast: 0.0,
        orf_reference_evaluated: false,
        orf_reference_matched: false,
        orf_reference_identity: 0.0,
        orf_reference_coverage: 0.0,
        orf_intron_length: 0,
        raw_num_taxa: num_taxa,
        raw_length: length,
        raw_gap_percent: gap_percent,
        retained_taxa: align.taxa.clone(),
        retained_taxon_basepairs: taxa_stats
            .iter()
            .map(|(taxon, basepairs, _)| (taxon.clone(), *basepairs))
            .collect(),
    };

    (summary, taxa_stats)
}

pub fn summarize_alignment(
    alignment: &Alignment,
    recipe: &TrimmingRecipe,
    total_dataset_taxa: usize,
) -> AlignmentSummary {
    let raw_length = alignment.length;
    let raw_num_taxa = alignment.num_taxa;

    // Apply trimming recipe to compute post-trimming metrics & quality assessment
    let (trimmed, diff) = apply_recipe(alignment, recipe, total_dataset_taxa);

    let (gap_count, total_chars, _) = calculate_gap_stats(&trimmed.sequences);
    let total_bp = total_chars.saturating_sub(gap_count);
    let gap_percent = diff.new_gap_percent;
    let variable_count = diff.new_variable;
    let variable_percent = if trimmed.length > 0 {
        (variable_count as f64 / trimmed.length as f64) * 100.0
    } else {
        0.0
    };
    let pis_count = diff.new_pis;
    let pis_percent = if trimmed.length > 0 {
        (pis_count as f64 / trimmed.length as f64) * 100.0
    } else {
        0.0
    };
    let mean_divergence = compute_mean_divergence(&trimmed.sequences);
    let gc_percent = calculate_gc_percent(&trimmed.sequences);

    let mut fail_reasons = diff.fail_reasons.clone();

    if trimmed.num_taxa == 0 {
        let r0 = "0 surviving taxa (all samples pruned)".to_string();
        if !fail_reasons.contains(&r0) {
            fail_reasons.push(r0);
        }
    }

    if recipe.assess_alignment {
        let (_, assess_reasons) = assess_alignment(
            trimmed.num_taxa,
            if total_dataset_taxa > 0 { total_dataset_taxa } else { raw_num_taxa },
            trimmed.length,
            gap_percent,
            recipe.min_taxa,
            recipe.min_taxa_occupancy_percent,
            recipe.min_length,
            recipe.max_gap_percent,
            pis_count,
            pis_percent,
            recipe.min_pis_count,
            recipe.min_pis_percent,
            variable_count,
            variable_percent,
            recipe.min_variable_count,
            recipe.min_variable_percent,
        );
        for r in assess_reasons {
            if !fail_reasons.contains(&r) {
                fail_reasons.push(r);
            }
        }
    }

    let pass = diff.pass && fail_reasons.is_empty() && trimmed.num_taxa > 0;

    AlignmentSummary {
        id: alignment.id.clone(),
        file_name: alignment.file_name.clone(),
        file_path: alignment.file_path.clone(),
        format: alignment.format,
        num_taxa: trimmed.num_taxa,
        length: trimmed.length,
        total_basepairs: total_bp,
        gap_count,
        gap_percent,
        variable_count,
        variable_percent,
        pis_count,
        pis_percent,
        mean_divergence,
        gc_percent,
        pass,
        fail_reasons,
        orf_valid: diff.found_valid_orf,
        orf_evaluated: diff.orf_evaluated,
        orf_candidate_found: diff.orf_candidate_found,
        orf_frame: diff.orf_frame,
        orf_start: diff.orf_start,
        orf_end: diff.orf_end,
        orf_support_count: diff.orf_support_count,
        orf_support_percent: diff.orf_support_percent,
        orf_retained_samples: diff.orf_retained_samples,
        orf_candidate_length_aa: diff.orf_candidate_length_aa,
        orf_coding_score: diff.orf_coding_score,
        orf_amino_acid_conservation: diff.orf_amino_acid_conservation,
        orf_frame_contrast: diff.orf_frame_contrast,
        orf_reference_evaluated: diff.orf_reference_evaluated,
        orf_reference_matched: diff.orf_reference_matched,
        orf_reference_identity: diff.orf_reference_identity,
        orf_reference_coverage: diff.orf_reference_coverage,
        orf_intron_length: diff.orf_intron_length,
        raw_num_taxa,
        raw_length,
        raw_gap_percent: diff.old_gap_percent,
        retained_taxa: trimmed.taxa.clone(),
        retained_taxon_basepairs: trimmed
            .taxa
            .iter()
            .zip(trimmed.sequences.iter())
            .map(|(taxon, sequence)| {
                let basepairs = sequence
                    .bytes()
                    .filter(|state| !matches!(state, b'-' | b'?' | b'N' | b'n'))
                    .count();
                (taxon.clone(), basepairs)
            })
            .collect(),
    }
}

pub fn compute_dataset_overview(summaries: &[AlignmentSummary]) -> DatasetOverview {
    let total_alignments = summaries.len();
    if total_alignments == 0 {
        return DatasetOverview {
            total_alignments: 0,
            passed_alignments: 0,
            discarded_alignments: 0,
            total_unique_taxa: 0,
            mean_taxa: 0.0,
            mean_length: 0.0,
            mean_gap_percent: 0.0,
            mean_pis: 0.0,
            total_matrix_basepairs: 0,
        };
    }

    let passed_alignments = summaries.iter().filter(|s| s.pass).count();
    let discarded_alignments = total_alignments - passed_alignments;

    let sum_taxa: usize = summaries.iter().map(|s| s.num_taxa).sum();
    let sum_length: usize = summaries.iter().map(|s| s.length).sum();
    let sum_gap_percent: f64 = summaries.iter().map(|s| s.gap_percent).sum();
    let sum_pis: usize = summaries.iter().map(|s| s.pis_count).sum();
    let total_matrix_basepairs: usize = summaries.iter().map(|s| s.total_basepairs).sum();

    DatasetOverview {
        total_alignments,
        passed_alignments,
        discarded_alignments,
        total_unique_taxa: 0, // Assigned in caller
        mean_taxa: sum_taxa as f64 / total_alignments as f64,
        mean_length: sum_length as f64 / total_alignments as f64,
        mean_gap_percent: sum_gap_percent / total_alignments as f64,
        mean_pis: sum_pis as f64 / total_alignments as f64,
        total_matrix_basepairs,
    }
}

/// Evaluates a TrimmingRecipe across pre-parsed in-memory Alignment objects in parallel
pub fn evaluate_recipe_on_alignments(
    alignments: &[Alignment],
    recipe: &TrimmingRecipe,
    total_unique_taxa: usize,
) -> (Vec<AlignmentSummary>, DatasetOverview) {
    evaluate_recipe_on_alignments_with_progress(
        alignments,
        recipe,
        total_unique_taxa,
        None::<fn(usize, usize, &str)>,
    )
}

/// Evaluates a recipe across cached alignments and reports completed loci as workers finish.
pub fn evaluate_recipe_on_alignments_with_progress<F>(
    alignments: &[Alignment],
    recipe: &TrimmingRecipe,
    total_unique_taxa: usize,
    progress_callback: Option<F>,
) -> (Vec<AlignmentSummary>, DatasetOverview)
where
    F: Fn(usize, usize, &str) + Send + Sync,
{
    evaluate_recipe_on_alignments_with_progress_and_cancel(
        alignments,
        recipe,
        total_unique_taxa,
        progress_callback,
        || false,
    )
    .expect("a non-cancellable catalog evaluation cannot be cancelled")
}

/// Evaluates cached alignments while allowing a newer UI request to stop stale work.
pub fn evaluate_recipe_on_alignments_with_progress_and_cancel<F, C>(
    alignments: &[Alignment],
    recipe: &TrimmingRecipe,
    total_unique_taxa: usize,
    progress_callback: Option<F>,
    should_cancel: C,
) -> Option<(Vec<AlignmentSummary>, DatasetOverview)>
where
    F: Fn(usize, usize, &str) + Send + Sync,
    C: Fn() -> bool + Send + Sync,
{
    let total = alignments.len();
    let completed = AtomicUsize::new(0);
    let runtime_recipe = recipe_with_dataset_sample_filter(recipe, alignments);
    let summaries: Vec<AlignmentSummary> = alignments
        .par_iter()
        .filter_map(|align| {
            if should_cancel() {
                return None;
            }
            let summary = summarize_alignment(align, &runtime_recipe, total_unique_taxa);
            if should_cancel() {
                return None;
            }
            let current = completed.fetch_add(1, Ordering::Relaxed) + 1;
            if let Some(ref callback) = progress_callback {
                callback(current, total, &align.file_name);
            }
            Some(summary)
        })
        .collect();

    if should_cancel() {
        return None;
    }
    let overview = compute_dataset_overview(&summaries);
    Some((summaries, overview))
}

/// Evaluates a TrimmingRecipe across all summary entries by re-parsing from disk in parallel (fallback)
pub fn evaluate_recipe_on_summaries(
    paths: &[String],
    recipe: &TrimmingRecipe,
    total_unique_taxa: usize,
) -> (Vec<AlignmentSummary>, DatasetOverview) {
    evaluate_recipe_on_summaries_with_progress(
        paths,
        recipe,
        total_unique_taxa,
        None::<fn(usize, usize, &str)>,
    )
}

/// Disk-backed fallback for catalog recalculation with completed-locus progress.
pub fn evaluate_recipe_on_summaries_with_progress<F>(
    paths: &[String],
    recipe: &TrimmingRecipe,
    total_unique_taxa: usize,
    progress_callback: Option<F>,
) -> (Vec<AlignmentSummary>, DatasetOverview)
where
    F: Fn(usize, usize, &str) + Send + Sync,
{
    let alignments: Vec<Alignment> = paths
        .par_iter()
        .filter_map(|path| parse_alignment(path).ok())
        .collect();
    evaluate_recipe_on_alignments_with_progress(
        &alignments,
        recipe,
        total_unique_taxa,
        progress_callback,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithms::orf::StopCodonAction;
    use crate::models::AlignmentFormat;

    #[test]
    fn test_scan_directory_test_data() {
        let test_dir = std::path::Path::new("../test_data");
        if test_dir.exists() {
            let (summaries, overview, occupancy, _alignments) =
                scan_alignment_directory(test_dir, None::<fn(usize, usize, &str)>).unwrap();
            assert!(!summaries.is_empty());
            assert!(overview.total_alignments >= 4);
            assert!(!occupancy.is_empty());
        }
    }

    #[test]
    fn test_catalog_recalculation_applies_orf_before_assessment() {
        let alignment = Alignment::new(
            "exon_001".to_string(),
            "exon_001.fa".to_string(),
            "/dummy/exon_001.fa".to_string(),
            AlignmentFormat::Fasta,
            vec![
                "Taxon_1".to_string(),
                "Taxon_2".to_string(),
                "Taxon_3".to_string(),
                "Taxon_Pseudogene".to_string(),
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
        recipe.exclude_uce = false;
        recipe.fail_if_no_orf = true;
        recipe.stop_codon_action = StopCodonAction::RemoveSample;
        recipe.trim_external = false;
        recipe.trim_columns = false;
        recipe.trim_coverage = false;
        recipe.min_taxa = 4;
        recipe.min_taxa_occupancy_percent = 10.0;
        recipe.min_length = 0;

        let progress_calls = AtomicUsize::new(0);
        let (summaries, overview) = evaluate_recipe_on_alignments_with_progress(
            &[alignment],
            &recipe,
            4,
            Some(|current, total, _file_name: &str| {
                assert_eq!(current, 1);
                assert_eq!(total, 1);
                progress_calls.fetch_add(1, Ordering::Relaxed);
            }),
        );

        assert_eq!(summaries.len(), 1);
        assert_eq!(progress_calls.load(Ordering::Relaxed), 1);
        assert_eq!(summaries[0].num_taxa, 3);
        assert!(!summaries[0].pass);
        assert!(summaries[0]
            .fail_reasons
            .iter()
            .any(|reason| reason.contains("Taxa count (3 < min 4)")));
        assert!(!summaries[0]
            .fail_reasons
            .iter()
            .any(|reason| reason.starts_with("ORF check failed")));
        assert_eq!(overview.discarded_alignments, 1);
    }

    #[test]
    fn test_dataset_sample_occupancy_drops_taxon_from_every_locus() {
        let alignments = vec![
            Alignment::new(
                "locus_1".to_string(),
                "locus_1.fa".to_string(),
                "/dummy/locus_1.fa".to_string(),
                AlignmentFormat::Fasta,
                vec!["Common".to_string(), "Rare".to_string()],
                vec!["ACGT".to_string(), "ACGT".to_string()],
            ),
            Alignment::new(
                "locus_2".to_string(),
                "locus_2.fa".to_string(),
                "/dummy/locus_2.fa".to_string(),
                AlignmentFormat::Fasta,
                vec!["Common".to_string()],
                vec!["ACGT".to_string()],
            ),
        ];

        let mut recipe = TrimmingRecipe::default();
        recipe.trim_similarity = false;
        recipe.trim_external = false;
        recipe.trim_coverage = false;
        recipe.assess_alignment = false;
        recipe.min_sample_locus_occupancy_percent = 60.0;

        let runtime_recipe = recipe_with_dataset_sample_filter(&recipe, &alignments);
        assert_eq!(runtime_recipe.excluded_taxa, vec!["Rare".to_string()]);

        let (summaries, _) = evaluate_recipe_on_alignments(&alignments, &recipe, 2);
        assert_eq!(summaries[0].num_taxa, 1);
        assert_eq!(summaries[1].num_taxa, 1);
    }

    #[test]
    fn test_superseded_catalog_evaluation_cancels_before_processing() {
        let alignment = Alignment::new(
            "locus_1".to_string(),
            "locus_1.fa".to_string(),
            "/dummy/locus_1.fa".to_string(),
            AlignmentFormat::Fasta,
            vec!["Taxon_1".to_string()],
            vec!["ACGT".to_string()],
        );

        let result = evaluate_recipe_on_alignments_with_progress_and_cancel(
            &[alignment],
            &TrimmingRecipe::default(),
            1,
            None::<fn(usize, usize, &str)>,
            || true,
        );

        assert!(result.is_none());
    }
}
