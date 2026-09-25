use std::borrow::Borrow;
use std::collections::{HashMap, HashSet};
#[cfg(not(target_arch = "wasm32"))]
use std::path::{Path, PathBuf};
#[cfg(not(target_arch = "wasm32"))]
use std::sync::atomic::{AtomicUsize, Ordering};
#[cfg(not(target_arch = "wasm32"))]
use std::sync::Arc;
#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;
#[cfg(not(target_arch = "wasm32"))]
use walkdir::WalkDir;

use crate::algorithms::assess::assess_alignment;
use crate::algorithms::orf::should_skip_orf_locus;
use crate::algorithms::stats::{calculate_gap_stats, calculate_gc_percent, compute_mean_divergence};
use crate::models::{Alignment, AlignmentSummary, DatasetOverview};
#[cfg(not(target_arch = "wasm32"))]
use crate::models::TaxonOccupancy;
#[cfg(not(target_arch = "wasm32"))]
use crate::parsers::parse_alignment;
use crate::pipeline::engine::apply_recipe;
use crate::pipeline::recipe::TrimmingRecipe;

/// Builds the runtime recipe for a dataset-wide sample occupancy threshold.
/// A taxon below the threshold is excluded from every processed alignment.
pub fn recipe_with_dataset_sample_filter<A: Borrow<Alignment>>(
    recipe: &TrimmingRecipe,
    alignments: &[A],
) -> TrimmingRecipe {
    // Borrow the taxon lists rather than copying the alignments.
    let presence: Vec<&[String]> = alignments
        .iter()
        .map(|alignment| alignment.borrow().taxa.as_slice())
        .collect();
    recipe_with_taxon_presence(recipe, &presence)
}

/// Same filter, driven by taxon lists alone. Callers that only hold sample
/// names avoid cloning every sequence in the dataset.
pub fn recipe_with_taxon_presence<T: AsRef<[String]>>(
    recipe: &TrimmingRecipe,
    loci_taxa: &[T],
) -> TrimmingRecipe {
    let mut runtime_recipe = recipe.clone();
    // Start from the samples the user discarded by hand, then add whatever the
    // occupancy threshold excludes. Clearing the list outright would throw the
    // manual choices away on every recalculation.
    runtime_recipe.excluded_taxa = recipe.discarded_taxa.clone();

    let threshold = recipe.min_sample_locus_occupancy_percent;
    if threshold <= 0.0 || loci_taxa.is_empty() {
        runtime_recipe.excluded_taxa.sort();
        runtime_recipe.excluded_taxa.dedup();
        return runtime_recipe;
    }

    let mut presence_counts: HashMap<&str, usize> = HashMap::new();
    for taxa in loci_taxa {
        let taxa_in_locus: HashSet<&str> = taxa.as_ref().iter().map(String::as_str).collect();
        for taxon in taxa_in_locus {
            *presence_counts.entry(taxon).or_insert(0) += 1;
        }
    }

    let total_loci = loci_taxa.len() as f64;
    runtime_recipe
        .excluded_taxa
        .extend(presence_counts.into_iter().filter_map(|(taxon, count)| {
            let occupancy_percent = (count as f64 / total_loci) * 100.0;
            (occupancy_percent < threshold).then_some(taxon.to_string())
        }));
    runtime_recipe.excluded_taxa.sort();
    runtime_recipe.excluded_taxa.dedup();
    runtime_recipe
}

/// Builds the ordinary alignment-processing branch used by Catalog QC and
/// general alignment export. ORF extraction is deliberately independent and
/// must never remove samples or columns before Catalog pass/fail is assessed.
pub fn recipe_without_orf_analysis(recipe: &TrimmingRecipe) -> TrimmingRecipe {
    let mut catalog_recipe = recipe.clone();
    catalog_recipe.enable_orf = false;
    catalog_recipe.orf_use_references = false;
    // Nothing reads the references without ORF analysis. Dropping them keeps
    // later copies of this recipe small.
    catalog_recipe.orf_reference_sequences = HashMap::new();
    catalog_recipe
}

/// The two recipes that measure each locus: the ordinary branch and the ORF
/// branch. A run builds them once, because a recipe can hold thousands of
/// reference sequences.
pub struct CatalogRecipes {
    catalog: TrimmingRecipe,
    orf: TrimmingRecipe,
}

impl CatalogRecipes {
    pub fn new(recipe: &TrimmingRecipe) -> Self {
        // Gating is applied later, so both branches run with it switched off and
        // record only the failures that processing itself produced.
        let mut orf = recipe.clone();
        orf.assess_alignment = false;
        let catalog = recipe_without_orf_analysis(&orf);
        Self { catalog, orf }
    }
}

/// How many directory levels below the chosen folder are searched. The browser
/// applies this limit and the skipped names in `isScannedPath` (src/tauriClient.ts).
pub const MAX_SCAN_DEPTH: usize = 4;

/// Extensions the scanner accepts. The browser file picker uses the same list
/// in `SUPPORTED_ALIGNMENT_EXTENSIONS` (src/tauriClient.ts); keep them in step.
pub const SUPPORTED_ALIGNMENT_EXTENSIONS: [&str; 10] = [
    "fa", "fasta", "fna", "ffn", "phy", "phylip", "nex", "nexus", "aln", "txt",
];

/// One file that could not be read, reported to the user instead of dropped.
#[cfg(not(target_arch = "wasm32"))]
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ParseFailure {
    pub file_name: String,
    pub file_path: String,
    pub error: String,
}

#[cfg(not(target_arch = "wasm32"))]
pub type ScanOutcome = (Vec<Alignment>, Vec<TaxonOccupancy>, Vec<ParseFailure>);

/// Reads a directory of alignments in parallel. Returns the alignments, the
/// dataset-wide sample occupancy, and the files that could not be parsed.
#[cfg(not(target_arch = "wasm32"))]
pub fn scan_alignment_directory<P: AsRef<Path>, F>(
    dir: P,
    progress_callback: Option<F>,
) -> Result<ScanOutcome, String>
where
    F: Fn(usize, usize, &str) + Send + Sync + 'static,
{
    let dir_path = dir.as_ref();
    if !dir_path.exists() || !dir_path.is_dir() {
        return Err(format!("Directory does not exist: {:?}", dir_path));
    }

    // Collect all candidate alignment files (skipping hidden and build dirs).
    // The depth limit keeps an accidental pick of a large tree, such as a home
    // directory, from walking the whole disk.
    let file_paths: Vec<PathBuf> = WalkDir::new(dir_path)
        .max_depth(MAX_SCAN_DEPTH)
        .into_iter()
        .filter_entry(|e| {
            // The chosen folder itself is always read, as in the browser build.
            let name = e.file_name().to_string_lossy();
            e.depth() == 0
                || (!name.starts_with('.') && name != "node_modules" && name != "target" && name != "__MACOSX")
        })
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .map(|e| e.path().to_path_buf())
        .filter(|p| {
            if let Some(ext) = p.extension().and_then(|s| s.to_str()) {
                let ext_lower = ext.to_ascii_lowercase();
                SUPPORTED_ALIGNMENT_EXTENSIONS.contains(&ext_lower.as_str())
            } else {
                false
            }
        })
        .collect();

    let total_files = file_paths.len();
    if total_files == 0 {
        return Err(format!(
            "No supported alignment files found in directory. AlignmentForge reads: {}",
            SUPPORTED_ALIGNMENT_EXTENSIONS
                .iter()
                .map(|ext| format!(".{ext}"))
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }

    // Notify initial count
    if let Some(ref cb) = progress_callback {
        cb(0, total_files, "Starting index...");
    }

    let counter = Arc::new(AtomicUsize::new(0));
    let progress_ref = progress_callback.map(Arc::new);

    type IndexedAlignment = (Vec<(String, usize, f64)>, Alignment);
    let processed: Vec<Result<IndexedAlignment, ParseFailure>> = file_paths
        .par_iter()
        .map(|path| {
            let parsed = match parse_alignment(path) {
                Ok(align) => Ok((sample_coverage(&align), align)),
                Err(error) => Err(ParseFailure {
                    file_name: path
                        .file_name()
                        .map(|s| s.to_string_lossy().to_string())
                        .unwrap_or_default(),
                    file_path: path.to_string_lossy().to_string(),
                    error,
                }),
            };

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

    let mut parse_failures = Vec::new();
    let mut succeeded = Vec::with_capacity(processed.len());
    for outcome in processed {
        match outcome {
            Ok(indexed) => succeeded.push(indexed),
            Err(failure) => parse_failures.push(failure),
        }
    }
    parse_failures.sort_by(|a, b| a.file_name.cmp(&b.file_name));

    if succeeded.is_empty() {
        let detail = parse_failures
            .iter()
            .take(3)
            .map(|failure| format!("{}: {}", failure.file_name, failure.error))
            .collect::<Vec<_>>()
            .join("; ");
        return Err(format!(
            "Failed to parse any alignment files in the selected folder. {detail}"
        ));
    }

    let processed = succeeded;
    let total_loci = processed.len();
    let mut alignments = Vec::with_capacity(total_loci);
    let mut taxa_stats_map: HashMap<String, (usize, usize, f64)> = HashMap::new();

    for (taxa_info, align) in processed {
        alignments.push(align);
        for (taxon, bp, gap_pct) in taxa_info {
            let entry = taxa_stats_map.entry(taxon).or_insert((0, 0, 0.0));
            entry.0 += 1;
            entry.1 += bp;
            entry.2 += gap_pct;
        }
    }

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

    Ok((alignments, occupancy, parse_failures))
}

/// Base pairs and missing-data percent of each sample in an unprocessed
/// alignment. `-`, `?`, and `N` count as missing data.
#[cfg(not(target_arch = "wasm32"))]
fn sample_coverage(align: &Alignment) -> Vec<(String, usize, f64)> {
    align
        .taxa
        .iter()
        .zip(align.sequences.iter())
        .map(|(taxon, sequence)| {
            let basepairs = count_observed_bases(sequence);
            let missing = sequence.len() - basepairs;
            let missing_percent = if sequence.is_empty() {
                100.0
            } else {
                (missing as f64 / sequence.len() as f64) * 100.0
            };
            (taxon.clone(), basepairs, missing_percent)
        })
        .collect()
}

/// Counts observed bases in a sequence. `-`, `?`, and `N` are missing data and
/// never count; every other state does. This is the single definition of a
/// "base pair" used by both the directory index and the post-trimming summary.
pub fn count_observed_bases(sequence: &str) -> usize {
    sequence
        .bytes()
        .filter(|state| !matches!(state, b'-' | b'?' | b'N' | b'n'))
        .count()
}

fn taxon_basepair_map(taxa: &[String], sequences: &[String]) -> HashMap<String, usize> {
    taxa.iter()
        .zip(sequences.iter())
        .map(|(taxon, sequence)| (taxon.clone(), count_observed_bases(sequence)))
        .collect()
}

/// A locus summary before the pass and fail thresholds are applied.
///
/// Everything expensive lives here: the alignment has been trimmed and every
/// metric measured. Only the gating decision is outstanding, so changing a
/// threshold can reuse this instead of running the pipeline again.
#[derive(Debug, Clone)]
pub struct UnassessedSummary {
    /// `fail_reasons` holds pipeline failures only, not threshold failures.
    /// The per-sample lists are empty; `retention` holds them.
    pub summary: AlignmentSummary,
    pub retention: SampleRetention,
    pub pipeline_pass: bool,
}

/// The samples that a processed locus keeps, with their base pairs. Only the QC
/// occupancy chart and the Matrix view read these lists, so the summaries sent
/// on each recipe change leave them out.
#[derive(Debug, Clone, serde::Serialize)]
pub struct SampleRetention {
    pub file_path: String,
    pub retained_taxa: Vec<String>,
    pub retained_taxon_basepairs: HashMap<String, usize>,
    pub orf_retained_taxa: Vec<String>,
    pub orf_retained_taxon_basepairs: HashMap<String, usize>,
}

/// Applies the pass and fail thresholds to an already measured locus. The
/// summary leaves out the per-sample lists.
pub fn apply_assessment(
    unassessed: &UnassessedSummary,
    recipe: &TrimmingRecipe,
    total_dataset_taxa: usize,
) -> AlignmentSummary {
    let mut summary = unassessed.summary.clone();
    let mut fail_reasons = summary.fail_reasons.clone();

    if summary.num_taxa == 0 {
        let reason = "0 surviving taxa (all samples pruned)".to_string();
        if !fail_reasons.contains(&reason) {
            fail_reasons.push(reason);
        }
    }
    if summary.length == 0 {
        let reason = "0 surviving columns (all alignment sites removed)".to_string();
        if !fail_reasons.contains(&reason) {
            fail_reasons.push(reason);
        }
    }

    if recipe.assess_alignment {
        let (_, assess_reasons) = assess_alignment(
            summary.num_taxa,
            if total_dataset_taxa > 0 {
                total_dataset_taxa
            } else {
                summary.raw_num_taxa
            },
            summary.length,
            summary.gap_percent,
            recipe.min_taxa,
            recipe.min_taxa_occupancy_percent,
            recipe.min_length,
            recipe.max_gap_percent,
            summary.pis_count,
            summary.pis_percent,
            recipe.min_pis_count,
            recipe.min_pis_percent,
            summary.variable_count,
            summary.variable_percent,
            recipe.min_variable_count,
            recipe.min_variable_percent,
        );
        for reason in assess_reasons {
            if !fail_reasons.contains(&reason) {
                fail_reasons.push(reason);
            }
        }
    }

    summary.pass =
        unassessed.pipeline_pass
            && fail_reasons.is_empty()
            && summary.num_taxa > 0
            && summary.length > 0;
    summary.fail_reasons = fail_reasons;
    summary
}

/// Adds the per-sample lists to a summary.
pub fn with_retention(
    mut summary: AlignmentSummary,
    retention: SampleRetention,
) -> AlignmentSummary {
    summary.retained_taxa = retention.retained_taxa;
    summary.retained_taxon_basepairs = retention.retained_taxon_basepairs;
    summary.orf_retained_taxa = retention.orf_retained_taxa;
    summary.orf_retained_taxon_basepairs = retention.orf_retained_taxon_basepairs;
    summary
}

/// Runs the pipeline and measures the locus, without applying the thresholds.
pub fn summarize_alignment_unassessed(
    alignment: &Alignment,
    recipes: &CatalogRecipes,
    total_dataset_taxa: usize,
) -> UnassessedSummary {
    let raw_length = alignment.length;
    let raw_num_taxa = alignment.num_taxa;

    // Catalog QC and ORF analysis are parallel outcomes. Catalog metrics come
    // from the ordinary alignment branch; ORF metadata comes from its own
    // branch and cannot change Catalog pass/fail or its filter measurements.
    let (trimmed, diff) = apply_recipe(alignment, &recipes.catalog, total_dataset_taxa);

    // The ORF branch only differs when ORF analysis actually runs on this locus.
    // A locus that ORF skips, such as an excluded UCE, would repeat the whole
    // pipeline for an identical result.
    let orf_recipe = &recipes.orf;
    let orf_branch_differs = orf_recipe.enable_orf
        && !(orf_recipe.exclude_uce
            && should_skip_orf_locus(&alignment.id, orf_recipe.orf_search_mode));
    let orf_run =
        orf_branch_differs.then(|| apply_recipe(alignment, orf_recipe, total_dataset_taxa));
    let (orf_alignment, orf_diff) = match &orf_run {
        Some((orf_alignment, orf_diff)) => (orf_alignment, orf_diff),
        None => (&trimmed, &diff),
    };

    let (gap_count, _, _) = calculate_gap_stats(&trimmed.sequences);
    let total_bp: usize = trimmed
        .sequences
        .iter()
        .map(|sequence| count_observed_bases(sequence))
        .sum();
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

    let fail_reasons = diff.fail_reasons.clone();
    let pipeline_pass = diff.pass;
    let pass = pipeline_pass;

    let summary = AlignmentSummary {
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
        orf_valid: orf_diff.found_valid_orf,
        orf_evaluated: orf_diff.orf_evaluated,
        orf_candidate_found: orf_diff.orf_candidate_found,
        orf_frame: orf_diff.orf_frame,
        orf_start: orf_diff.orf_start,
        orf_end: orf_diff.orf_end,
        orf_support_count: orf_diff.orf_support_count,
        orf_support_percent: orf_diff.orf_support_percent,
        orf_retained_samples: orf_diff.orf_retained_samples,
        orf_candidate_length_aa: orf_diff.orf_candidate_length_aa,
        orf_coding_score: orf_diff.orf_coding_score,
        orf_amino_acid_conservation: orf_diff.orf_amino_acid_conservation,
        orf_frame_contrast: orf_diff.orf_frame_contrast,
        orf_reference_evaluated: orf_diff.orf_reference_evaluated,
        orf_reference_matched: orf_diff.orf_reference_matched,
        orf_reference_identity: orf_diff.orf_reference_identity,
        orf_reference_coverage: orf_diff.orf_reference_coverage,
        orf_intron_length: orf_diff.orf_intron_length,
        raw_num_taxa,
        raw_length,
        raw_gap_percent: diff.old_gap_percent,
        retained_taxa: Vec::new(),
        retained_taxon_basepairs: HashMap::new(),
        orf_retained_taxa: Vec::new(),
        orf_retained_taxon_basepairs: HashMap::new(),
    };
    let retention = SampleRetention {
        file_path: alignment.file_path.clone(),
        retained_taxa: trimmed.taxa.clone(),
        retained_taxon_basepairs: taxon_basepair_map(&trimmed.taxa, &trimmed.sequences),
        orf_retained_taxa: orf_alignment.taxa.clone(),
        orf_retained_taxon_basepairs: taxon_basepair_map(
            &orf_alignment.taxa,
            &orf_alignment.sequences,
        ),
    };

    UnassessedSummary {
        summary,
        retention,
        pipeline_pass,
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

/// Measures every locus in parallel, without applying the thresholds, and
/// reports each completed locus. Returns `None` when `should_cancel` stops the
/// run.
#[cfg(not(target_arch = "wasm32"))]
pub fn measure_alignments<A, F, C>(
    alignments: &[A],
    recipe: &TrimmingRecipe,
    total_unique_taxa: usize,
    progress_callback: Option<F>,
    should_cancel: C,
) -> Option<Vec<UnassessedSummary>>
where
    A: Borrow<Alignment> + Sync,
    F: Fn(usize, usize, &str) + Send + Sync,
    C: Fn() -> bool + Send + Sync,
{
    let total = alignments.len();
    let completed = AtomicUsize::new(0);
    let runtime_recipe = recipe_with_dataset_sample_filter(recipe, alignments);
    let recipes = CatalogRecipes::new(&runtime_recipe);
    let measured: Vec<UnassessedSummary> = alignments
        .par_iter()
        .filter_map(|alignment| {
            let alignment = alignment.borrow();
            if should_cancel() {
                return None;
            }
            let unassessed = summarize_alignment_unassessed(alignment, &recipes, total_unique_taxa);
            if should_cancel() {
                return None;
            }
            let current = completed.fetch_add(1, Ordering::Relaxed) + 1;
            if let Some(ref callback) = progress_callback {
                callback(current, total, &alignment.file_name);
            }
            Some(unassessed)
        })
        .collect();

    if should_cancel() {
        return None;
    }
    Some(measured)
}

/// Applies the thresholds to measured loci. The summaries leave out the
/// per-sample lists.
pub fn assess_measured(
    measured: &[UnassessedSummary],
    recipe: &TrimmingRecipe,
    total_unique_taxa: usize,
) -> (Vec<AlignmentSummary>, DatasetOverview) {
    let summaries: Vec<AlignmentSummary> = measured
        .iter()
        .map(|locus| apply_assessment(locus, recipe, total_unique_taxa))
        .collect();
    let overview = compute_dataset_overview(&summaries);
    (summaries, overview)
}

/// Evaluates a recipe across alignments. The summaries include the per-sample
/// lists.
#[cfg(not(target_arch = "wasm32"))]
pub fn evaluate_recipe_on_alignments_with_progress<F>(
    alignments: &[Alignment],
    recipe: &TrimmingRecipe,
    total_unique_taxa: usize,
    progress_callback: Option<F>,
) -> (Vec<AlignmentSummary>, DatasetOverview)
where
    F: Fn(usize, usize, &str) + Send + Sync,
{
    let measured =
        measure_alignments(alignments, recipe, total_unique_taxa, progress_callback, || false)
            .expect("a run that cannot be cancelled always finishes");
    let summaries: Vec<AlignmentSummary> = measured
        .into_iter()
        .map(|locus| {
            with_retention(apply_assessment(&locus, recipe, total_unique_taxa), locus.retention)
        })
        .collect();
    let overview = compute_dataset_overview(&summaries);
    (summaries, overview)
}

#[cfg(all(test, not(target_arch = "wasm32")))]
mod tests {
    use super::*;
    use crate::algorithms::orf::StopCodonAction;
    use crate::models::AlignmentFormat;

    #[test]
    fn scan_reads_each_format_and_reports_unreadable_files() {
        let root = std::env::temp_dir().join("alignmentforge_scan_formats_test");
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(&root).unwrap();
        std::fs::write(root.join("locus1.fa"), ">a\nACGTACGT\n>b\nACGTACGA\n").unwrap();
        std::fs::write(root.join("locus2.phy"), "2 6\na ACGTAC\nc ACGTAA\n").unwrap();
        std::fs::write(
            root.join("locus3.nex"),
            "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=3 NCHAR=4;\nFORMAT DATATYPE=DNA;\nMATRIX\na ACGT\nb ACGA\nc ACGG\n;\nEND;\n",
        )
        .unwrap();
        std::fs::write(root.join("broken.fa"), "ACGT\n>a\nACGT\n").unwrap();

        let (alignments, occupancy, failures) =
            scan_alignment_directory(&root, None::<fn(usize, usize, &str)>).unwrap();

        let mut ids: Vec<&str> = alignments.iter().map(|alignment| alignment.id.as_str()).collect();
        ids.sort();
        assert_eq!(ids, vec!["locus1", "locus2", "locus3"]);
        assert_eq!(occupancy.len(), 3);
        assert_eq!(failures.len(), 1, "{failures:?}");
        assert_eq!(failures[0].file_name, "broken.fa");
        let loci_with = |taxon: &str| {
            occupancy
                .iter()
                .find(|entry| entry.taxon_name == taxon)
                .map_or(0, |entry| entry.present_loci_count)
        };
        assert_eq!([loci_with("a"), loci_with("b"), loci_with("c")], [3, 2, 2]);

        let _ = std::fs::remove_dir_all(&root);
    }

    /// The browser build reads each file as text and applies the pass thresholds
    /// as a separate step (`wasm_api.rs`). The desktop build reads files from disk
    /// and assesses each locus in one pass. Both must give the same summaries.
    #[test]
    fn browser_and_desktop_paths_give_the_same_summaries() {
        use crate::parsers::parse_alignment_text;

        let mut paths: Vec<PathBuf> = std::fs::read_dir("../public/example_data/all_markers")
            .expect("the bundled example data is missing")
            .map(|entry| entry.unwrap().path())
            .filter(|path| path.extension().is_some_and(|extension| extension == "phy"))
            .collect();
        paths.sort();
        assert!(!paths.is_empty());

        let mut desktop_alignments = Vec::new();
        let mut browser_alignments = Vec::new();
        for path in &paths {
            let file_name = path.file_name().unwrap().to_string_lossy().to_string();
            // The browser names a locus by its file name without the last
            // extension (src/tauriClient.ts).
            let id = file_name.rsplit_once('.').map_or(file_name.as_str(), |(stem, _)| stem);
            let content = std::fs::read_to_string(path).unwrap();
            let browser =
                parse_alignment_text(&content, id, &file_name, &path.to_string_lossy()).unwrap();
            let desktop = parse_alignment(path).unwrap();
            assert_eq!(
                serde_json::to_value(&desktop).unwrap(),
                serde_json::to_value(&browser).unwrap(),
                "{file_name}"
            );
            desktop_alignments.push(desktop);
            browser_alignments.push(browser);
        }
        let dataset_taxa: Vec<Vec<String>> =
            browser_alignments.iter().map(|alignment| alignment.taxa.clone()).collect();
        let total_unique_taxa = dataset_taxa.iter().flatten().collect::<HashSet<_>>().len();

        let desktop_summaries = |recipe: &TrimmingRecipe| {
            evaluate_recipe_on_alignments_with_progress(
                &desktop_alignments,
                recipe,
                total_unique_taxa,
                None::<fn(usize, usize, &str)>,
            )
            .0
        };
        // The browser measures the loci with one recipe and can apply the pass
        // thresholds of another recipe to the stored measurements. It fetches
        // the per-sample lists with a separate call.
        let browser_summaries = |measured_with: &TrimmingRecipe, assessed_with: &TrimmingRecipe| {
            let measure =
                CatalogRecipes::new(&recipe_with_taxon_presence(measured_with, &dataset_taxa));
            let assess = recipe_with_taxon_presence(assessed_with, &dataset_taxa);
            browser_alignments
                .iter()
                .map(|alignment| {
                    let measured = summarize_alignment_unassessed(alignment, &measure, total_unique_taxa);
                    let summary = apply_assessment(&measured, &assess, total_unique_taxa);
                    with_retention(summary, measured.retention)
                })
                .collect::<Vec<_>>()
        };

        let default = TrimmingRecipe::default();
        let with_orf = TrimmingRecipe { enable_orf: true, ..TrimmingRecipe::default() };
        for recipe in [&default, &TrimmingRecipe::preset_strict(), &with_orf] {
            assert_eq!(
                serde_json::to_value(desktop_summaries(recipe)).unwrap(),
                serde_json::to_value(browser_summaries(recipe, recipe)).unwrap(),
                "{}",
                recipe.name
            );
        }

        // When only the thresholds change, the browser uses the stored
        // measurements again (the gating fields of `recipeCacheKey` in
        // src/engine/enginePool.ts).
        let mut stricter = default.clone();
        stricter.min_taxa = 40;
        stricter.min_taxa_occupancy_percent = 80.0;
        stricter.min_length = 2000;
        stricter.max_gap_percent = 30.0;
        stricter.min_pis_count = 20;
        stricter.min_pis_percent = 2.0;
        stricter.min_variable_count = 40;
        stricter.min_variable_percent = 4.0;
        let passing = |summaries: &[AlignmentSummary]| summaries.iter().filter(|summary| summary.pass).count();
        let stricter_desktop = desktop_summaries(&stricter);
        assert!(passing(&stricter_desktop) < passing(&desktop_summaries(&default)));
        assert_eq!(
            serde_json::to_value(&stricter_desktop).unwrap(),
            serde_json::to_value(browser_summaries(&default, &stricter)).unwrap()
        );
    }

    #[test]
    fn scan_reads_three_subfolder_levels_and_skips_hidden_entries() {
        // The chosen folder itself starts with `.`, which must not stop the scan.
        let root = std::env::temp_dir().join(".af_scan_depth_test");
        let _ = std::fs::remove_dir_all(&root);
        for relative in [
            "top.fa",
            "a/b/c/deep.fa",
            "a/b/c/d/too_deep.fa",
            ".hidden/skipped.fa",
            "__MACOSX/a/skipped.fa",
            "a/._skipped.fa",
        ] {
            let path = root.join(relative);
            std::fs::create_dir_all(path.parent().unwrap()).unwrap();
            std::fs::write(&path, ">x\nACGT\n>y\nACGA\n").unwrap();
        }

        let (alignments, _, failures) =
            scan_alignment_directory(&root, None::<fn(usize, usize, &str)>).unwrap();
        let mut ids: Vec<&str> = alignments.iter().map(|alignment| alignment.id.as_str()).collect();
        ids.sort();
        assert_eq!(ids, vec!["deep", "top"]);
        assert!(failures.is_empty(), "{failures:?}");

        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn test_orf_sample_pruning_does_not_change_catalog_assessment() {
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
        assert_eq!(summaries[0].num_taxa, 4);
        assert!(summaries[0].pass);
        assert!(summaries[0].fail_reasons.is_empty());
        assert_eq!(summaries[0].orf_retained_samples, 3);
        assert_eq!(summaries[0].orf_retained_taxa.len(), 3);
        assert_eq!(overview.passed_alignments, 1);
    }

    #[test]
    fn discarded_samples_survive_the_occupancy_filter() {
        let alignments = vec![Alignment::new(
            "locus_1".to_string(),
            "locus_1.fa".to_string(),
            "/dummy/locus_1.fa".to_string(),
            AlignmentFormat::Fasta,
            vec!["Keep".to_string(), "DropByHand".to_string()],
            vec!["ACGT".to_string(), "ACGT".to_string()],
        )];

        let mut recipe = TrimmingRecipe {
            discarded_taxa: vec!["DropByHand".to_string()],
            ..Default::default()
        };
        // With no occupancy threshold the manual list must still be applied.
        let runtime = recipe_with_dataset_sample_filter(&recipe, &alignments);
        assert_eq!(runtime.excluded_taxa, vec!["DropByHand".to_string()]);

        // And it must survive a run that also derives exclusions from occupancy.
        recipe.min_sample_locus_occupancy_percent = 60.0;
        let runtime = recipe_with_dataset_sample_filter(&recipe, &alignments);
        assert!(runtime.excluded_taxa.contains(&"DropByHand".to_string()));
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

        let (summaries, _) = evaluate_recipe_on_alignments_with_progress(
            &alignments,
            &recipe,
            2,
            None::<fn(usize, usize, &str)>,
        );
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

        let result = measure_alignments(
            &[alignment],
            &TrimmingRecipe::default(),
            1,
            None::<fn(usize, usize, &str)>,
            || true,
        );

        assert!(result.is_none());
    }
}
