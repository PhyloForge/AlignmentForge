use serde::{Deserialize, Serialize};
use std::path::PathBuf;
use tauri::Emitter;

use crate::export::batch::{execute_batch_export, BatchExportConfig, BatchExportResult};
use crate::export::concatenate::{concatenate_alignments, ConcatenateConfig, ConcatenateResult};
use crate::export::group::{concatenate_alignments_by_gene, GroupedConcatenateConfig, GroupedConcatenateResult};
use crate::filter_config::{load_filter_config as read_filter_config, save_filter_config as write_filter_config};
use crate::models::{
    Alignment, AlignmentFormat, AlignmentSummary, DatasetOverview, TaxonOccupancy, TrimmingDiff,
};
use crate::parsers::parse_alignment;
use crate::pipeline::catalog::{
    evaluate_recipe_on_alignments_with_progress_and_cancel,
    evaluate_recipe_on_summaries_with_progress, recipe_with_taxon_presence,
    scan_alignment_directory, ParseFailure,
};
use crate::pipeline::engine::apply_recipe;
use crate::pipeline::recipe::TrimmingRecipe;
use crate::state::AlignmentCache;

/// Runs pipeline work and turns a panic into an error the user can see.
///
/// The release profile unwinds, so one malformed alignment reports a failed
/// command instead of closing the application.
fn guard_panics<T, F>(what: &str, work: F) -> Result<T, String>
where
    F: FnOnce() -> Result<T, String>,
{
    // The Tauri handle is not `UnwindSafe`, but nothing here observes state
    // after a panic: the result is discarded and an error is returned instead.
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(work)) {
        Ok(result) => result,
        Err(payload) => {
            let detail = payload
                .downcast_ref::<&str>()
                .map(|s| (*s).to_string())
                .or_else(|| payload.downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "unknown error".to_string());
            Err(format!("{what} failed unexpectedly: {detail}"))
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProgressPayload {
    pub current: usize,
    pub total: usize,
    pub percent: f64,
    pub file_name: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct CatalogProgressPayload {
    pub job_id: String,
    pub current: usize,
    pub total: usize,
    pub percent: f64,
    pub file_name: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScanResponse {
    pub summaries: Vec<AlignmentSummary>,
    pub overview: DatasetOverview,
    pub occupancy: Vec<TaxonOccupancy>,
    /// Files found in the folder that could not be parsed.
    #[serde(default)]
    pub parse_failures: Vec<ParseFailure>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentViewResponse {
    pub raw_alignment: Alignment,
    pub trimmed_alignment: Alignment,
    pub diff: TrimmingDiff,
    pub pis_mask: Vec<bool>,
    pub majority_consensus: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CatalogUpdateResponse {
    pub summaries: Vec<AlignmentSummary>,
    pub overview: DatasetOverview,
}

#[tauri::command]
pub async fn scan_directory(
    app: tauri::AppHandle,
    state: tauri::State<'_, AlignmentCache>,
    dir_path: String,
) -> Result<ScanResponse, String> {
    let cache = state.inner().clone();
    tokio::task::spawn_blocking(move || {
      guard_panics("Directory scan", move || {
        let app_clone = app.clone();
        let (summaries, overview, occupancy, alignments, parse_failures) =
            scan_alignment_directory(&dir_path, Some(move |cur, total, name: &str| {
                let pct = if total > 0 {
                    (cur as f64 / total as f64) * 100.0
                } else {
                    0.0
                };
                let _ = app_clone.emit(
                    "scan_progress",
                    ProgressPayload {
                        current: cur,
                        total,
                        percent: pct,
                        file_name: name.to_string(),
                    },
                );
            }))?;

        cache.store(alignments);

        Ok(ScanResponse {
            summaries,
            overview,
            occupancy,
            parse_failures,
        })
      })
    })
    .await
    .map_err(|e| e.to_string())?
}

#[tauri::command]
pub async fn get_alignment(
    state: tauri::State<'_, AlignmentCache>,
    file_path: String,
    recipe: TrimmingRecipe,
    total_unique_taxa: usize,
) -> Result<AlignmentViewResponse, String> {
    let cache = state.inner().clone();
    tokio::task::spawn_blocking(move || {
      guard_panics("Opening the alignment", move || {
        let raw = match cache.get(&file_path) {
            Some(align) => align,
            None => parse_alignment(PathBuf::from(&file_path))?,
        };
        // Only the taxon lists matter here, so avoid copying the whole dataset's
        // sequences every time the user opens an alignment.
        let mut dataset_taxa = cache.taxon_presence();
        if !cache.contains(&raw.file_path) {
            dataset_taxa.push(raw.taxa.clone());
        }
        let runtime_recipe = recipe_with_taxon_presence(&recipe, &dataset_taxa);
        let (trimmed, diff) = apply_recipe(&raw, &runtime_recipe, total_unique_taxa);

        let (_, _, pis_mask) =
            crate::algorithms::informative::calculate_parsimony_informative_sites(&raw.sequences, true);
        let majority_consensus =
            crate::algorithms::stats::compute_majority_consensus(&raw.sequences, false);

        Ok(AlignmentViewResponse {
            raw_alignment: raw,
            trimmed_alignment: trimmed,
            diff,
            pis_mask,
            majority_consensus,
        })
      })
    })
    .await
    .map_err(|e| e.to_string())?
}

#[tauri::command]
pub async fn recalculate_catalog(
    app: tauri::AppHandle,
    state: tauri::State<'_, AlignmentCache>,
    paths: Vec<String>,
    recipe: TrimmingRecipe,
    total_unique_taxa: usize,
    job_id: String,
) -> Result<CatalogUpdateResponse, String> {
    let cache = state.inner().clone();
    let generation = cache.begin_catalog_job();
    tokio::task::spawn_blocking(move || {
      guard_panics("Catalog recalculation", move || {
        let cached = cache.get_by_paths(&paths);
        let app_clone = app.clone();
        let progress_job_id = job_id.clone();
        let progress_callback = move |current: usize, total: usize, file_name: &str| {
            let emit_every = (total / 100).max(1);
            if current == total || current % emit_every == 0 {
                let percent = if total > 0 {
                    (current as f64 / total as f64) * 100.0
                } else {
                    100.0
                };
                let _ = app_clone.emit(
                    "catalog_recalculation_progress",
                    CatalogProgressPayload {
                        job_id: progress_job_id.clone(),
                        current,
                        total,
                        percent,
                        file_name: file_name.to_string(),
                    },
                );
            }
        };
        let (summaries, mut overview) = if cached.len() == paths.len() && !cached.is_empty() {
            evaluate_recipe_on_alignments_with_progress_and_cancel(
                &cached,
                &recipe,
                total_unique_taxa,
                Some(progress_callback),
                || !cache.is_catalog_job_current(generation),
            )
            .ok_or_else(|| "Catalog recalculation superseded by a newer request".to_string())?
        } else {
            if !cache.is_catalog_job_current(generation) {
                return Err("Catalog recalculation superseded by a newer request".to_string());
            }
            evaluate_recipe_on_summaries_with_progress(
                &paths,
                &recipe,
                total_unique_taxa,
                Some(progress_callback),
            )
        };
        if !cache.is_catalog_job_current(generation) {
            return Err("Catalog recalculation superseded by a newer request".to_string());
        }
        overview.total_unique_taxa = total_unique_taxa;
        Ok(CatalogUpdateResponse { summaries, overview })
      })
    })
    .await
    .map_err(|e| e.to_string())?
}

/// Gives an export the dataset context of the loaded folder: the sample
/// exclusions and the number of unique samples. With no folder loaded, it
/// returns `None`, and the export reads the context from its own input files.
fn export_context(cache: &AlignmentCache, recipe: &TrimmingRecipe) -> (TrimmingRecipe, Option<usize>) {
    let taxon_presence = cache.taxon_presence();
    if taxon_presence.is_empty() {
        return (recipe.clone(), None);
    }
    (
        recipe_with_taxon_presence(recipe, &taxon_presence),
        Some(cache.unique_taxon_count()),
    )
}

#[tauri::command]
pub async fn run_batch_export(
    state: tauri::State<'_, AlignmentCache>,
    config: BatchExportConfig,
    recipe: TrimmingRecipe,
) -> Result<BatchExportResult, String> {
    let cache = state.inner().clone();
    tokio::task::spawn_blocking(move || {
      guard_panics("Batch export", move || {
        let (runtime_recipe, dataset_taxa) = export_context(&cache, &recipe);
        execute_batch_export(&config, &runtime_recipe, dataset_taxa)
      })
    })
        .await
        .map_err(|e| e.to_string())?
}

#[tauri::command]
pub async fn run_concatenate(
    state: tauri::State<'_, AlignmentCache>,
    config: ConcatenateConfig,
    recipe: TrimmingRecipe,
) -> Result<ConcatenateResult, String> {
    let cache = state.inner().clone();
    tokio::task::spawn_blocking(move || {
      guard_panics("Supermatrix assembly", move || {
        let (runtime_recipe, dataset_taxa) = export_context(&cache, &recipe);
        concatenate_alignments(&config, &runtime_recipe, dataset_taxa)
      })
    })
        .await
        .map_err(|e| e.to_string())?
}

#[tauri::command]
pub async fn run_grouped_concatenate(
    state: tauri::State<'_, AlignmentCache>,
    config: GroupedConcatenateConfig,
    recipe: TrimmingRecipe,
) -> Result<GroupedConcatenateResult, String> {
    let cache = state.inner().clone();
    tokio::task::spawn_blocking(move || {
      guard_panics("Gene concatenation", move || {
        let (runtime_recipe, dataset_taxa) = export_context(&cache, &recipe);
        concatenate_alignments_by_gene(&config, &runtime_recipe, dataset_taxa)
      })
    })
        .await
        .map_err(|e| e.to_string())?
}

#[tauri::command]
pub fn save_filter_config(
    file_path: String,
    recipe: TrimmingRecipe,
) -> Result<String, String> {
    write_filter_config(PathBuf::from(file_path).as_path(), &recipe)
}

#[tauri::command]
pub fn load_filter_config(file_path: String) -> Result<TrimmingRecipe, String> {
    read_filter_config(PathBuf::from(file_path).as_path())
}

#[tauri::command]
pub fn save_alignment_stats_csv(
    file_path: String,
    summaries: Vec<AlignmentSummary>,
) -> Result<String, String> {
    let mut writer = csv::Writer::from_path(&file_path)
        .map_err(|error| format!("Could not create alignment statistics CSV: {error}"))?;

    writer
        .write_record([
            "locus_id",
            "file_name",
            "file_path",
            "format",
            "num_taxa",
            "length_bp",
            "total_basepairs",
            "gap_count",
            "gap_percent",
            "variable_sites",
            "variable_sites_percent",
            "parsimony_informative_sites",
            "parsimony_informative_percent",
            "mean_divergence",
            "gc_percent",
            "orf_evaluated",
            "orf_candidate_found",
            "orf_valid",
            "orf_frame",
            "orf_raw_start_1based",
            "orf_raw_end",
            "orf_support_count",
            "orf_support_percent",
            "orf_retained_samples",
            "orf_candidate_length_aa",
            "orf_coding_score",
            "orf_amino_acid_conservation_percent",
            "orf_frame_contrast_percent",
            "orf_reference_evaluated",
            "orf_reference_matched",
            "orf_reference_identity_percent",
            "orf_reference_coverage_percent",
            "orf_intron_length_bp",
            "pass",
            "fail_reasons",
            "raw_num_taxa",
            "raw_length_bp",
            "raw_gap_percent",
        ])
        .map_err(|error| format!("Could not write alignment statistics header: {error}"))?;

    for summary in summaries {
        let format_name = match summary.format {
            AlignmentFormat::Fasta => "fasta",
            AlignmentFormat::Phylip => "phylip",
            AlignmentFormat::Nexus => "nexus",
        };
        writer
            .write_record([
                summary.id,
                summary.file_name,
                summary.file_path,
                format_name.to_string(),
                summary.num_taxa.to_string(),
                summary.length.to_string(),
                summary.total_basepairs.to_string(),
                summary.gap_count.to_string(),
                format!("{:.3}", summary.gap_percent),
                summary.variable_count.to_string(),
                format!("{:.3}", summary.variable_percent),
                summary.pis_count.to_string(),
                format!("{:.3}", summary.pis_percent),
                format!("{:.6}", summary.mean_divergence),
                format!("{:.3}", summary.gc_percent),
                summary.orf_evaluated.to_string(),
                summary.orf_candidate_found.to_string(),
                summary.orf_valid.to_string(),
                summary.orf_frame.map(|value| value.to_string()).unwrap_or_default(),
                summary
                    .orf_start
                    .map(|value| (value + 1).to_string())
                    .unwrap_or_default(),
                summary.orf_end.map(|value| value.to_string()).unwrap_or_default(),
                summary.orf_support_count.to_string(),
                format!("{:.3}", summary.orf_support_percent),
                summary.orf_retained_samples.to_string(),
                summary.orf_candidate_length_aa.to_string(),
                format!("{:.3}", summary.orf_coding_score),
                format!("{:.3}", summary.orf_amino_acid_conservation),
                format!("{:.3}", summary.orf_frame_contrast),
                summary.orf_reference_evaluated.to_string(),
                summary.orf_reference_matched.to_string(),
                format!("{:.3}", summary.orf_reference_identity),
                format!("{:.3}", summary.orf_reference_coverage),
                summary.orf_intron_length.to_string(),
                summary.pass.to_string(),
                summary.fail_reasons.join(" | "),
                summary.raw_num_taxa.to_string(),
                summary.raw_length.to_string(),
                format!("{:.3}", summary.raw_gap_percent),
            ])
            .map_err(|error| format!("Could not write alignment statistics row: {error}"))?;
    }

    writer
        .flush()
        .map_err(|error| format!("Could not finish alignment statistics CSV: {error}"))?;
    Ok(file_path)
}

#[tauri::command]
pub fn get_presets() -> Vec<TrimmingRecipe> {
    vec![
        TrimmingRecipe::default(),
        TrimmingRecipe::preset_strict(),
        TrimmingRecipe::preset_relaxed_uce(),
        TrimmingRecipe::preset_exon_codon(),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn export_context_comes_from_the_loaded_folder_when_there_is_one() {
        let recipe = TrimmingRecipe {
            min_sample_locus_occupancy_percent: 60.0,
            ..TrimmingRecipe::default()
        };
        let cache = AlignmentCache::new();

        let (runtime_recipe, dataset_taxa) = export_context(&cache, &recipe);
        assert_eq!(dataset_taxa, None);
        assert!(runtime_recipe.excluded_taxa.is_empty());

        let locus = |id: &str, taxa: [&str; 2]| {
            Alignment::new(
                id.to_string(),
                format!("{id}.fa"),
                format!("/data/{id}.fa"),
                AlignmentFormat::Fasta,
                taxa.iter().map(|taxon| taxon.to_string()).collect(),
                vec!["ACGT".to_string(); 2],
            )
        };
        cache.store(vec![locus("locus1", ["a", "b"]), locus("locus2", ["a", "c"])]);

        let (runtime_recipe, dataset_taxa) = export_context(&cache, &recipe);
        assert_eq!(dataset_taxa, Some(3));
        // b and c are each in 50% of the loci, below the 60% threshold.
        assert_eq!(runtime_recipe.excluded_taxa, vec!["b", "c"]);
    }
}
