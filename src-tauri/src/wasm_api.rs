//! JavaScript boundary for the browser build.
//!
//! The browser runs the same trimming engine as the desktop app, compiled to
//! WebAssembly. Each call takes an alignment as text and returns JSON, so a
//! worker can own a shard of the dataset without sharing Rust state.

use wasm_bindgen::prelude::*;

use crate::models::Alignment;
use crate::parsers::parse_alignment_text;
use crate::pipeline::catalog::{
    apply_assessment, recipe_with_taxon_presence, summarize_alignment_unassessed, CatalogRecipes,
    SampleRetention, UnassessedSummary,
};
use crate::pipeline::engine::apply_recipe;
use crate::pipeline::recipe::TrimmingRecipe;

/// A previously computed shard result, kept so that toggling a filter back to
/// an earlier setting returns immediately instead of recomputing every locus.
struct CachedRun {
    key: String,
    summaries: Vec<UnassessedSummary>,
}

#[wasm_bindgen]
pub struct EngineSession {
    /// Parsed alignments held for the lifetime of a worker, in load order.
    alignments: Vec<Alignment>,
    /// Measured loci accumulated by the current run, before gating is applied.
    pending: Vec<UnassessedSummary>,
    /// The full recipe for the current run, including its pass thresholds.
    runtime_recipe: Option<TrimmingRecipe>,
    /// The recipes that measure each locus in the current run. Built once per
    /// run, because a recipe can hold thousands of reference sequences.
    measure_recipes: Option<CatalogRecipes>,
    cache: Vec<CachedRun>,
    /// Every locus's taxon list across the whole dataset, set once after
    /// loading. Passing it on each call would copy it into every worker on
    /// every recipe change, which dominates the cost on a large dataset.
    dataset_taxa: Vec<Vec<String>>,
    total_unique_taxa: usize,
}

fn to_js_error(message: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&message.to_string())
}

/// Serialises to plain JavaScript objects.
///
/// The default serializer turns a Rust map into a JavaScript `Map`, which the
/// interface cannot index like an object, so every map-valued field has to go
/// through the JSON compatible serializer instead.
fn to_js<T: serde::Serialize + ?Sized>(value: &T) -> Result<JsValue, JsValue> {
    let serializer = serde_wasm_bindgen::Serializer::json_compatible();
    value.serialize(&serializer).map_err(to_js_error)
}

#[wasm_bindgen]
impl EngineSession {
    #[wasm_bindgen(constructor)]
    pub fn new() -> EngineSession {
        EngineSession {
            alignments: Vec::new(),
            pending: Vec::new(),
            runtime_recipe: None,
            measure_recipes: None,
            cache: Vec::new(),
            dataset_taxa: Vec::new(),
            total_unique_taxa: 0,
        }
    }

    /// Parses one alignment and keeps it. Returns its taxon names so the caller
    /// can build dataset-wide occupancy without shipping the sequences back.
    #[wasm_bindgen(js_name = addAlignment)]
    pub fn add_alignment(
        &mut self,
        content: &str,
        id: &str,
        file_name: &str,
        file_path: &str,
    ) -> Result<JsValue, JsValue> {
        let alignment =
            parse_alignment_text(content, id, file_name, file_path).map_err(to_js_error)?;
        let taxa = alignment.taxa.clone();
        self.alignments.push(alignment);
        to_js(&taxa)
    }

    #[wasm_bindgen(js_name = alignmentCount)]
    pub fn alignment_count(&self) -> usize {
        self.alignments.len()
    }

    /// Records the dataset-wide taxon lists once, after every shard has loaded.
    #[wasm_bindgen(js_name = setDataset)]
    pub fn set_dataset(
        &mut self,
        dataset_taxa: JsValue,
        total_unique_taxa: usize,
    ) -> Result<(), JsValue> {
        self.dataset_taxa = serde_wasm_bindgen::from_value(dataset_taxa).map_err(to_js_error)?;
        self.total_unique_taxa = total_unique_taxa;
        Ok(())
    }

    /// Prepares a summarisation run.
    ///
    /// Returns true when this exact recipe was computed before, in which case
    /// the caller can skip straight to `summarizeFinish`. Recomputing every
    /// locus to return to a setting the user just left is pure waste.
    #[wasm_bindgen(js_name = summarizeBegin)]
    pub fn summarize_begin(&mut self, recipe: JsValue, cache_key: &str) -> Result<bool, JsValue> {
        self.pending.clear();
        // The recipe is always needed, because the thresholds are applied even
        // when the measured loci come from the cache.
        let recipe: TrimmingRecipe =
            serde_wasm_bindgen::from_value(recipe).map_err(to_js_error)?;
        let runtime_recipe = recipe_with_taxon_presence(&recipe, &self.dataset_taxa);
        let cached = self.cache.iter().any(|entry| entry.key == cache_key);
        self.measure_recipes = (!cached).then(|| CatalogRecipes::new(&runtime_recipe));
        self.runtime_recipe = Some(runtime_recipe);
        Ok(cached)
    }

    /// Summarises one range of this shard into the pending buffer.
    ///
    /// Nothing crosses into JavaScript here, so the caller can use small ranges
    /// to report progress without paying to serialise each one.
    #[wasm_bindgen(js_name = summarizeChunk)]
    pub fn summarize_chunk(&mut self, start: usize, count: usize) -> Result<(), JsValue> {
        let recipes = self
            .measure_recipes
            .as_ref()
            .ok_or_else(|| to_js_error("summarizeBegin was not called"))?;
        let end = start.saturating_add(count).min(self.alignments.len());
        for alignment in &self.alignments[start.min(end)..end] {
            self.pending.push(summarize_alignment_unassessed(
                alignment,
                recipes,
                self.total_unique_taxa,
            ));
        }
        Ok(())
    }

    /// Applies the thresholds and returns the run's summaries.
    ///
    /// Gating runs here, not in the pipeline. When only a pass threshold
    /// changes, the measured loci are used again and only pass and fail change.
    /// The result does not include the per-sample retention lists, because only
    /// the QC occupancy chart uses them.
    #[wasm_bindgen(js_name = summarizeFinish)]
    pub fn summarize_finish(&mut self, cache_key: &str) -> Result<JsValue, JsValue> {
        let recipe = self
            .runtime_recipe
            .take()
            .ok_or_else(|| to_js_error("summarizeBegin was not called"))?;
        self.measure_recipes = None;

        if let Some(index) = self.cache.iter().position(|entry| entry.key == cache_key) {
            let hit = self.cache.remove(index);
            let value = to_js(&self.assess_all(&hit.summaries, &recipe))?;
            self.cache.insert(0, hit);
            return Ok(value);
        }

        let summaries = std::mem::take(&mut self.pending);
        let value = to_js(&self.assess_all(&summaries, &recipe))?;
        self.cache.insert(
            0,
            CachedRun {
                key: cache_key.to_string(),
                summaries,
            },
        );
        // A handful of recent recipes is enough to make toggling a filter feel
        // immediate without holding the whole history in memory.
        self.cache.truncate(4);
        Ok(value)
    }

    /// Per-sample retention for the QC occupancy chart, fetched only when that
    /// view needs it.
    #[wasm_bindgen(js_name = retentionDetails)]
    pub fn retention_details(&self, cache_key: &str) -> Result<JsValue, JsValue> {
        let entry = self
            .cache
            .iter()
            .find(|entry| entry.key == cache_key)
            .ok_or_else(|| to_js_error("No summaries are cached for that recipe"))?;
        let details: Vec<&SampleRetention> =
            entry.summaries.iter().map(|locus| &locus.retention).collect();
        to_js(&details)
    }

    /// Per-taxon presence and coverage for this worker's alignments.
    ///
    /// The caller merges the shards into dataset-wide occupancy. Returning
    /// counts rather than sequences keeps the alignments inside the worker.
    #[wasm_bindgen(js_name = taxonStats)]
    pub fn taxon_stats(&self) -> Result<JsValue, JsValue> {
        #[derive(serde::Serialize)]
        struct TaxonStat {
            taxon: String,
            basepairs: usize,
            gap_percent: f64,
        }

        let mut stats = Vec::new();
        for alignment in &self.alignments {
            for (taxon, sequence) in alignment.taxa.iter().zip(alignment.sequences.iter()) {
                let gaps = sequence
                    .bytes()
                    .filter(|state| matches!(state, b'-' | b'?' | b'N' | b'n'))
                    .count();
                let length = sequence.len();
                stats.push(TaxonStat {
                    taxon: taxon.clone(),
                    basepairs: length.saturating_sub(gaps),
                    gap_percent: if length > 0 {
                        (gaps as f64 / length as f64) * 100.0
                    } else {
                        100.0
                    },
                });
            }
        }
        to_js(&stats)
    }

    /// One alignment for the viewer: trimmed, and the diff. The unprocessed
    /// locus does not change with the recipe, so it is added only on request.
    #[wasm_bindgen(js_name = viewAlignment)]
    pub fn view_alignment(
        &self,
        file_path: &str,
        recipe: JsValue,
        include_raw: bool,
    ) -> Result<JsValue, JsValue> {
        #[derive(serde::Serialize)]
        struct RawView<'a> {
            raw_alignment: &'a Alignment,
            pis_mask: Vec<bool>,
            majority_consensus: String,
        }
        #[derive(serde::Serialize)]
        struct View<'a> {
            trimmed_alignment: Alignment,
            /// A JSON value, which writes the column numbers of `column_reasons`
            /// as strings. The JavaScript serializer accepts only string keys.
            diff: serde_json::Value,
            #[serde(skip_serializing_if = "Option::is_none")]
            raw: Option<RawView<'a>>,
        }

        let recipe: TrimmingRecipe =
            serde_wasm_bindgen::from_value(recipe).map_err(to_js_error)?;
        let runtime = recipe_with_taxon_presence(&recipe, &self.dataset_taxa);
        let total_unique_taxa = self.total_unique_taxa;

        let raw = self
            .alignments
            .iter()
            .find(|alignment| alignment.file_path == file_path)
            .ok_or_else(|| to_js_error(format!("Alignment not found: {file_path}")))?;

        let (trimmed, diff) = apply_recipe(raw, &runtime, total_unique_taxa);
        let raw_view = include_raw.then(|| {
            let (_, _, pis_mask) =
                crate::algorithms::informative::calculate_parsimony_informative_sites(
                    &raw.sequences,
                    true,
                );
            RawView {
                raw_alignment: raw,
                pis_mask,
                majority_consensus: crate::algorithms::stats::compute_majority_consensus(
                    &raw.sequences,
                    false,
                ),
            }
        });
        to_js(&View {
            trimmed_alignment: trimmed,
            diff: serde_json::to_value(&diff).map_err(to_js_error)?,
            raw: raw_view,
        })
    }
}

impl Default for EngineSession {
    fn default() -> Self {
        Self::new()
    }
}

impl EngineSession {
    /// Applies the thresholds. The summaries leave out the per-sample lists,
    /// which `retentionDetails` returns.
    fn assess_all(
        &self,
        summaries: &[UnassessedSummary],
        recipe: &TrimmingRecipe,
    ) -> Vec<crate::models::AlignmentSummary> {
        summaries
            .iter()
            .map(|entry| apply_assessment(entry, recipe, self.total_unique_taxa))
            .collect()
    }
}

/// The engine's built-in presets, so the browser shows the same list as desktop.
#[wasm_bindgen(js_name = enginePresets)]
pub fn engine_presets() -> Result<JsValue, JsValue> {
    let presets = vec![
        TrimmingRecipe::default(),
        TrimmingRecipe::preset_strict(),
        TrimmingRecipe::preset_relaxed_uce(),
        TrimmingRecipe::preset_exon_codon(),
    ];
    to_js(&presets)
}
