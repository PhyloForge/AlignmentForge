use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::export::write_alignment_atomic;
use crate::models::{Alignment, AlignmentFormat};
use crate::parsers::parse_alignment;
use crate::pipeline::catalog::{
    recipe_with_dataset_sample_filter, recipe_without_orf_analysis,
};
use crate::pipeline::engine::apply_recipe;
use crate::pipeline::recipe::TrimmingRecipe;
use crate::algorithms::orf::OrfSearchMode;
use crate::algorithms::reference::intron_alignment_outside;

fn default_true() -> bool {
    true
}

fn default_general_alignment_directory_name() -> String {
    "all_alignments".to_string()
}

fn default_orf_alignment_directory_name() -> String {
    "orf_alignments".to_string()
}

fn default_intron_directory_name() -> String {
    "intron_alignments".to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchExportConfig {
    pub input_paths: Vec<String>,
    pub output_directory: String,
    #[serde(default = "default_general_alignment_directory_name")]
    pub general_alignment_directory_name: String,
    #[serde(default = "default_orf_alignment_directory_name")]
    pub orf_alignment_directory_name: String,
    #[serde(default = "default_intron_directory_name")]
    pub intron_directory_name: String,
    pub output_format: AlignmentFormat,
    pub only_passing: bool,
    #[serde(default = "default_true")]
    pub export_general_alignments: bool,
    #[serde(default = "default_true")]
    pub export_orf_alignments: bool,
    pub save_recipe_json: bool,
    pub save_summary_csv: bool,
    #[serde(default)]
    pub export_introns: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchExportResult {
    pub total_processed: usize,
    pub total_exported: usize,
    pub total_discarded: usize,
    pub total_orfs_exported: usize,
    pub alignment_directory_path: Option<String>,
    pub orf_directory_path: Option<String>,
    pub summary_csv_path: Option<String>,
    pub recipe_json_path: Option<String>,
    /// FASTA file with the recipe's reference sequences, which `recipe.json` leaves out.
    pub reference_fasta_path: Option<String>,
    pub total_introns_exported: usize,
    pub intron_directory_path: Option<String>,
    pub total_failed: usize,
    pub errors: Vec<BatchExportError>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchExportError {
    pub input_path: String,
    pub output_path: String,
    pub error: String,
}

struct ExportRecord {
    id: String,
    catalog_pass: bool,
    general_exported: bool,
    orf_accepted: bool,
    orf_exported: bool,
    old_taxa: usize,
    catalog_taxa: usize,
    orf_taxa: usize,
    old_length: usize,
    catalog_length: usize,
    orf_length: usize,
    old_gap_percent: f64,
    catalog_gap_percent: f64,
    intron_exported: bool,
    errors: Vec<BatchExportError>,
}

const REFERENCE_FASTA_NAME: &str = "recipe_references.fasta";

/// The contents of `recipe.json`. The recipe does not serialize its reference
/// sequences, so the record names the FASTA file that holds them.
#[derive(Serialize)]
struct RecipeRecord<'a> {
    #[serde(flatten)]
    recipe: &'a TrimmingRecipe,
    #[serde(skip_serializing_if = "Option::is_none")]
    orf_reference_file: Option<&'a str>,
}

/// Writes the recipe's reference sequences as FASTA, which the ORF sidebar can
/// load again. Returns `None` when the recipe has no references.
fn write_reference_fasta(out_dir: &Path, recipe: &TrimmingRecipe) -> Result<Option<PathBuf>, String> {
    if recipe.orf_reference_sequences.is_empty() {
        return Ok(None);
    }
    let mut references: Vec<(String, String)> = recipe
        .orf_reference_sequences
        .iter()
        .map(|(id, sequence)| (id.clone(), sequence.clone()))
        .collect();
    references.sort();
    let (ids, sequences): (Vec<String>, Vec<String>) = references.into_iter().unzip();
    let path = out_dir.join(REFERENCE_FASTA_NAME);
    write_alignment_atomic(&path, &ids, &sequences, AlignmentFormat::Fasta)
        .map_err(|error| format!("Failed to write reference sequences: {error}"))?;
    Ok(Some(path))
}

fn output_path(directory: &Path, alignment: &Alignment, format: AlignmentFormat) -> PathBuf {
    directory.join(format!("{}.{}", alignment.id, format.extension()))
}

fn validate_directory_name<'a>(label: &str, name: &'a str) -> Result<&'a str, String> {
    let name = name.trim();
    if name.is_empty() || name == "." || name == ".." || name.contains('/') || name.contains('\\') {
        return Err(format!(
            "{label} must be a single non-empty folder name without / or \\\\."
        ));
    }
    Ok(name)
}

/// Intron export needs the exon span that a reference ORF mode finds.
fn exports_introns(config: &BatchExportConfig, recipe: &TrimmingRecipe) -> bool {
    config.export_introns
        && recipe.enable_orf
        && recipe.orf_use_references
        && matches!(
            recipe.orf_search_mode,
            OrfSearchMode::ReferenceGuided | OrfSearchMode::ReferenceCandidateOrf
        )
}

pub fn execute_batch_export(
    config: &BatchExportConfig,
    recipe: &TrimmingRecipe,
    resolved_dataset_taxa: Option<usize>,
) -> Result<BatchExportResult, String> {
    let general_directory_name = validate_directory_name(
        "General alignment folder name",
        &config.general_alignment_directory_name,
    )?;
    let orf_directory_name = validate_directory_name(
        "ORF alignment folder name",
        &config.orf_alignment_directory_name,
    )?;
    let intron_directory_name =
        validate_directory_name("Intron folder name", &config.intron_directory_name)?;
    let export_introns = exports_introns(config, recipe);
    let mut active_directory_names = Vec::new();
    if config.export_general_alignments {
        active_directory_names.push(general_directory_name.to_lowercase());
    }
    if config.export_orf_alignments && recipe.enable_orf {
        active_directory_names.push(orf_directory_name.to_lowercase());
    }
    if export_introns {
        active_directory_names.push(intron_directory_name.to_lowercase());
    }
    active_directory_names.sort();
    if active_directory_names
        .windows(2)
        .any(|names| names[0] == names[1])
    {
        return Err("Exported alignment folders must have different names.".to_string());
    }

    let out_dir = Path::new(&config.output_directory);
    if out_dir.exists() && !out_dir.is_dir() {
        return Err(format!(
            "Export destination is not a directory: {}",
            out_dir.display()
        ));
    }
    if !out_dir.exists() {
        fs::create_dir_all(out_dir)
            .map_err(|e| format!("Failed to create output directory: {}", e))?;
    }

    let parsed: Vec<(String, Result<Alignment, String>)> = config
        .input_paths
        .par_iter()
        .map(|input_path| (input_path.clone(), parse_alignment(Path::new(input_path))))
        .collect();
    let mut errors = Vec::new();
    let mut raw_alignments = Vec::new();
    for (input_path, result) in parsed {
        match result {
            Ok(alignment) => raw_alignments.push(alignment),
            Err(error) => errors.push(BatchExportError {
                input_path,
                output_path: String::new(),
                error,
            }),
        }
    }
    let runtime_recipe = if resolved_dataset_taxa.is_none() {
        recipe_with_dataset_sample_filter(recipe, &raw_alignments)
    } else {
        recipe.clone()
    };
    let total_dataset_taxa = resolved_dataset_taxa.unwrap_or_else(|| {
        raw_alignments
            .iter()
            .flat_map(|alignment| alignment.taxa.iter())
            .collect::<std::collections::HashSet<_>>()
            .len()
    });
    let alignment_dir = out_dir.join(general_directory_name);
    let orf_dir = out_dir.join(orf_directory_name);
    let intron_dir = out_dir.join(intron_directory_name);

    let mut destinations: HashMap<String, String> = HashMap::new();
    for alignment in &raw_alignments {
        for directory in [&alignment_dir, &orf_dir, &intron_dir] {
            let destination = output_path(directory, alignment, config.output_format);
            let key = destination.to_string_lossy().to_lowercase();
            if let Some(previous) = destinations.insert(key, alignment.file_path.clone()) {
                return Err(format!(
                    "Inputs '{}' and '{}' resolve to the same output path: {}",
                    previous,
                    alignment.file_path,
                    destination.display()
                ));
            }
        }
    }
    if config.export_general_alignments {
        fs::create_dir_all(&alignment_dir)
            .map_err(|e| format!("Failed to create general alignment output directory: {e}"))?;
    }
    if config.export_orf_alignments && runtime_recipe.enable_orf {
        fs::create_dir_all(&orf_dir)
            .map_err(|e| format!("Failed to create ORF output directory: {e}"))?;
    }
    if export_introns {
        fs::create_dir_all(&intron_dir)
            .map_err(|e| format!("Failed to create intron output directory: {e}"))?;
    }

    let catalog_recipe = recipe_without_orf_analysis(&runtime_recipe);

    // Process all files in parallel using one dataset-wide sample exclusion
    // list, but keep general alignment and ORF outcomes independent.
    let results: Vec<ExportRecord> = raw_alignments
        .par_iter()
        .map(|raw_align| {
            let mut record_errors = Vec::new();
            let (catalog_alignment, catalog_diff) =
                apply_recipe(raw_align, &catalog_recipe, total_dataset_taxa);

            let should_export_general = config.export_general_alignments
                && (!config.only_passing || catalog_diff.pass)
                && !catalog_alignment.sequences.is_empty()
                && catalog_alignment.length > 0;
            let general_path = output_path(&alignment_dir, &catalog_alignment, config.output_format);
            let general_exported = if should_export_general {
                match write_alignment_atomic(
                    &general_path,
                    &catalog_alignment.taxa,
                    &catalog_alignment.sequences,
                    config.output_format,
                ) {
                    Ok(()) => true,
                    Err(error) => {
                        record_errors.push(BatchExportError {
                            input_path: raw_align.file_path.clone(),
                            output_path: general_path.to_string_lossy().to_string(),
                            error,
                        });
                        false
                    }
                }
            } else {
                false
            };

            // One ORF run gives the ORF output and the exon span for the intron.
            let orf_run = (runtime_recipe.enable_orf
                && (config.export_orf_alignments || export_introns))
                .then(|| apply_recipe(raw_align, &runtime_recipe, total_dataset_taxa));

            let (orf_accepted, orf_exported, orf_taxa, orf_length) =
                if let Some((orf_alignment, orf_diff)) =
                    orf_run.as_ref().filter(|_| config.export_orf_alignments)
                {
                    let accepted = orf_diff.orf_evaluated
                        && orf_diff.orf_candidate_found
                        && orf_diff.found_valid_orf
                        && !orf_alignment.sequences.is_empty()
                        && orf_alignment.num_taxa > 0
                        && orf_alignment.length > 0;
                    let orf_path = output_path(&orf_dir, &orf_alignment, config.output_format);
                    let exported = if accepted {
                        match write_alignment_atomic(
                            &orf_path,
                            &orf_alignment.taxa,
                            &orf_alignment.sequences,
                            config.output_format,
                        ) {
                            Ok(()) => true,
                            Err(error) => {
                                record_errors.push(BatchExportError {
                                    input_path: raw_align.file_path.clone(),
                                    output_path: orf_path.to_string_lossy().to_string(),
                                    error,
                                });
                                false
                            }
                        }
                    } else {
                        false
                    };
                    (
                        accepted,
                        exported,
                        orf_alignment.num_taxa,
                        orf_alignment.length,
                    )
                } else {
                    (false, false, 0, 0)
                };

            let mut intron_exported = false;
            // No span means no reference match, or a fallback ORF that can
            // include intron columns; the locus then has no intron file.
            let exon_span = orf_run.as_ref().and_then(|(_, orf_diff)| {
                Some((orf_diff.reference_exon_start?, orf_diff.reference_exon_end?))
            });
            if export_introns {
                if let Some((exon_start, exon_end)) = exon_span {
                    let intron_alignment = intron_alignment_outside(raw_align, exon_start, exon_end);
                    if intron_alignment.length > 0 {
                        let (filtered_intron, intron_diff) =
                            apply_recipe(&intron_alignment, &catalog_recipe, total_dataset_taxa);
                        let should_export_intron = !config.only_passing || intron_diff.pass;
                        if should_export_intron
                            && !filtered_intron.sequences.is_empty()
                            && filtered_intron.length > 0
                        {
                            let intron_path = output_path(
                                &intron_dir,
                                &filtered_intron,
                                config.output_format,
                            );
                            intron_exported = match write_alignment_atomic(
                                &intron_path,
                                &filtered_intron.taxa,
                                &filtered_intron.sequences,
                                config.output_format,
                            ) {
                                Ok(()) => true,
                                Err(error) => {
                                    record_errors.push(BatchExportError {
                                        input_path: raw_align.file_path.clone(),
                                        output_path: intron_path.to_string_lossy().to_string(),
                                        error,
                                    });
                                    false
                                }
                            };
                        }
                    }
                }
            }

            ExportRecord {
                id: catalog_alignment.id,
                catalog_pass: catalog_diff.pass,
                general_exported,
                orf_accepted,
                orf_exported,
                old_taxa: catalog_diff.old_taxa_count,
                catalog_taxa: catalog_diff.new_taxa_count,
                orf_taxa,
                old_length: catalog_diff.old_length,
                catalog_length: catalog_diff.new_length,
                orf_length,
                old_gap_percent: catalog_diff.old_gap_percent,
                catalog_gap_percent: catalog_diff.new_gap_percent,
                intron_exported,
                errors: record_errors,
            }
        })
        .collect();

    let total_processed = config.input_paths.len();
    let total_exported = results.iter().filter(|result| result.general_exported).count();
    let total_discarded = results.iter().filter(|result| !result.catalog_pass).count();
    let total_orfs_exported = results.iter().filter(|result| result.orf_exported).count();
    let total_introns_exported = results.iter().filter(|result| result.intron_exported).count();
    errors.extend(results.iter().flat_map(|result| result.errors.iter().cloned()));

    // Write summary CSV
    let summary_csv_path = if config.save_summary_csv {
        let csv_path = out_dir.join("alignment-trimming_summary.csv");
        {
            let mut writer = csv::Writer::from_path(&csv_path)
                .map_err(|error| format!("Failed to create export summary CSV: {error}"))?;
            writer.write_record([
                "Alignment", "CatalogPass", "GeneralExported", "ORFAccepted",
                "ORFExported", "startSamples", "catalogSamples", "orfSamples",
                "startLength", "catalogLength", "orfLength", "startPerGaps",
                "catalogPerGaps", "intronExported",
            ]).map_err(|error| format!("Failed to write export summary CSV: {error}"))?;
            for result in &results {
                writer.serialize((
                    &result.id, result.catalog_pass, result.general_exported,
                    result.orf_accepted, result.orf_exported, result.old_taxa,
                    result.catalog_taxa, result.orf_taxa, result.old_length,
                    result.catalog_length, result.orf_length,
                    format!("{:.2}", result.old_gap_percent),
                    format!("{:.2}", result.catalog_gap_percent), result.intron_exported,
                )).map_err(|error| format!("Failed to write export summary CSV: {error}"))?;
            }
            writer.flush()
                .map_err(|error| format!("Failed to finish export summary CSV: {error}"))?;
            Some(csv_path.to_string_lossy().to_string())
        }
    } else {
        None
    };

    // Write recipe JSON, and the reference sequences that it names
    let (recipe_json_path, reference_fasta_path) = if config.save_recipe_json {
        let reference_path = write_reference_fasta(out_dir, &runtime_recipe)?;
        let record = RecipeRecord {
            recipe: &runtime_recipe,
            orf_reference_file: reference_path.as_ref().map(|_| REFERENCE_FASTA_NAME),
        };
        let recipe_path = out_dir.join("recipe.json");
        let json = serde_json::to_string_pretty(&record)
            .map_err(|error| format!("Failed to serialize export recipe: {error}"))?;
        fs::write(&recipe_path, json)
            .map_err(|error| format!("Failed to write export recipe: {error}"))?;
        (
            Some(recipe_path.to_string_lossy().to_string()),
            reference_path.map(|path| path.to_string_lossy().to_string()),
        )
    } else {
        (None, None)
    };

    Ok(BatchExportResult {
        total_processed,
        total_exported,
        total_discarded,
        total_orfs_exported,
        alignment_directory_path: config
            .export_general_alignments
            .then(|| alignment_dir.to_string_lossy().to_string()),
        orf_directory_path: (config.export_orf_alignments && runtime_recipe.enable_orf)
            .then(|| orf_dir.to_string_lossy().to_string()),
        summary_csv_path,
        recipe_json_path,
        reference_fasta_path,
        total_introns_exported,
        intron_directory_path: (total_introns_exported > 0)
            .then(|| intron_dir.to_string_lossy().to_string()),
        total_failed: errors.len(),
        errors,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parsers::write_alignment;

    #[test]
    fn test_batch_export_folder_name_defaults() {
        let config: BatchExportConfig = serde_json::from_value(serde_json::json!({
            "input_paths": [],
            "output_directory": "trimmed_alignments",
            "output_format": "fasta",
            "only_passing": true,
            "save_recipe_json": false,
            "save_summary_csv": false
        }))
        .unwrap();

        assert_eq!(config.general_alignment_directory_name, "all_alignments");
        assert_eq!(config.orf_alignment_directory_name, "orf_alignments");
        assert_eq!(config.intron_directory_name, "intron_alignments");
    }

    #[test]
    fn reports_missing_inputs() {
        let output = std::env::temp_dir().join("alignmentforge_missing_export_input");
        let _ = fs::remove_dir_all(&output);
        let config = BatchExportConfig {
            input_paths: vec![output.join("missing.fa").to_string_lossy().to_string()],
            output_directory: output.to_string_lossy().to_string(),
            general_alignment_directory_name: "all_alignments".to_string(),
            orf_alignment_directory_name: "orf_alignments".to_string(),
            intron_directory_name: "intron_alignments".to_string(),
            output_format: AlignmentFormat::Fasta,
            only_passing: false,
            export_general_alignments: true,
            export_orf_alignments: false,
            save_recipe_json: false,
            save_summary_csv: false,
            export_introns: false,
        };

        let result = execute_batch_export(&config, &TrimmingRecipe::default(), None).unwrap();
        assert_eq!(result.total_exported, 0);
        assert_eq!(result.total_failed, 1);
        assert!(result.errors[0].error.contains("Failed to"));
        let _ = fs::remove_dir_all(output);
    }

    #[test]
    fn rejects_duplicate_output_names_before_writing() {
        let root = std::env::temp_dir().join("alignmentforge_duplicate_export_names");
        let _ = fs::remove_dir_all(&root);
        let first_dir = root.join("first");
        let second_dir = root.join("second");
        fs::create_dir_all(&first_dir).unwrap();
        fs::create_dir_all(&second_dir).unwrap();
        let first = first_dir.join("same.fa");
        let second = second_dir.join("SAME.fa");
        fs::write(&first, ">a\nAAAA\n").unwrap();
        fs::write(&second, ">a\nCCCC\n").unwrap();
        let config = BatchExportConfig {
            input_paths: vec![
                first.to_string_lossy().to_string(),
                second.to_string_lossy().to_string(),
            ],
            output_directory: root.join("output").to_string_lossy().to_string(),
            general_alignment_directory_name: "all_alignments".to_string(),
            orf_alignment_directory_name: "orf_alignments".to_string(),
            intron_directory_name: "intron_alignments".to_string(),
            output_format: AlignmentFormat::Fasta,
            only_passing: false,
            export_general_alignments: true,
            export_orf_alignments: false,
            save_recipe_json: false,
            save_summary_csv: false,
            export_introns: false,
        };

        let error = execute_batch_export(&config, &TrimmingRecipe::default(), None).unwrap_err();
        assert!(error.contains("same output path"));
        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn batch_export_writes_alignment_summary_and_recipe() {
        let root = std::env::temp_dir().join("alignmentforge_batch_export_test");
        let output = root.join("output");
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let input = root.join("locus1.fa");
        fs::write(&input, ">a\nACGTACGT\n>b\nACGTACGA\n>c\nACGTACGC\n").unwrap();

        let config = BatchExportConfig {
            input_paths: vec![input.to_string_lossy().to_string()],
            output_directory: output.to_string_lossy().to_string(),
            general_alignment_directory_name: "all_alignments".to_string(),
            orf_alignment_directory_name: "orf_alignments".to_string(),
            intron_directory_name: "intron_alignments".to_string(),
            output_format: AlignmentFormat::Phylip,
            only_passing: false,
            export_general_alignments: true,
            export_orf_alignments: false,
            save_recipe_json: true,
            save_summary_csv: true,
            export_introns: false,
        };
        let recipe = TrimmingRecipe {
            trim_similarity: false,
            trim_external: false,
            trim_coverage: false,
            assess_alignment: false,
            ..TrimmingRecipe::default()
        };

        let result = execute_batch_export(&config, &recipe, None).unwrap();

        assert_eq!(result.total_processed, 1);
        assert_eq!(result.total_exported, 1);
        assert_eq!(result.total_failed, 0);
        let exported = parse_alignment(output.join("all_alignments/locus1.phy")).unwrap();
        assert_eq!(exported.taxa, vec!["a", "b", "c"]);
        assert_eq!(exported.sequences, vec!["ACGTACGT", "ACGTACGA", "ACGTACGC"]);

        let summary = fs::read_to_string(output.join("alignment-trimming_summary.csv")).unwrap();
        let rows: Vec<&str> = summary.lines().collect();
        assert_eq!(rows.len(), 2, "{summary}");
        assert!(rows[0].starts_with("Alignment,CatalogPass,GeneralExported,"), "{summary}");
        assert!(rows[1].starts_with("locus1,true,true,"), "{summary}");

        let saved: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(output.join("recipe.json")).unwrap()).unwrap();
        assert!(saved.get("orf_reference_file").is_none());
        let saved_recipe: TrimmingRecipe = serde_json::from_value(saved).unwrap();
        assert_eq!(saved_recipe.name, recipe.name);
        assert!(!saved_recipe.trim_similarity);
        assert!(result.reference_fasta_path.is_none());

        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn recipe_json_names_a_fasta_file_with_the_reference_sequences() {
        let root = std::env::temp_dir().join("alignmentforge_recipe_reference_test");
        let output = root.join("output");
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let input = root.join("exon1.fa");
        fs::write(&input, ">a\nATGAAATAA\n>b\nATGAAATAA\n").unwrap();

        let config = BatchExportConfig {
            input_paths: vec![input.to_string_lossy().to_string()],
            output_directory: output.to_string_lossy().to_string(),
            general_alignment_directory_name: "all_alignments".to_string(),
            orf_alignment_directory_name: "orf_alignments".to_string(),
            intron_directory_name: "intron_alignments".to_string(),
            output_format: AlignmentFormat::Fasta,
            only_passing: false,
            export_general_alignments: true,
            export_orf_alignments: false,
            save_recipe_json: true,
            save_summary_csv: false,
            export_introns: false,
        };
        let recipe = TrimmingRecipe {
            orf_use_references: true,
            orf_search_mode: OrfSearchMode::ReferenceGuided,
            orf_reference_sequences: HashMap::from([
                ("exon2".to_string(), "ATGCCCCCCTAA".to_string()),
                ("exon1".to_string(), "ATGAAATAA".to_string()),
            ]),
            ..TrimmingRecipe::default()
        };

        let result = execute_batch_export(&config, &recipe, None).unwrap();

        // Sorted by locus, in the form that the ORF sidebar loads.
        let reference_path = result.reference_fasta_path.expect("no reference file was written");
        assert_eq!(
            fs::read_to_string(&reference_path).unwrap(),
            ">exon1\nATGAAATAA\n>exon2\nATGCCCCCCTAA\n"
        );
        let saved: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(output.join("recipe.json")).unwrap()).unwrap();
        assert_eq!(saved["orf_reference_file"], "recipe_references.fasta");
        assert!(saved.get("orf_reference_sequences").is_none());
        assert_eq!(saved["orf_search_mode"], "referenceguided");

        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn test_batch_export_writes_reference_introns_separately() {
        let root = std::env::temp_dir().join("alignmentforge_intron_export_test");
        let input_dir = root.join("input");
        let output_dir = root.join("output");
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&input_dir).unwrap();
        let input_path = input_dir.join("exon_export.fa");
        let exon = "ATGGCTGCTGCTGCTGCTGCTGCTGCTTAA";
        let taxa = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let sequences = vec![format!("CCCCCC{exon}GGGGGG"); 3];
        write_alignment(&input_path, &taxa, &sequences, AlignmentFormat::Fasta).unwrap();

        let mut recipe = TrimmingRecipe::default();
        recipe.enable_orf = true;
        recipe.orf_use_references = true;
        recipe.orf_search_mode = OrfSearchMode::ReferenceGuided;
        recipe
            .orf_reference_sequences
            .insert("exon_export".to_string(), exon.to_string());
        recipe.trim_similarity = false;
        recipe.trim_hmm = false;
        recipe.trim_external = false;
        recipe.trim_coverage = false;
        recipe.assess_alignment = false;
        let config = BatchExportConfig {
            input_paths: vec![input_path.to_string_lossy().to_string()],
            output_directory: output_dir.to_string_lossy().to_string(),
            general_alignment_directory_name: "custom_general".to_string(),
            orf_alignment_directory_name: "custom_orfs".to_string(),
            intron_directory_name: "custom_introns".to_string(),
            output_format: AlignmentFormat::Fasta,
            only_passing: false,
            export_general_alignments: true,
            export_orf_alignments: true,
            save_recipe_json: false,
            save_summary_csv: false,
            export_introns: true,
        };

        let result = execute_batch_export(&config, &recipe, None).unwrap();
        assert_eq!(result.total_introns_exported, 1);
        let intron_path = output_dir.join("custom_introns/exon_export_intron.fa");
        let intron = fs::read_to_string(&intron_path).unwrap();
        assert!(intron.contains("CCCCCCGGGGGG"), "{intron}");

        // Other ORF modes do not use the reference span, so they write no introns.
        let _ = fs::remove_dir_all(&output_dir);
        recipe.orf_search_mode = OrfSearchMode::ContinuousCds;
        let result = execute_batch_export(&config, &recipe, None).unwrap();
        assert_eq!(result.total_introns_exported, 0);

        let _ = fs::remove_dir_all(root);
    }

    #[test]
    fn test_general_and_orf_exports_use_independent_statuses() {
        use crate::algorithms::orf::StopCodonAction;

        let root = std::env::temp_dir().join("alignmentforge_independent_export_test");
        let input_dir = root.join("input");
        let output_dir = root.join("output");
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&input_dir).unwrap();
        let input_path = input_dir.join("exon_independent.fa");
        write_alignment(
            &input_path,
            &["A".into(), "B".into(), "C".into(), "Pseudogene".into()],
            &[
                "ATGAAAGGG".into(),
                "ATGAAAGGG".into(),
                "ATGAAAGGG".into(),
                "ATGTAAGGG".into(),
            ],
            AlignmentFormat::Fasta,
        )
        .unwrap();

        let mut recipe = TrimmingRecipe::default();
        recipe.enable_orf = true;
        recipe.exclude_uce = false;
        recipe.stop_codon_action = StopCodonAction::RemoveSample;
        recipe.trim_similarity = false;
        recipe.trim_hmm = false;
        recipe.trim_external = false;
        recipe.trim_coverage = false;
        recipe.min_taxa = 4;
        recipe.min_taxa_occupancy_percent = 0.0;
        recipe.min_length = 0;

        let config = BatchExportConfig {
            input_paths: vec![input_path.to_string_lossy().to_string()],
            output_directory: output_dir.to_string_lossy().to_string(),
            general_alignment_directory_name: "all_alignments".to_string(),
            orf_alignment_directory_name: "orf_alignments".to_string(),
            intron_directory_name: "intron_alignments".to_string(),
            output_format: AlignmentFormat::Fasta,
            only_passing: true,
            export_general_alignments: true,
            export_orf_alignments: true,
            save_recipe_json: false,
            save_summary_csv: true,
            export_introns: false,
        };

        let result = execute_batch_export(&config, &recipe, None).unwrap();
        assert_eq!(result.total_exported, 1);
        assert_eq!(result.total_orfs_exported, 1);

        let general =
            parse_alignment(output_dir.join("all_alignments/exon_independent.fa")).unwrap();
        let orf = parse_alignment(output_dir.join("orf_alignments/exon_independent.fa")).unwrap();
        assert_eq!(general.num_taxa, 4);
        assert_eq!(orf.num_taxa, 3);

        // The inverse must also hold: Catalog QC can reject the general
        // alignment without suppressing an independently accepted ORF.
        let catalog_fail_output = root.join("catalog_fail_output");
        let mut catalog_fail_recipe = recipe.clone();
        catalog_fail_recipe.min_length = 100;
        let mut catalog_fail_config = config.clone();
        catalog_fail_config.output_directory =
            catalog_fail_output.to_string_lossy().to_string();
        let catalog_fail_result =
            execute_batch_export(&catalog_fail_config, &catalog_fail_recipe, None).unwrap();
        assert_eq!(catalog_fail_result.total_exported, 0);
        assert_eq!(catalog_fail_result.total_discarded, 1);
        assert_eq!(catalog_fail_result.total_orfs_exported, 1);
        assert!(!catalog_fail_output
            .join("all_alignments/exon_independent.fa")
            .exists());
        assert!(catalog_fail_output
            .join("orf_alignments/exon_independent.fa")
            .exists());

        let _ = fs::remove_dir_all(root);
    }
}
