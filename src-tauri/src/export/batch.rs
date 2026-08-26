use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::models::{Alignment, AlignmentFormat};
use crate::parsers::{parse_alignment, write_alignment};
use crate::pipeline::catalog::{
    recipe_with_dataset_sample_filter, recipe_without_orf_analysis,
};
use crate::pipeline::engine::apply_recipe;
use crate::pipeline::recipe::TrimmingRecipe;
use crate::algorithms::reference::intron_alignment_from_reference;

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
    pub total_introns_exported: usize,
    pub intron_directory_path: Option<String>,
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

pub fn execute_batch_export(
    config: &BatchExportConfig,
    recipe: &TrimmingRecipe,
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
    let mut active_directory_names = Vec::new();
    if config.export_general_alignments {
        active_directory_names.push(general_directory_name.to_lowercase());
    }
    if config.export_orf_alignments && recipe.enable_orf {
        active_directory_names.push(orf_directory_name.to_lowercase());
    }
    if config.export_introns && recipe.enable_orf && recipe.orf_use_references {
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
    if !out_dir.exists() {
        fs::create_dir_all(out_dir)
            .map_err(|e| format!("Failed to create output directory: {}", e))?;
    }

    let raw_alignments: Vec<Alignment> = config
        .input_paths
        .par_iter()
        .filter_map(|input_path| parse_alignment(Path::new(input_path)).ok())
        .collect();
    let runtime_recipe = if recipe.excluded_taxa.is_empty() {
        recipe_with_dataset_sample_filter(recipe, &raw_alignments)
    } else {
        recipe.clone()
    };
    let alignment_dir = out_dir.join(general_directory_name);
    let orf_dir = out_dir.join(orf_directory_name);
    let intron_dir = out_dir.join(intron_directory_name);
    if config.export_general_alignments {
        fs::create_dir_all(&alignment_dir)
            .map_err(|e| format!("Failed to create general alignment output directory: {e}"))?;
    }
    if config.export_orf_alignments && runtime_recipe.enable_orf {
        fs::create_dir_all(&orf_dir)
            .map_err(|e| format!("Failed to create ORF output directory: {e}"))?;
    }
    if config.export_introns && runtime_recipe.enable_orf && runtime_recipe.orf_use_references {
        fs::create_dir_all(&intron_dir)
            .map_err(|e| format!("Failed to create intron output directory: {e}"))?;
    }

    let catalog_recipe = recipe_without_orf_analysis(&runtime_recipe);

    // Process all files in parallel using one dataset-wide sample exclusion
    // list, but keep general alignment and ORF outcomes independent.
    let results: Vec<ExportRecord> = raw_alignments
        .par_iter()
        .map(|raw_align| {
            let (catalog_alignment, catalog_diff) = apply_recipe(raw_align, &catalog_recipe, 0);

            let should_export_general = config.export_general_alignments
                && (!config.only_passing || catalog_diff.pass)
                && !catalog_alignment.sequences.is_empty()
                && catalog_alignment.length > 0;
            let general_exported = should_export_general
                && write_alignment(
                    alignment_dir.join(format!(
                        "{}.{}",
                        catalog_alignment.id,
                        config.output_format.extension()
                    )),
                    &catalog_alignment.taxa,
                    &catalog_alignment.sequences,
                    config.output_format,
                )
                .is_ok();

            let (orf_accepted, orf_exported, orf_taxa, orf_length) =
                if config.export_orf_alignments && runtime_recipe.enable_orf {
                    let (orf_alignment, orf_diff) = apply_recipe(raw_align, &runtime_recipe, 0);
                    let accepted = orf_diff.orf_evaluated
                        && orf_diff.orf_candidate_found
                        && orf_diff.found_valid_orf
                        && !orf_alignment.sequences.is_empty()
                        && orf_alignment.num_taxa > 0
                        && orf_alignment.length > 0;
                    let exported = accepted
                        && write_alignment(
                            orf_dir.join(format!(
                                "{}.{}",
                                orf_alignment.id,
                                config.output_format.extension()
                            )),
                            &orf_alignment.taxa,
                            &orf_alignment.sequences,
                            config.output_format,
                        )
                        .is_ok();
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
            if config.export_introns
                && runtime_recipe.enable_orf
                && runtime_recipe.orf_use_references
            {
                if let Some(reference) = runtime_recipe.orf_reference_sequences.get(&raw_align.id) {
                    if let Some((intron_alignment, _)) =
                        intron_alignment_from_reference(raw_align, reference)
                    {
                        if intron_alignment.length > 0 {
                            let mut intron_recipe = catalog_recipe.clone();
                            intron_recipe.orf_reference_sequences.clear();
                            let (filtered_intron, intron_diff) =
                                apply_recipe(&intron_alignment, &intron_recipe, 0);
                            let should_export_intron = !config.only_passing || intron_diff.pass;
                            if should_export_intron
                                && !filtered_intron.sequences.is_empty()
                                && filtered_intron.length > 0
                            {
                                intron_exported = write_alignment(
                                    intron_dir.join(format!(
                                        "{}.{}",
                                        filtered_intron.id,
                                        config.output_format.extension()
                                    )),
                                    &filtered_intron.taxa,
                                    &filtered_intron.sequences,
                                    config.output_format,
                                )
                                .is_ok();
                            }
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
            }
        })
        .collect();

    let total_processed = results.len();
    let total_exported = results.iter().filter(|result| result.general_exported).count();
    let total_discarded = results.iter().filter(|result| !result.catalog_pass).count();
    let total_orfs_exported = results.iter().filter(|result| result.orf_exported).count();
    let total_introns_exported = results.iter().filter(|result| result.intron_exported).count();

    // Write summary CSV
    let summary_csv_path = if config.save_summary_csv {
        let csv_path = out_dir.join("alignment-trimming_summary.csv");
        if let Ok(mut file) = File::create(&csv_path) {
            let _ = writeln!(
                file,
                "Alignment,CatalogPass,GeneralExported,ORFAccepted,ORFExported,startSamples,catalogSamples,orfSamples,startLength,catalogLength,orfLength,startPerGaps,catalogPerGaps,intronExported"
            );
            for result in &results {
                let _ = writeln!(
                    file,
                    "{},{},{},{},{},{},{},{},{},{},{},{:.2},{:.2},{}",
                    result.id,
                    result.catalog_pass,
                    result.general_exported,
                    result.orf_accepted,
                    result.orf_exported,
                    result.old_taxa,
                    result.catalog_taxa,
                    result.orf_taxa,
                    result.old_length,
                    result.catalog_length,
                    result.orf_length,
                    result.old_gap_percent,
                    result.catalog_gap_percent,
                    result.intron_exported
                );
            }
            Some(csv_path.to_string_lossy().to_string())
        } else {
            None
        }
    } else {
        None
    };

    // Write recipe JSON
    let recipe_json_path = if config.save_recipe_json {
        let recipe_path = out_dir.join("recipe.json");
        if let Ok(json_str) = serde_json::to_string_pretty(recipe) {
            let _ = fs::write(&recipe_path, json_str);
            Some(recipe_path.to_string_lossy().to_string())
        } else {
            None
        }
    } else {
        None
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
        total_introns_exported,
        intron_directory_path: (total_introns_exported > 0)
            .then(|| intron_dir.to_string_lossy().to_string()),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn test_batch_export() {
        let temp_dir = std::env::temp_dir().join("test_batch_export_dir");
        let _ = fs::remove_dir_all(&temp_dir);

        let test_input = "../test_data/uce-1001.phy";
        if std::path::Path::new(test_input).exists() {
            let config = BatchExportConfig {
                input_paths: vec![test_input.to_string()],
                output_directory: temp_dir.to_string_lossy().to_string(),
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
            let recipe = TrimmingRecipe::default();
            let res = execute_batch_export(&config, &recipe).unwrap();

            assert_eq!(res.total_processed, 1);
            assert!(temp_dir.join("all_alignments/uce-1001.phy").exists());
            assert!(temp_dir.join("alignment-trimming_summary.csv").exists());
            assert!(temp_dir.join("recipe.json").exists());
        }

        let _ = fs::remove_dir_all(temp_dir);
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

        let result = execute_batch_export(&config, &recipe).unwrap();
        assert_eq!(result.total_introns_exported, 1);
        assert!(output_dir
            .join("custom_introns/exon_export_intron.fa")
            .exists());

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

        let result = execute_batch_export(&config, &recipe).unwrap();
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
            execute_batch_export(&catalog_fail_config, &catalog_fail_recipe).unwrap();
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
