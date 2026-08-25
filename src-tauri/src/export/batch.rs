use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::models::{Alignment, AlignmentFormat};
use crate::parsers::{parse_alignment, write_alignment};
use crate::pipeline::catalog::recipe_with_dataset_sample_filter;
use crate::pipeline::engine::apply_recipe;
use crate::pipeline::recipe::TrimmingRecipe;
use crate::algorithms::reference::intron_alignment_from_reference;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BatchExportConfig {
    pub input_paths: Vec<String>,
    pub output_directory: String,
    pub output_format: AlignmentFormat,
    pub only_passing: bool,
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
    pub summary_csv_path: Option<String>,
    pub recipe_json_path: Option<String>,
    pub total_introns_exported: usize,
    pub intron_directory_path: Option<String>,
}

pub fn execute_batch_export(
    config: &BatchExportConfig,
    recipe: &TrimmingRecipe,
) -> Result<BatchExportResult, String> {
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
    let intron_dir = out_dir.join("introns");
    if config.export_introns && runtime_recipe.enable_orf && runtime_recipe.orf_use_references {
        fs::create_dir_all(&intron_dir)
            .map_err(|e| format!("Failed to create intron output directory: {e}"))?;
    }

    // Process all files in parallel using one dataset-wide sample exclusion list.
    let results: Vec<(String, bool, usize, usize, usize, usize, f64, f64, bool)> = raw_alignments
        .par_iter()
        .map(|raw_align| {
            let (transformed, diff) = apply_recipe(raw_align, &runtime_recipe, 0);

            let should_export = !config.only_passing || diff.pass;
            if should_export && !transformed.sequences.is_empty() {
                let out_filename = format!("{}.{}", transformed.id, config.output_format.extension());
                let out_file_path = out_dir.join(out_filename);
                let _ = write_alignment(
                    &out_file_path,
                    &transformed.taxa,
                    &transformed.sequences,
                    config.output_format,
                );
            }

            let mut intron_exported = false;
            if config.export_introns && runtime_recipe.enable_orf && runtime_recipe.orf_use_references {
                if let Some(reference) = runtime_recipe.orf_reference_sequences.get(&raw_align.id) {
                    if let Some((intron_alignment, _)) =
                        intron_alignment_from_reference(raw_align, reference)
                    {
                        if intron_alignment.length > 0 {
                            let mut intron_recipe = runtime_recipe.clone();
                            intron_recipe.enable_orf = false;
                            intron_recipe.orf_use_references = false;
                            intron_recipe.orf_reference_sequences.clear();
                            let (filtered_intron, intron_diff) =
                                apply_recipe(&intron_alignment, &intron_recipe, 0);
                            let should_export_intron =
                                !config.only_passing || intron_diff.pass;
                            if should_export_intron
                                && !filtered_intron.sequences.is_empty()
                                && filtered_intron.length > 0
                            {
                                let intron_filename = format!(
                                    "{}.{}",
                                    filtered_intron.id,
                                    config.output_format.extension()
                                );
                                intron_exported = write_alignment(
                                    intron_dir.join(intron_filename),
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

            (
                transformed.id,
                diff.pass,
                diff.old_taxa_count,
                diff.new_taxa_count,
                diff.old_length,
                diff.new_length,
                diff.old_gap_percent,
                diff.new_gap_percent,
                intron_exported,
            )
        })
        .collect();

    let total_processed = results.len();
    let total_exported = results.iter().filter(|r| !config.only_passing || r.1).count();
    let total_discarded = total_processed - results.iter().filter(|r| r.1).count();
    let total_introns_exported = results.iter().filter(|result| result.8).count();

    // Write summary CSV
    let summary_csv_path = if config.save_summary_csv {
        let csv_path = out_dir.join("alignment-trimming_summary.csv");
        if let Ok(mut file) = File::create(&csv_path) {
            let _ = writeln!(
                file,
                "Alignment,Pass,startSamples,trimmedSamples,startLength,trimmedLength,startPerGaps,trimmedPerGaps,intronExported"
            );
            for (id, pass, old_taxa, new_taxa, old_len, new_len, old_gap, new_gap, intron_exported) in &results {
                let _ = writeln!(
                    file,
                    "{},{},{},{},{},{},{:.2},{:.2},{}",
                    id, pass, old_taxa, new_taxa, old_len, new_len, old_gap, new_gap,
                    intron_exported
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
    fn test_batch_export() {
        let temp_dir = std::env::temp_dir().join("test_batch_export_dir");
        let _ = fs::remove_dir_all(&temp_dir);

        let test_input = "../test_data/uce-1001.phy";
        if std::path::Path::new(test_input).exists() {
            let config = BatchExportConfig {
                input_paths: vec![test_input.to_string()],
                output_directory: temp_dir.to_string_lossy().to_string(),
                output_format: AlignmentFormat::Phylip,
                only_passing: false,
                save_recipe_json: true,
                save_summary_csv: true,
                export_introns: false,
            };
            let recipe = TrimmingRecipe::default();
            let res = execute_batch_export(&config, &recipe).unwrap();

            assert_eq!(res.total_processed, 1);
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
            output_format: AlignmentFormat::Fasta,
            only_passing: false,
            save_recipe_json: false,
            save_summary_csv: false,
            export_introns: true,
        };

        let result = execute_batch_export(&config, &recipe).unwrap();
        assert_eq!(result.total_introns_exported, 1);
        assert!(output_dir.join("introns/exon_export_intron.fa").exists());

        let _ = fs::remove_dir_all(root);
    }
}
