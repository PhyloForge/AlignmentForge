use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use serde::{Deserialize, Serialize};

use crate::models::{Alignment, AlignmentFormat};
use crate::parsers::{parse_alignment, write_alignment};
use crate::pipeline::catalog::recipe_with_dataset_sample_filter;
use crate::pipeline::engine::apply_recipe;
use crate::pipeline::recipe::TrimmingRecipe;
use crate::export::concatenate::{LocusPartition};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroupedConcatenateConfig {
    pub input_paths: Vec<String>,
    pub output_directory: String,
    pub gene_mapping_csv_path: String,
    pub output_format: AlignmentFormat,
    pub only_passing: bool,
    pub write_raxml_partitions: bool,
    pub write_nexus_partitions: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroupedConcatenateResult {
    pub total_genes: usize,
    pub total_exons_processed: usize,
    pub output_directory: String,
}

pub fn concatenate_alignments_by_gene(
    config: &GroupedConcatenateConfig,
    recipe: &TrimmingRecipe,
) -> Result<GroupedConcatenateResult, String> {
    // 1. Parse the metadata file (CSV, TSV, or TXT)
    let file = File::open(&config.gene_mapping_csv_path)
        .map_err(|e| format!("Could not open gene mapping file: {}", e))?;
    let reader = BufReader::new(file);

    // exon_name -> gene_name
    let mut exon_to_gene: BTreeMap<String, String> = BTreeMap::new();
    
    for line in reader.lines() {
        let line = line.unwrap_or_default();
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        // Try to split by tab first, if no tab, try comma
        let separator = if trimmed.contains('\t') { '\t' } else { ',' };
        let parts: Vec<&str> = trimmed.split(separator).collect();
        
        if parts.len() >= 2 {
            let exon = parts[0].trim().to_string();
            let gene = parts[1].trim().to_string();
            // Skip headers if they happen to be obvious, or just include them (they won't match a file)
            exon_to_gene.insert(exon, gene);
        }
    }

    if exon_to_gene.is_empty() {
        return Err("No valid mapping found in the metadata document.".to_string());
    }

    // 2. Group input paths by gene
    let mut gene_to_paths: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for path_str in &config.input_paths {
        let path = Path::new(path_str);
        if let Some(stem) = path.file_stem().and_then(|s| s.to_str()) {
            if let Some(gene) = exon_to_gene.get(stem) {
                gene_to_paths.entry(gene.clone()).or_default().push(path_str.clone());
            }
        }
    }

    if gene_to_paths.is_empty() {
        return Err("No input files matched the exons in the metadata document.".to_string());
    }

    // Ensure output directory exists
    std::fs::create_dir_all(&config.output_directory)
        .map_err(|e| format!("Failed to create output directory: {}", e))?;

    let mut total_genes = 0;
    let mut total_exons_processed = 0;

    // 3. Process each gene group
    for (gene, paths) in &gene_to_paths {
        let raw_alignments: Vec<Alignment> = paths
            .iter()
            .filter_map(|p| parse_alignment(Path::new(p)).ok())
            .collect();
            
        if raw_alignments.is_empty() {
            continue;
        }

        let runtime_recipe = if recipe.excluded_taxa.is_empty() {
            recipe_with_dataset_sample_filter(recipe, &raw_alignments)
        } else {
            recipe.clone()
        };

        let mut passing_alignments = Vec::new();
        let mut all_taxa_set = BTreeSet::new();

        for raw in &raw_alignments {
            let (transformed, diff) = apply_recipe(raw, &runtime_recipe, 0);
            if (!config.only_passing || diff.pass) && transformed.length > 0 && !transformed.taxa.is_empty() {
                for taxon in &transformed.taxa {
                    all_taxa_set.insert(taxon.clone());
                }
                passing_alignments.push(transformed);
            }
        }

        if passing_alignments.is_empty() {
            continue;
        }

        let all_taxa: Vec<String> = all_taxa_set.into_iter().collect();
        let mut concatenated_seqs: BTreeMap<String, String> = BTreeMap::new();
        for taxon in &all_taxa {
            concatenated_seqs.insert(taxon.clone(), String::new());
        }

        let mut partitions = Vec::new();
        let mut current_offset = 0usize;

        for align in &passing_alignments {
            let locus_len = align.length;
            let start_1based = current_offset + 1;
            let end_1based = current_offset + locus_len;

            partitions.push(LocusPartition {
                name: align.id.clone(),
                start: start_1based,
                end: end_1based,
                length: locus_len,
            });

            let taxon_seq_map: BTreeMap<&str, &str> = align
                .taxa
                .iter()
                .zip(align.sequences.iter())
                .map(|(t, s)| (t.as_str(), s.as_str()))
                .collect();

            let gap_pad = "-".repeat(locus_len);
            for taxon in &all_taxa {
                let seq_chunk = taxon_seq_map.get(taxon.as_str()).copied().unwrap_or(&gap_pad);
                concatenated_seqs.get_mut(taxon).unwrap().push_str(seq_chunk);
            }
            current_offset += locus_len;
            total_exons_processed += 1;
        }

        let final_taxa: Vec<String> = concatenated_seqs.keys().cloned().collect();
        let final_seqs: Vec<String> = concatenated_seqs.values().cloned().collect();
        let out_prefix = Path::new(&config.output_directory).join(gene);
        let out_prefix_str = out_prefix.to_string_lossy();
        let supermatrix_path = format!("{}.{}", out_prefix_str, config.output_format.extension());

        if let Err(e) = write_alignment(
            &supermatrix_path,
            &final_taxa,
            &final_seqs,
            config.output_format,
        ) {
            eprintln!("Failed to write gene alignment {}: {}", gene, e);
        }

        if config.write_raxml_partitions {
            let raxml_path = format!("{}_partitions.txt", out_prefix_str);
            if let Ok(mut file) = File::create(&raxml_path) {
                use std::io::Write;
                for part in &partitions {
                    let _ = writeln!(file, "DNA, {} = {}-{}", part.name, part.start, part.end);
                }
            }
        }

        if config.write_nexus_partitions {
            let nex_path = format!("{}_partitions.nex", out_prefix_str);
            if let Ok(mut file) = File::create(&nex_path) {
                use std::io::Write;
                let _ = writeln!(file, "#NEXUS\nBEGIN SETS;");
                for part in &partitions {
                    let _ = writeln!(file, "  CHARSET {} = {}-{};", part.name, part.start, part.end);
                }
                let _ = writeln!(file, "END;");
            }
        }

        total_genes += 1;
    }

    Ok(GroupedConcatenateResult {
        total_genes,
        total_exons_processed,
        output_directory: config.output_directory.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::io::Write;
    use crate::models::AlignmentFormat;

    #[test]
    fn test_concatenate_alignments_by_gene() {
        let temp_dir = std::env::temp_dir().join("test_group_concat_dir");
        let _ = fs::create_dir_all(&temp_dir);

        // 1. Create a mock mapping file
        let mapping_path = temp_dir.join("mapping.csv");
        let mut map_file = fs::File::create(&mapping_path).unwrap();
        writeln!(map_file, "Exon1,GeneA").unwrap();
        writeln!(map_file, "Exon2,GeneA").unwrap();
        writeln!(map_file, "Exon3,GeneB").unwrap();

        // 2. Mock input alignments
        let test_input1 = "../test_data/uce-1002.phy";
        let test_input2 = "../test_data/uce-1003.fa";
        
        // We will symlink or copy these into our temp_dir to give them names "Exon1", "Exon2", etc.
        let exon1_path = temp_dir.join("Exon1.phy");
        let exon2_path = temp_dir.join("Exon2.fa");
        let exon3_path = temp_dir.join("Exon3.phy");

        if Path::new(test_input1).exists() && Path::new(test_input2).exists() {
            fs::copy(test_input1, &exon1_path).unwrap();
            fs::copy(test_input2, &exon2_path).unwrap();
            fs::copy(test_input1, &exon3_path).unwrap();

            let config = GroupedConcatenateConfig {
                input_paths: vec![
                    exon1_path.to_string_lossy().to_string(),
                    exon2_path.to_string_lossy().to_string(),
                    exon3_path.to_string_lossy().to_string(),
                ],
                output_directory: temp_dir.to_string_lossy().to_string(),
                gene_mapping_csv_path: mapping_path.to_string_lossy().to_string(),
                output_format: AlignmentFormat::Phylip,
                only_passing: false,
                write_raxml_partitions: true,
                write_nexus_partitions: false,
            };

            let mut recipe = TrimmingRecipe::default();
            recipe.trim_coverage = false;

            let res = concatenate_alignments_by_gene(&config, &recipe).unwrap();
            
            assert_eq!(res.total_genes, 2);
            assert_eq!(res.total_exons_processed, 3);
            assert!(temp_dir.join("GeneA.phy").exists());
            assert!(temp_dir.join("GeneA_partitions.txt").exists());
            assert!(temp_dir.join("GeneB.phy").exists());
            assert!(temp_dir.join("GeneB_partitions.txt").exists());
        }

        let _ = fs::remove_dir_all(temp_dir);
    }
}
