use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use serde::{Deserialize, Serialize};

use crate::models::{Alignment, AlignmentFormat};
use crate::parsers::{parse_alignment, write_alignment};
use crate::pipeline::catalog::recipe_with_dataset_sample_filter;
use crate::pipeline::engine::apply_recipe;
use crate::pipeline::recipe::TrimmingRecipe;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConcatenateConfig {
    pub input_paths: Vec<String>,
    pub output_file_prefix: String,
    pub output_format: AlignmentFormat,
    pub only_passing: bool,
    pub write_raxml_partitions: bool,
    pub write_nexus_partitions: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConcatenateResult {
    pub total_taxa: usize,
    pub total_length: usize,
    pub total_loci: usize,
    pub supermatrix_path: String,
    pub raxml_partition_path: Option<String>,
    pub nexus_partition_path: Option<String>,
}

#[derive(Debug, Clone)]
pub struct LocusPartition {
    pub name: String,
    pub start: usize, // 1-based inclusive
    pub end: usize,   // 1-based inclusive
    pub length: usize,
}

/// Joins the loci into one supermatrix in the given order. Taxa are sorted by
/// name, and a taxon that is absent from a locus gets gaps for that locus.
pub fn build_supermatrix(
    alignments: &[Alignment],
) -> (Vec<String>, Vec<String>, Vec<LocusPartition>) {
    // Map each taxon to its string buffer
    let mut concatenated_seqs: BTreeMap<String, String> = BTreeMap::new();
    for align in alignments {
        for taxon in &align.taxa {
            concatenated_seqs.entry(taxon.clone()).or_default();
        }
    }

    let mut partitions = Vec::new();
    let mut current_offset = 0usize;

    for align in alignments {
        let locus_len = align.length;
        partitions.push(LocusPartition {
            name: align.id.clone(),
            start: current_offset + 1,
            end: current_offset + locus_len,
            length: locus_len,
        });

        // Fast lookup map for this locus
        let taxon_seq_map: BTreeMap<&str, &str> = align
            .taxa
            .iter()
            .zip(align.sequences.iter())
            .map(|(t, s)| (t.as_str(), s.as_str()))
            .collect();

        let gap_pad = "-".repeat(locus_len);
        for (taxon, sequence) in concatenated_seqs.iter_mut() {
            let seq_chunk = taxon_seq_map.get(taxon.as_str()).copied().unwrap_or(&gap_pad);
            sequence.push_str(seq_chunk);
        }

        current_offset += locus_len;
    }

    let (taxa, sequences): (Vec<String>, Vec<String>) = concatenated_seqs.into_iter().unzip();
    (taxa, sequences, partitions)
}

/// Writes a RAxML-style partition file with one line per locus.
pub fn write_raxml_partitions(path: &str, partitions: &[LocusPartition]) -> Result<(), String> {
    let mut file = File::create(path)
        .map_err(|error| format!("Failed to create RAxML partition file '{path}': {error}"))?;
    for part in partitions {
        writeln!(file, "DNA, {} = {}-{}", part.name, part.start, part.end)
            .map_err(|error| format!("Failed to write RAxML partition file '{path}': {error}"))?;
    }
    Ok(())
}

/// Writes a NEXUS / IQ-TREE partition file with one CHARSET per locus.
pub fn write_nexus_partitions(path: &str, partitions: &[LocusPartition]) -> Result<(), String> {
    let mut file = File::create(path)
        .map_err(|error| format!("Failed to create NEXUS partition file '{path}': {error}"))?;
    writeln!(file, "#NEXUS\nBEGIN SETS;")
        .map_err(|error| format!("Failed to write NEXUS partition file '{path}': {error}"))?;
    for part in partitions {
        writeln!(file, "  CHARSET {} = {}-{};", part.name, part.start, part.end)
            .map_err(|error| format!("Failed to write NEXUS partition file '{path}': {error}"))?;
    }
    writeln!(file, "END;")
        .map_err(|error| format!("Failed to finish NEXUS partition file '{path}': {error}"))?;
    Ok(())
}

pub fn concatenate_alignments(
    config: &ConcatenateConfig,
    recipe: &TrimmingRecipe,
    resolved_dataset_taxa: Option<usize>,
) -> Result<ConcatenateResult, String> {
    let mut passing_alignments = Vec::new();

    let raw_alignments: Vec<Alignment> = config
        .input_paths
        .iter()
        .map(|path| {
            parse_alignment(Path::new(path))
                .map_err(|error| format!("Failed to read input '{path}': {error}"))
        })
        .collect::<Result<_, _>>()?;
    let runtime_recipe = if resolved_dataset_taxa.is_none() {
        recipe_with_dataset_sample_filter(recipe, &raw_alignments)
    } else {
        recipe.clone()
    };
    let total_dataset_taxa = resolved_dataset_taxa.unwrap_or_else(|| {
        raw_alignments
            .iter()
            .flat_map(|alignment| alignment.taxa.iter())
            .collect::<BTreeSet<_>>()
            .len()
    });

    for raw in &raw_alignments {
        let (transformed, diff) = apply_recipe(raw, &runtime_recipe, total_dataset_taxa);
        if (!config.only_passing || diff.pass) && transformed.length > 0 && !transformed.taxa.is_empty() {
            passing_alignments.push(transformed);
        }
    }

    if passing_alignments.is_empty() {
        return Err("No passing alignments available for concatenation".to_string());
    }

    let (final_taxa, final_seqs, partitions) = build_supermatrix(&passing_alignments);
    let total_length: usize = partitions.iter().map(|part| part.length).sum();

    // Write supermatrix
    let out_prefix = &config.output_file_prefix;
    let supermatrix_path = format!("{}.{}", out_prefix, config.output_format.extension());
    write_alignment(
        &supermatrix_path,
        &final_taxa,
        &final_seqs,
        config.output_format,
    )?;

    // Write RAxML partition file
    let raxml_partition_path = if config.write_raxml_partitions {
        let raxml_path = format!("{}_partitions.txt", out_prefix);
        write_raxml_partitions(&raxml_path, &partitions)?;
        Some(raxml_path)
    } else {
        None
    };

    // Write NEXUS / IQ-TREE partition file
    let nexus_partition_path = if config.write_nexus_partitions {
        let nex_path = format!("{}_partitions.nex", out_prefix);
        write_nexus_partitions(&nex_path, &partitions)?;
        Some(nex_path)
    } else {
        None
    };

    Ok(ConcatenateResult {
        total_taxa: final_taxa.len(),
        total_length,
        total_loci: passing_alignments.len(),
        supermatrix_path,
        raxml_partition_path,
        nexus_partition_path,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn test_concatenate_alignments() {
        let temp_dir = std::env::temp_dir().join("test_concat_dir");
        let _ = fs::create_dir_all(&temp_dir);
        let out_prefix = temp_dir.join("supermatrix");

        let test_input1 = "../test_data/uce-1002.phy";
        let test_input2 = "../test_data/uce-1003.fa";

        if std::path::Path::new(test_input1).exists() && std::path::Path::new(test_input2).exists() {
            let config = ConcatenateConfig {
                input_paths: vec![test_input1.to_string(), test_input2.to_string()],
                output_file_prefix: out_prefix.to_string_lossy().to_string(),
                output_format: AlignmentFormat::Phylip,
                only_passing: false,
                write_raxml_partitions: true,
                write_nexus_partitions: true,
            };
            // This test exercises concatenation, not sample-filter ordering.
            let mut recipe = TrimmingRecipe::default();
            recipe.trim_coverage = false;
            let res = concatenate_alignments(&config, &recipe, None).unwrap();

            assert_eq!(res.total_loci, 2);
            assert_eq!(res.total_taxa, 4);
            assert!(temp_dir.join("supermatrix_partitions.txt").exists());
            assert!(temp_dir.join("supermatrix_partitions.nex").exists());
        }

        let _ = fs::remove_dir_all(temp_dir);
    }
}
