use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use serde::{Deserialize, Serialize};

use crate::models::{Alignment, AlignmentFormat};
use crate::parsers::{parse_alignment, write_alignment};
use crate::pipeline::catalog::recipe_with_dataset_sample_filter;
use crate::pipeline::engine::apply_recipe;
use crate::pipeline::recipe::TrimmingRecipe;
use crate::export::concatenate::{build_supermatrix, write_nexus_partitions, write_raxml_partitions};

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

/// Chooses the gene map delimiter. A tab in the first line marks a
/// tab-separated file, whatever its extension.
fn gene_map_delimiter(path: &str, content: &str) -> u8 {
    let is_tsv = Path::new(path)
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("tsv"));
    let first_line = content.lines().find(|line| !line.trim().is_empty()).unwrap_or_default();
    if is_tsv || first_line.contains('\t') {
        b'\t'
    } else {
        b','
    }
}

pub fn concatenate_alignments_by_gene(
    config: &GroupedConcatenateConfig,
    recipe: &TrimmingRecipe,
    resolved_dataset_taxa: Option<usize>,
) -> Result<GroupedConcatenateResult, String> {
    // 1. Parse the metadata file (CSV, TSV, or TXT)
    let content = std::fs::read_to_string(&config.gene_mapping_csv_path)
        .map_err(|e| format!("Could not open gene mapping file: {e}"))?;
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(gene_map_delimiter(&config.gene_mapping_csv_path, &content))
        .from_reader(content.as_bytes());

    // exon_name -> gene_name
    let mut exon_to_gene: BTreeMap<String, String> = BTreeMap::new();
    
    for record in reader.records() {
        let record = record.map_err(|error| format!("Could not read gene mapping: {error}"))?;
        if record.len() < 2 {
            continue;
        }
        let exon = record[0].trim().to_string();
        let gene = record[1].trim().to_string();
        if exon.is_empty() || gene.is_empty() {
            continue;
        }
        let gene_path = Path::new(&gene);
        if gene == "."
            || gene == ".."
            || gene_path.components().count() != 1
            || gene.contains('/')
            || gene.contains('\\')
        {
            return Err(format!(
                "Gene name '{gene}' must be a single filename component"
            ));
        }
        exon_to_gene.insert(exon, gene);
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

    let all_raw_alignments: Vec<Alignment> = config
        .input_paths
        .iter()
        .map(|path| parse_alignment(Path::new(path)))
        .collect::<Result<_, _>>()?;
    let runtime_recipe = if resolved_dataset_taxa.is_none() {
        recipe_with_dataset_sample_filter(recipe, &all_raw_alignments)
    } else {
        recipe.clone()
    };
    let total_dataset_taxa = resolved_dataset_taxa.unwrap_or_else(|| {
        all_raw_alignments
            .iter()
            .flat_map(|alignment| alignment.taxa.iter())
            .collect::<BTreeSet<_>>()
            .len()
    });

    // 3. Process each gene group
    for (gene, paths) in &gene_to_paths {
        let raw_alignments: Vec<&Alignment> = paths
            .iter()
            .filter_map(|path| all_raw_alignments.iter().find(|alignment| &alignment.file_path == path))
            .collect();
            
        if raw_alignments.is_empty() {
            continue;
        }

        let mut passing_alignments = Vec::new();

        for raw in &raw_alignments {
            let (transformed, diff) = apply_recipe(raw, &runtime_recipe, total_dataset_taxa);
            if (!config.only_passing || diff.pass) && transformed.length > 0 && !transformed.taxa.is_empty() {
                passing_alignments.push(transformed);
            }
        }

        if passing_alignments.is_empty() {
            continue;
        }

        let (final_taxa, final_seqs, partitions) = build_supermatrix(&passing_alignments);
        total_exons_processed += passing_alignments.len();
        let out_prefix = Path::new(&config.output_directory).join(gene);
        let out_prefix_str = out_prefix.to_string_lossy();
        let supermatrix_path = format!("{}.{}", out_prefix_str, config.output_format.extension());

        write_alignment(
            &supermatrix_path,
            &final_taxa,
            &final_seqs,
            config.output_format,
        )
        .map_err(|error| format!("Failed to write gene alignment '{gene}': {error}"))?;

        if config.write_raxml_partitions {
            write_raxml_partitions(&format!("{}_partitions.txt", out_prefix_str), &partitions)?;
        }

        if config.write_nexus_partitions {
            write_nexus_partitions(&format!("{}_partitions.nex", out_prefix_str), &partitions)?;
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
    fn tab_separated_gene_map_is_detected_from_its_content() {
        assert_eq!(gene_map_delimiter("map.txt", "exon1\tgeneA\n"), b'\t');
        assert_eq!(gene_map_delimiter("map.tsv", "exon1,geneA\n"), b'\t');
        assert_eq!(gene_map_delimiter("map.txt", "\nexon1,geneA\n"), b',');
    }

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

            let res = concatenate_alignments_by_gene(&config, &recipe, None).unwrap();
            
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
