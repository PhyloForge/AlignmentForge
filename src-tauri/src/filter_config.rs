use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

use crate::pipeline::recipe::TrimmingRecipe;

const FILTER_CONFIG_VERSION: u32 = 1;

#[derive(Debug, Clone, Serialize, Deserialize)]
struct FilterConfigDocument {
    #[serde(default = "default_config_version")]
    config_version: u32,
    #[serde(flatten)]
    recipe: TrimmingRecipe,
}

fn default_config_version() -> u32 {
    FILTER_CONFIG_VERSION
}

fn section_comment(key: &str) -> Option<&'static str> {
    match key {
        "name" => Some("# Preset metadata"),
        "replace_n_with_gap" => Some(
            "# Sanitation\n# ambiguity_strategy: keep | majoritybase | converttogap | fixedstandard\n# remove_gap_only_columns removes columns containing only -, N, or ?",
        ),
        "trim_similarity" => Some("# Paralog / divergence pruning"),
        "trim_hmm" => Some("# Profile HMM cleaner"),
        "trim_segments" => Some("# Sliding-window segment masking"),
        "enable_orf" => Some(
            "# Candidate open reading frame extraction / codon optimization\n# orf_search_mode: continuouscds | bestsharedsegment | referenceguided | referencecandidateorf\n# stop_codon_action: removesample | maskcodon | keep\n# genetic_code: standard | vertebratemitochondrial | invertebratemitochondrial\n# bestsharedsegment uses shared stop-free length plus an internal protein-profile coding score\n# referencecandidateorf tries a matched reference first, then candidate extraction\n# UCE and explicitly non-coding IDs remain skipped",
        ),
        "trim_external" => Some("# Ragged-edge trimming"),
        "trim_columns" => Some("# Column gap filter"),
        "enable_statistical_columns" => Some(
            "# Statistical column trimming\n# stat_col_method: none | trimalsimilarity | gblocksblocks | entropy\n# stat_col_heuristic: custom | gappyout | strict | strictplus\n# stat_col_gap_treatment: none | half | all",
        ),
        "trim_coverage" => Some(
            "# Sample coverage filter\n# relative_width: sample | alignment",
        ),
        "assess_alignment" => Some("# Locus pass / fail criteria"),
        _ => None,
    }
}

fn validate_fraction(name: &str, value: f64) -> Result<(), String> {
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        Ok(())
    } else {
        Err(format!("{name} must be between 0 and 1"))
    }
}

fn validate_percent(name: &str, value: f64) -> Result<(), String> {
    if value.is_finite() && (0.0..=100.0).contains(&value) {
        Ok(())
    } else {
        Err(format!("{name} must be between 0 and 100"))
    }
}

fn validate_recipe(recipe: &TrimmingRecipe) -> Result<(), String> {
    if recipe.name.trim().is_empty() {
        return Err("name cannot be empty".to_string());
    }

    validate_fraction("similarity_threshold", recipe.similarity_threshold)?;
    validate_fraction("hmm_min_posterior", recipe.hmm_min_posterior)?;
    validate_fraction("segment_threshold", recipe.segment_threshold)?;
    validate_fraction(
        "stat_col_similarity_threshold",
        recipe.stat_col_similarity_threshold,
    )?;

    validate_percent("min_external_percent", recipe.min_external_percent)?;
    validate_percent("min_column_gap_percent", recipe.min_column_gap_percent)?;
    validate_percent("min_coverage_percent", recipe.min_coverage_percent)?;
    validate_percent(
        "orf_min_shared_support_percent",
        recipe.orf_min_shared_support_percent,
    )?;
    validate_percent("orf_min_coding_score", recipe.orf_min_coding_score)?;
    validate_percent(
        "min_sample_locus_occupancy_percent",
        recipe.min_sample_locus_occupancy_percent,
    )?;
    validate_percent(
        "min_taxa_occupancy_percent",
        recipe.min_taxa_occupancy_percent,
    )?;
    validate_percent("max_gap_percent", recipe.max_gap_percent)?;
    validate_percent("min_pis_percent", recipe.min_pis_percent)?;
    validate_percent("min_variable_percent", recipe.min_variable_percent)?;

    if recipe.hmm_min_segment_length == 0 {
        return Err("hmm_min_segment_length must be at least 1".to_string());
    }
    if recipe.orf_min_segment_aa == 0 {
        return Err("orf_min_segment_aa must be at least 1".to_string());
    }
    if recipe.segment_window_size == 0 {
        return Err("segment_window_size must be at least 1".to_string());
    }
    if recipe.stat_col_window_size == 0 {
        return Err("stat_col_window_size must be at least 1".to_string());
    }
    if recipe.min_taxa == 0 {
        return Err("min_taxa must be at least 1".to_string());
    }
    if !recipe.stat_col_entropy_threshold.is_finite() || recipe.stat_col_entropy_threshold < 0.0 {
        return Err("stat_col_entropy_threshold must be zero or greater".to_string());
    }

    Ok(())
}

pub fn serialize_filter_config(recipe: &TrimmingRecipe) -> Result<String, String> {
    validate_recipe(recipe)?;
    let document = FilterConfigDocument {
        config_version: FILTER_CONFIG_VERSION,
        recipe: recipe.clone(),
    };
    let body = toml::to_string_pretty(&document)
        .map_err(|error| format!("Could not serialize filter config: {error}"))?;

    let mut output = String::from(
        "# AlignmentForge filter configuration\n\
# Edit values after '=' and keep string values inside quotes.\n\
# Fractions use 0 to 1; percentage fields use 0 to 100.\n",
    );

    for line in body.lines() {
        let key = line
            .split_once('=')
            .map(|(key, _)| key.trim())
            .unwrap_or("");
        if let Some(comment) = section_comment(key) {
            output.push('\n');
            output.push_str(comment);
            output.push('\n');
        }
        output.push_str(line);
        output.push('\n');
    }

    Ok(output)
}

pub fn parse_filter_config(contents: &str) -> Result<TrimmingRecipe, String> {
    let document: FilterConfigDocument = toml::from_str(contents)
        .map_err(|error| format!("Invalid AlignmentForge filter config: {error}"))?;

    if document.config_version != FILTER_CONFIG_VERSION {
        return Err(format!(
            "Unsupported filter config version {} (this app supports version {})",
            document.config_version, FILTER_CONFIG_VERSION
        ));
    }

    validate_recipe(&document.recipe)
        .map_err(|error| format!("Invalid AlignmentForge filter config: {error}"))?;
    Ok(document.recipe)
}

pub fn save_filter_config(path: &Path, recipe: &TrimmingRecipe) -> Result<String, String> {
    let contents = serialize_filter_config(recipe)?;
    fs::write(path, contents).map_err(|error| {
        format!(
            "Could not write filter config '{}': {error}",
            path.display()
        )
    })?;
    Ok(path.to_string_lossy().to_string())
}

pub fn load_filter_config(path: &Path) -> Result<TrimmingRecipe, String> {
    let contents = fs::read_to_string(path)
        .map_err(|error| format!("Could not read filter config '{}': {error}", path.display()))?;
    parse_filter_config(&contents)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn filter_config_round_trip_preserves_recipe() {
        let recipe = TrimmingRecipe::default();
        let contents = serialize_filter_config(&recipe).unwrap();
        let parsed = parse_filter_config(&contents).unwrap();

        assert!(contents.contains("# Sample coverage filter"));
        assert!(contents.contains("config_version = 1"));
        assert_eq!(parsed.name, recipe.name);
        assert_eq!(parsed.trim_hmm, recipe.trim_hmm);
        assert_eq!(parsed.min_coverage_bp, recipe.min_coverage_bp);
        assert_eq!(parsed.min_taxa_occupancy_percent, 50.0);
    }

    #[test]
    fn filter_config_rejects_out_of_range_values() {
        let recipe = TrimmingRecipe::default();
        let contents = serialize_filter_config(&recipe)
            .unwrap()
            .replace("similarity_threshold = 0.4", "similarity_threshold = 1.4");
        let error = parse_filter_config(&contents).unwrap_err();

        assert!(error.contains("similarity_threshold must be between 0 and 1"));
    }
}
