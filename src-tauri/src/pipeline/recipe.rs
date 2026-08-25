use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use crate::algorithms::coverage::RelativeWidth;
use crate::algorithms::orf::{GeneticCode, OrfSearchMode, StopCodonAction};
use crate::algorithms::sanitize::AmbiguityStrategy;
use crate::algorithms::statistical_columns::{
    StatColGapTreatment, StatisticalColumnMethod, TrimalHeuristic,
};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrimmingRecipe {
    pub name: String,
    pub description: String,

    // Step 1: Sanitation
    pub replace_n_with_gap: bool,
    pub ambiguity_strategy: AmbiguityStrategy,
    #[serde(default = "default_true")]
    pub remove_gap_only_columns: bool,

    // Step 2: Outlier / Paralog Similarity Filter
    pub trim_similarity: bool,
    pub similarity_threshold: f64,

    // Step 3a: Profile HMM Segment Cleaner (TAPIR-Style)
    #[serde(default)]
    pub trim_hmm: bool,
    #[serde(default = "default_hmm_posterior")]
    pub hmm_min_posterior: f64,
    #[serde(default = "default_hmm_seg_len")]
    pub hmm_min_segment_length: usize,
    #[serde(default = "default_hmm_island_len")]
    pub hmm_min_island_length: usize,

    // Step 3b: Sliding Window Segment Masking
    pub trim_segments: bool,
    pub segment_window_size: usize,
    pub segment_threshold: f64,

    // Step 4: Candidate Open Reading Frame Extraction & Codon QC (Exons)
    #[serde(default)]
    pub enable_orf: bool,
    #[serde(default = "default_true")]
    pub auto_shift_frame: bool,
    #[serde(default = "default_true")]
    pub auto_flip_reverse: bool,
    #[serde(default)]
    pub stop_codon_action: StopCodonAction,
    #[serde(default = "default_true")]
    pub macse_trim_terminal: bool,
    #[serde(default = "default_macse_sample")]
    pub macse_max_internal_sample: usize,
    #[serde(default = "default_macse_locus")]
    pub macse_max_internal_locus: usize,
    #[serde(default = "default_stop_sample")]
    pub max_stop_codons_sample: usize,
    #[serde(default = "default_stop_locus")]
    pub max_stop_codons_locus: usize,

    #[serde(default)]
    pub genetic_code: GeneticCode,
    #[serde(default)]
    pub orf_search_mode: OrfSearchMode,
    #[serde(default = "default_orf_shared_support")]
    pub orf_min_shared_support_percent: f64,
    #[serde(default = "default_orf_segment_aa")]
    pub orf_min_segment_aa: usize,
    #[serde(default = "default_orf_coding_score")]
    pub orf_min_coding_score: f64,
    #[serde(default = "default_true")]
    pub exclude_uce: bool,
    #[serde(default = "default_true")]
    pub fail_if_no_orf: bool,
    #[serde(default, skip_serializing)]
    pub orf_use_references: bool,
    #[serde(default, skip_serializing)]
    pub orf_reference_sequences: HashMap<String, String>,

    // Step 5: External Edge Trimming
    pub trim_external: bool,
    pub min_external_percent: f64,
    pub codon_preserving: bool,

    // Step 6a: Column Gap Trimming
    pub trim_columns: bool,
    pub min_column_gap_percent: f64,
    pub count_n_as_gap: bool,

    // Step 6b: Statistical Column Trimming (trimAl & Gblocks)
    #[serde(default)]
    pub enable_statistical_columns: bool,
    #[serde(default)]
    pub stat_col_method: StatisticalColumnMethod,
    #[serde(default = "default_similarity_thresh")]
    pub stat_col_similarity_threshold: f64,
    #[serde(default = "default_window_3")]
    pub stat_col_window_size: usize,
    #[serde(default)]
    pub stat_col_heuristic: TrimalHeuristic,
    #[serde(default = "default_min_block_len")]
    pub stat_col_min_block_length: usize,
    #[serde(default = "default_max_nonconserved")]
    pub stat_col_max_nonconserved: usize,
    #[serde(default)]
    pub stat_col_gap_treatment: StatColGapTreatment,
    #[serde(default = "default_entropy_thresh")]
    pub stat_col_entropy_threshold: f64,

    // Step 7: Sample Coverage Trimming
    pub trim_coverage: bool,
    pub min_coverage_bp: usize,
    pub min_coverage_percent: f64,
    pub relative_width: RelativeWidth,
    #[serde(default)]
    pub min_sample_locus_occupancy_percent: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub excluded_taxa: Vec<String>,

    // Step 8: Locus Assessment Gating
    pub assess_alignment: bool,
    pub min_taxa: usize,
    #[serde(default)]
    pub min_taxa_occupancy_percent: f64,
    pub min_length: usize,
    pub max_gap_percent: f64,
    #[serde(default)]
    pub min_pis_count: usize,
    #[serde(default)]
    pub min_pis_percent: f64,
    #[serde(default)]
    pub min_variable_count: usize,
    #[serde(default)]
    pub min_variable_percent: f64,
}

fn default_true() -> bool {
    true
}
fn default_macse_sample() -> usize { 3 }
fn default_macse_locus() -> usize { 10 }
fn default_stop_sample() -> usize { 2 }
fn default_stop_locus() -> usize { 5 }


fn default_orf_shared_support() -> f64 {
    90.0
}

fn default_orf_segment_aa() -> usize {
    35
}

fn default_orf_coding_score() -> f64 {
    40.0
}

fn default_hmm_posterior() -> f64 {
    0.45
}

fn default_hmm_seg_len() -> usize {
    8
}

fn default_hmm_island_len() -> usize {
    20
}

fn default_similarity_thresh() -> f64 {
    0.35
}

fn default_window_3() -> usize {
    3
}

fn default_min_block_len() -> usize {
    5
}

fn default_max_nonconserved() -> usize {
    4
}

fn default_entropy_thresh() -> f64 {
    1.5
}

impl Default for TrimmingRecipe {
    fn default() -> Self {
        Self {
            name: "AlignmentForge Default".to_string(),
            description: "Standard balanced phylogenomic filtering pipeline".to_string(),
            replace_n_with_gap: true,
            ambiguity_strategy: AmbiguityStrategy::Keep,
            remove_gap_only_columns: true,
            trim_similarity: true,
            similarity_threshold: 0.40,
            trim_hmm: false,
            hmm_min_posterior: 0.45,
            hmm_min_segment_length: 8,
            hmm_min_island_length: 20,
            trim_segments: false,
            segment_window_size: 100,
            segment_threshold: 0.45,
            enable_orf: false,
            auto_shift_frame: true,
            auto_flip_reverse: true,
            stop_codon_action: StopCodonAction::RemoveSample,
            macse_trim_terminal: true,
            macse_max_internal_sample: 3,
            macse_max_internal_locus: 10,
            max_stop_codons_sample: 2,
            max_stop_codons_locus: 5,
            genetic_code: GeneticCode::Standard,
            orf_search_mode: OrfSearchMode::ContinuousCds,
            orf_min_shared_support_percent: 90.0,
            orf_min_segment_aa: 35,
            orf_min_coding_score: 40.0,
            exclude_uce: true,
            fail_if_no_orf: false,
            orf_use_references: false,
            orf_reference_sequences: HashMap::new(),
            trim_external: true,
            min_external_percent: 50.0,
            codon_preserving: false,
            trim_columns: false,
            min_column_gap_percent: 60.0,
            count_n_as_gap: true,
            enable_statistical_columns: false,
            stat_col_method: StatisticalColumnMethod::TrimalSimilarity,
            stat_col_similarity_threshold: 0.35,
            stat_col_window_size: 3,
            stat_col_heuristic: TrimalHeuristic::Custom,
            stat_col_min_block_length: 5,
            stat_col_max_nonconserved: 4,
            stat_col_gap_treatment: StatColGapTreatment::Half,
            stat_col_entropy_threshold: 1.5,
            trim_coverage: true,
            min_coverage_bp: 60,
            min_coverage_percent: 50.0,
            relative_width: RelativeWidth::Sample,
            min_sample_locus_occupancy_percent: 0.0,
            excluded_taxa: Vec::new(),
            assess_alignment: true,
            min_taxa: 4,
            min_taxa_occupancy_percent: 50.0,
            min_length: 100,
            max_gap_percent: 50.0,
            min_pis_count: 0,
            min_pis_percent: 0.0,
            min_variable_count: 0,
            min_variable_percent: 0.0,
        }
    }
}

impl TrimmingRecipe {
    pub fn preset_strict() -> Self {
        Self {
            name: "Strict Phylogenomics".to_string(),
            description: "Aggressive filtering for high-confidence core ortholog sets".to_string(),
            replace_n_with_gap: true,
            ambiguity_strategy: AmbiguityStrategy::FixedStandard,
            remove_gap_only_columns: true,
            trim_similarity: true,
            similarity_threshold: 0.35,
            trim_hmm: true,
            hmm_min_posterior: 0.55,
            hmm_min_segment_length: 6,
            hmm_min_island_length: 20,
            trim_segments: true,
            segment_window_size: 80,
            segment_threshold: 0.40,
            enable_orf: false,
            auto_shift_frame: true,
            auto_flip_reverse: true,
            stop_codon_action: StopCodonAction::RemoveSample,
            macse_trim_terminal: true,
            macse_max_internal_sample: 3,
            macse_max_internal_locus: 10,
            max_stop_codons_sample: 2,
            max_stop_codons_locus: 5,
            genetic_code: GeneticCode::Standard,
            orf_search_mode: OrfSearchMode::ContinuousCds,
            orf_min_shared_support_percent: 90.0,
            orf_min_segment_aa: 35,
            orf_min_coding_score: 40.0,
            exclude_uce: true,
            fail_if_no_orf: false,
            orf_use_references: false,
            orf_reference_sequences: HashMap::new(),
            trim_external: true,
            min_external_percent: 70.0,
            codon_preserving: false,
            trim_columns: true,
            min_column_gap_percent: 50.0,
            count_n_as_gap: true,
            enable_statistical_columns: true,
            stat_col_method: StatisticalColumnMethod::TrimalSimilarity,
            stat_col_similarity_threshold: 0.40,
            stat_col_window_size: 3,
            stat_col_heuristic: TrimalHeuristic::Strict,
            stat_col_min_block_length: 6,
            stat_col_max_nonconserved: 3,
            stat_col_gap_treatment: StatColGapTreatment::None,
            stat_col_entropy_threshold: 1.2,
            trim_coverage: true,
            min_coverage_bp: 150,
            min_coverage_percent: 65.0,
            relative_width: RelativeWidth::Alignment,
            min_sample_locus_occupancy_percent: 0.0,
            excluded_taxa: Vec::new(),
            assess_alignment: true,
            min_taxa: 4,
            min_taxa_occupancy_percent: 60.0,
            min_length: 200,
            max_gap_percent: 35.0,
            min_pis_count: 0,
            min_pis_percent: 0.0,
            min_variable_count: 0,
            min_variable_percent: 0.0,
        }
    }

    pub fn preset_relaxed_uce() -> Self {
        Self {
            name: "Relaxed UCE / Sequence Capture".to_string(),
            description: "Permissive filtering preserving maximum informative flanking sequence".to_string(),
            replace_n_with_gap: true,
            ambiguity_strategy: AmbiguityStrategy::Keep,
            remove_gap_only_columns: true,
            trim_similarity: true,
            similarity_threshold: 0.45,
            trim_hmm: true,
            hmm_min_posterior: 0.35,
            hmm_min_segment_length: 12,
            hmm_min_island_length: 30,
            trim_segments: false,
            segment_window_size: 100,
            segment_threshold: 0.45,
            enable_orf: false,
            auto_shift_frame: true,
            auto_flip_reverse: true,
            stop_codon_action: StopCodonAction::RemoveSample,
            macse_trim_terminal: true,
            macse_max_internal_sample: 3,
            macse_max_internal_locus: 10,
            max_stop_codons_sample: 2,
            max_stop_codons_locus: 5,
            genetic_code: GeneticCode::Standard,
            orf_search_mode: OrfSearchMode::ContinuousCds,
            orf_min_shared_support_percent: 90.0,
            orf_min_segment_aa: 35,
            orf_min_coding_score: 40.0,
            exclude_uce: true,
            fail_if_no_orf: false,
            orf_use_references: false,
            orf_reference_sequences: HashMap::new(),
            trim_external: true,
            min_external_percent: 40.0,
            codon_preserving: false,
            trim_columns: true,
            min_column_gap_percent: 80.0,
            count_n_as_gap: true,
            enable_statistical_columns: false,
            stat_col_method: StatisticalColumnMethod::TrimalSimilarity,
            stat_col_similarity_threshold: 0.30,
            stat_col_window_size: 3,
            stat_col_heuristic: TrimalHeuristic::Gappyout,
            stat_col_min_block_length: 4,
            stat_col_max_nonconserved: 5,
            stat_col_gap_treatment: StatColGapTreatment::All,
            stat_col_entropy_threshold: 1.7,
            trim_coverage: true,
            min_coverage_bp: 40,
            min_coverage_percent: 40.0,
            relative_width: RelativeWidth::Sample,
            min_sample_locus_occupancy_percent: 0.0,
            excluded_taxa: Vec::new(),
            assess_alignment: true,
            min_taxa: 4,
            min_taxa_occupancy_percent: 35.0,
            min_length: 80,
            max_gap_percent: 60.0,
            min_pis_count: 0,
            min_pis_percent: 0.0,
            min_variable_count: 0,
            min_variable_percent: 0.0,
        }
    }

    pub fn preset_exon_codon() -> Self {
        Self {
            name: "Exon / Coding Region".to_string(),
            description: "Frame-ready exon filtering; enable ORF analysis explicitly in its workspace".to_string(),
            replace_n_with_gap: true,
            ambiguity_strategy: AmbiguityStrategy::MajorityBase,
            remove_gap_only_columns: true,
            trim_similarity: true,
            similarity_threshold: 0.35,
            trim_hmm: true,
            hmm_min_posterior: 0.45,
            hmm_min_segment_length: 9,
            hmm_min_island_length: 20,
            trim_segments: false,
            segment_window_size: 100,
            segment_threshold: 0.45,
            enable_orf: false,
            auto_shift_frame: true,
            auto_flip_reverse: true,
            stop_codon_action: StopCodonAction::RemoveSample,
            macse_trim_terminal: true,
            macse_max_internal_sample: 3,
            macse_max_internal_locus: 10,
            max_stop_codons_sample: 2,
            max_stop_codons_locus: 5,
            genetic_code: GeneticCode::Standard,
            orf_search_mode: OrfSearchMode::ContinuousCds,
            orf_min_shared_support_percent: 90.0,
            orf_min_segment_aa: 35,
            orf_min_coding_score: 40.0,
            exclude_uce: true,
            fail_if_no_orf: false,
            orf_use_references: false,
            orf_reference_sequences: HashMap::new(),
            trim_external: true,
            min_external_percent: 60.0,
            codon_preserving: true,
            trim_columns: true,
            min_column_gap_percent: 50.0,
            count_n_as_gap: true,
            enable_statistical_columns: false,
            stat_col_method: StatisticalColumnMethod::TrimalSimilarity,
            stat_col_similarity_threshold: 0.35,
            stat_col_window_size: 3,
            stat_col_heuristic: TrimalHeuristic::Custom,
            stat_col_min_block_length: 5,
            stat_col_max_nonconserved: 4,
            stat_col_gap_treatment: StatColGapTreatment::Half,
            stat_col_entropy_threshold: 1.5,
            trim_coverage: true,
            min_coverage_bp: 90,
            min_coverage_percent: 50.0,
            relative_width: RelativeWidth::Sample,
            min_sample_locus_occupancy_percent: 0.0,
            excluded_taxa: Vec::new(),
            assess_alignment: true,
            min_taxa: 4,
            min_taxa_occupancy_percent: 50.0,
            min_length: 120,
            max_gap_percent: 40.0,
            min_pis_count: 0,
            min_pis_percent: 0.0,
            min_variable_count: 0,
            min_variable_percent: 0.0,
        }
    }
}
