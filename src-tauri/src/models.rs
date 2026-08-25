use serde::{Deserialize, Serialize};

fn default_true() -> bool {
    true
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum AlignmentFormat {
    Fasta,
    Phylip,
    Nexus,
}

impl AlignmentFormat {
    pub fn extension(&self) -> &'static str {
        match self {
            AlignmentFormat::Fasta => "fa",
            AlignmentFormat::Phylip => "phy",
            AlignmentFormat::Nexus => "nex",
        }
    }
}

/// A parsed multiple sequence alignment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Alignment {
    pub id: String,
    pub file_name: String,
    pub file_path: String,
    pub format: AlignmentFormat,
    pub taxa: Vec<String>,
    pub sequences: Vec<String>,
    pub length: usize,
    pub num_taxa: usize,
}

impl Alignment {
    pub fn new(
        id: String,
        file_name: String,
        file_path: String,
        format: AlignmentFormat,
        taxa: Vec<String>,
        sequences: Vec<String>,
    ) -> Self {
        let length = sequences.first().map_or(0, |s| s.len());
        let num_taxa = taxa.len();
        Self {
            id,
            file_name,
            file_path,
            format,
            taxa,
            sequences,
            length,
            num_taxa,
        }
    }
}

/// Fast summary metrics for an alignment used in Catalog and QC views (calculated post-trimming)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentSummary {
    pub id: String,
    pub file_name: String,
    pub file_path: String,
    pub format: AlignmentFormat,
    pub num_taxa: usize,
    pub length: usize,
    pub total_basepairs: usize,
    pub gap_count: usize,
    pub gap_percent: f64,
    #[serde(default)]
    pub variable_count: usize,
    #[serde(default)]
    pub variable_percent: f64,
    pub pis_count: usize,
    pub pis_percent: f64,
    pub mean_divergence: f64,
    pub gc_percent: f64,
    pub pass: bool,
    pub fail_reasons: Vec<String>,
    #[serde(default = "default_true")]
    pub orf_valid: bool,
    #[serde(default)]
    pub orf_evaluated: bool,
    #[serde(default)]
    pub orf_candidate_found: bool,
    #[serde(default)]
    pub orf_frame: Option<i8>,
    #[serde(default)]
    pub orf_start: Option<usize>,
    #[serde(default)]
    pub orf_end: Option<usize>,
    #[serde(default)]
    pub orf_support_count: usize,
    #[serde(default)]
    pub orf_support_percent: f64,
    #[serde(default)]
    pub orf_retained_samples: usize,
    #[serde(default)]
    pub orf_candidate_length_aa: usize,
    #[serde(default)]
    pub orf_coding_score: f64,
    #[serde(default)]
    pub orf_amino_acid_conservation: f64,
    #[serde(default)]
    pub orf_frame_contrast: f64,
    #[serde(default)]
    pub orf_reference_evaluated: bool,
    #[serde(default)]
    pub orf_reference_matched: bool,
    #[serde(default)]
    pub orf_reference_identity: f64,
    #[serde(default)]
    pub orf_reference_coverage: f64,
    #[serde(default)]
    pub orf_intron_length: usize,
    #[serde(default)]
    pub raw_num_taxa: usize,
    #[serde(default)]
    pub raw_length: usize,
    #[serde(default)]
    pub raw_gap_percent: f64,
    /// Compact per-sample retention metadata used by the QC occupancy chart.
    #[serde(default)]
    pub retained_taxa: Vec<String>,
    #[serde(default)]
    pub retained_taxon_basepairs: std::collections::HashMap<String, usize>,
}

/// Overall dataset statistics across all alignments
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatasetOverview {
    pub total_alignments: usize,
    pub passed_alignments: usize,
    pub discarded_alignments: usize,
    pub total_unique_taxa: usize,
    pub mean_taxa: f64,
    pub mean_length: f64,
    pub mean_gap_percent: f64,
    pub mean_pis: f64,
    pub total_matrix_basepairs: usize,
}

/// Visual diff describing what a trimming recipe changes on an alignment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrimmingDiff {
    pub id: String,
    pub old_taxa_count: usize,
    pub new_taxa_count: usize,
    pub dropped_taxa: Vec<String>,
    pub old_length: usize,
    pub new_length: usize,
    pub trimmed_columns: Vec<usize>,
    pub masked_segments: Vec<MaskedSegment>,
    #[serde(default)]
    pub column_reasons: std::collections::HashMap<usize, String>,
    #[serde(default)]
    pub dropped_taxa_reasons: std::collections::HashMap<String, String>,
    #[serde(default)]
    pub stop_codons: Vec<StopCodonPos>,
    #[serde(default)]
    pub final_stop_codons: Vec<StopCodonPos>,
    pub old_gap_percent: f64,
    pub new_gap_percent: f64,
    #[serde(default)]
    pub old_variable: usize,
    #[serde(default)]
    pub new_variable: usize,
    pub old_pis: usize,
    pub new_pis: usize,
    #[serde(default = "default_true")]
    pub found_valid_orf: bool,
    #[serde(default)]
    pub orf_evaluated: bool,
    #[serde(default)]
    pub orf_candidate_found: bool,
    #[serde(default)]
    pub orf_frame: Option<i8>,
    #[serde(default)]
    pub orf_start: Option<usize>,
    #[serde(default)]
    pub orf_end: Option<usize>,
    #[serde(default)]
    pub orf_support_count: usize,
    #[serde(default)]
    pub orf_support_percent: f64,
    #[serde(default)]
    pub orf_retained_samples: usize,
    #[serde(default)]
    pub orf_candidate_length_aa: usize,
    #[serde(default)]
    pub orf_coding_score: f64,
    #[serde(default)]
    pub orf_amino_acid_conservation: f64,
    #[serde(default)]
    pub orf_frame_contrast: f64,
    #[serde(default)]
    pub orf_reference_evaluated: bool,
    #[serde(default)]
    pub orf_reference_matched: bool,
    #[serde(default)]
    pub orf_reference_identity: f64,
    #[serde(default)]
    pub orf_reference_coverage: f64,
    #[serde(default)]
    pub orf_intron_length: usize,
    pub pass: bool,
    pub fail_reasons: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StopCodonPos {
    pub taxon: String,
    pub start: usize,
    pub end: usize,
    pub codon: String,
    pub is_terminal: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaskedSegment {
    pub taxon: String,
    pub start: usize,
    pub end: usize,
}

/// Taxon presence and occupancy metrics across dataset
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TaxonOccupancy {
    pub taxon_name: String,
    pub present_loci_count: usize,
    pub present_loci_percent: f64,
    pub mean_gap_percent: f64,
    pub total_bp: usize,
}
