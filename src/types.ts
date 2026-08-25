export type AlignmentFormat = 'fasta' | 'phylip' | 'nexus';

export type AmbiguityStrategy = 'keep' | 'majoritybase' | 'converttogap' | 'fixedstandard';

export type RelativeWidth = 'alignment' | 'sample';

export type GeneticCode = 'standard' | 'vertebratemitochondrial' | 'invertebratemitochondrial';

export type StopCodonAction = 'removesample' | 'maskcodon' | 'keep';

export type OrfSearchMode =
  | 'continuouscds'
  | 'bestsharedsegment'
  | 'referenceguided'
  | 'referencecandidateorf';

export type StatisticalColumnMethod = 'none' | 'trimalsimilarity' | 'gblocksblocks' | 'entropy';

export type StatColGapTreatment = 'none' | 'half' | 'all';

export type TrimalHeuristic = 'custom' | 'gappyout' | 'strict' | 'strictplus';

export type ViewMode = 'catalog' | 'msa' | 'orf' | 'matrix' | 'qc';

export type ColorScheme =
  | 'clustalx'
  | 'nucleotide'
  | 'nucleotide-colorblind'
  | 'nucleotide-high-contrast'
  | 'nucleotide-pastel'
  | 'nucleotide-ocean'
  | 'nucleotide-sunset'
  | 'nucleotide-forest'
  | 'zappo'
  | 'difference'
  | 'monochrome';

export type AminoAcidColorScheme =
  | 'chemistry'
  | 'colorblind'
  | 'zappo'
  | 'taylor'
  | 'hydrophobicity'
  | 'monochrome';

export interface AminoAcidViewerSettings {
  enabled: boolean;
  hideTerminalStops: boolean;
  colorScheme: AminoAcidColorScheme;
  dimConsensusMatches: boolean;
}

export type CatalogSortField =
  | 'id'
  | 'num_taxa'
  | 'length'
  | 'variable_count'
  | 'variable_percent'
  | 'pis_count'
  | 'pis_percent'
  | 'gap_percent'
  | 'mean_divergence'
  | 'gc_percent'
  | 'orf_candidate_found'
  | 'orf_frame'
  | 'orf_support_percent'
  | 'orf_retained_samples'
  | 'orf_candidate_length_aa'
  | 'orf_coding_score'
  | 'pass';

export interface Alignment {
  id: string;
  file_name: string;
  file_path: string;
  format: AlignmentFormat;
  taxa: string[];
  sequences: string[];
  length: usize;
  num_taxa: usize;
}

export type usize = number;

export interface AlignmentSummary {
  id: string;
  file_name: string;
  file_path: string;
  format: AlignmentFormat;
  num_taxa: number;
  length: number;
  total_basepairs: number;
  gap_count: number;
  gap_percent: number;
  variable_count: number;
  variable_percent: number;
  pis_count: number;
  pis_percent: number;
  mean_divergence: number;
  gc_percent: number;
  pass: boolean;
  fail_reasons: string[];
  orf_valid?: boolean;
  orf_evaluated?: boolean;
  orf_candidate_found?: boolean;
  orf_frame?: number;
  orf_start?: number;
  orf_end?: number;
  orf_support_count?: number;
  orf_support_percent?: number;
  orf_retained_samples?: number;
  orf_candidate_length_aa?: number;
  orf_coding_score?: number;
  orf_amino_acid_conservation?: number;
  orf_frame_contrast?: number;
  orf_reference_evaluated?: boolean;
  orf_reference_matched?: boolean;
  orf_reference_identity?: number;
  orf_reference_coverage?: number;
  orf_intron_length?: number;
  raw_num_taxa?: number;
  raw_length?: number;
  raw_gap_percent?: number;
  /** Samples and usable base-pair totals surviving the current processing recipe. */
  retained_taxa?: string[];
  retained_taxon_basepairs?: Record<string, number>;
}

export interface DatasetOverview {
  total_alignments: number;
  passed_alignments: number;
  discarded_alignments: number;
  total_unique_taxa: number;
  mean_taxa: number;
  mean_length: number;
  mean_gap_percent: number;
  mean_pis: number;
  total_matrix_basepairs: number;
}

export interface TaxonOccupancy {
  taxon_name: string;
  present_loci_count: number;
  present_loci_percent: number;
  mean_gap_percent: number;
  total_bp: number;
}

export interface MaskedSegment {
  taxon: string;
  start: number;
  end: number;
}

export interface StopCodonPos {
  taxon: string;
  start: number;
  end: number;
  codon: string;
  is_terminal: boolean;
}

export interface TrimmingDiff {
  id: string;
  old_taxa_count: number;
  new_taxa_count: number;
  dropped_taxa: string[];
  old_length: number;
  new_length: number;
  trimmed_columns: number[];
  masked_segments: MaskedSegment[];
  column_reasons?: Record<number, string>;
  dropped_taxa_reasons?: Record<string, string>;
  stop_codons?: StopCodonPos[];
  final_stop_codons?: StopCodonPos[];
  old_gap_percent: number;
  new_gap_percent: number;
  old_variable: number;
  new_variable: number;
  old_pis: number;
  new_pis: number;
  found_valid_orf?: boolean;
  orf_evaluated?: boolean;
  orf_candidate_found?: boolean;
  orf_frame?: number;
  orf_start?: number;
  orf_end?: number;
  orf_support_count?: number;
  orf_support_percent?: number;
  orf_retained_samples?: number;
  orf_candidate_length_aa?: number;
  orf_coding_score?: number;
  orf_amino_acid_conservation?: number;
  orf_frame_contrast?: number;
  orf_reference_evaluated?: boolean;
  orf_reference_matched?: boolean;
  orf_reference_identity?: number;
  orf_reference_coverage?: number;
  orf_intron_length?: number;
  pass: boolean;
  fail_reasons: string[];
}

export interface TrimmingRecipe {
  name: string;
  description: string;

  replace_n_with_gap: boolean;
  ambiguity_strategy: AmbiguityStrategy;
  remove_gap_only_columns: boolean;

  trim_similarity: boolean;
  similarity_threshold: number;

  trim_hmm: boolean;
  hmm_min_posterior: number;
  hmm_min_segment_length: number;
  hmm_min_island_length: number;

  trim_segments: boolean;
  segment_window_size: number;
  segment_threshold: number;

  enable_orf: boolean;
  auto_shift_frame: boolean;
  auto_flip_reverse: boolean;
  stop_codon_action: StopCodonAction;
  macse_trim_terminal: boolean;
  macse_max_internal_sample: number;
  macse_max_internal_locus: number;
  max_stop_codons_sample: number;
  max_stop_codons_locus: number;
  genetic_code: GeneticCode;
  orf_search_mode: OrfSearchMode;
  orf_min_shared_support_percent: number;
  orf_min_segment_aa: number;
  orf_min_coding_score: number;
  exclude_uce: boolean;
  fail_if_no_orf: boolean;
  orf_use_references: boolean;
  orf_reference_sequences: Record<string, string>;

  trim_external: boolean;
  min_external_percent: number;
  codon_preserving: boolean;

  trim_columns: boolean;
  min_column_gap_percent: number;
  count_n_as_gap: boolean;

  enable_statistical_columns: boolean;
  stat_col_method: StatisticalColumnMethod;
  stat_col_similarity_threshold: number;
  stat_col_window_size: number;
  stat_col_heuristic: TrimalHeuristic;
  stat_col_min_block_length: number;
  stat_col_max_nonconserved: number;
  stat_col_gap_treatment: StatColGapTreatment;
  stat_col_entropy_threshold: number;

  trim_coverage: boolean;
  min_coverage_bp: number;
  min_coverage_percent: number;
  relative_width: RelativeWidth;
  min_sample_locus_occupancy_percent: number;
  excluded_taxa?: string[];

  assess_alignment: boolean;
  min_taxa: number;
  min_taxa_occupancy_percent: number;
  min_length: number;
  max_gap_percent: number;
  min_pis_count: number;
  min_pis_percent: number;
  min_variable_count: number;
  min_variable_percent: number;
}

export interface ScanResponse {
  summaries: AlignmentSummary[];
  overview: DatasetOverview;
  occupancy: TaxonOccupancy[];
}

export interface CatalogUpdateResponse {
  summaries: AlignmentSummary[];
  overview: DatasetOverview;
}

export interface AlignmentViewResponse {
  raw_alignment: Alignment;
  trimmed_alignment: Alignment;
  diff: TrimmingDiff;
  pis_mask: boolean[];
  majority_consensus: string;
}

export interface BatchExportConfig {
  input_paths: string[];
  output_directory: string;
  output_format: AlignmentFormat;
  only_passing: boolean;
  save_recipe_json: boolean;
  save_summary_csv: boolean;
  export_introns: boolean;
}

export interface BatchExportResult {
  total_processed: number;
  total_exported: number;
  total_discarded: number;
  summary_csv_path?: string;
  recipe_json_path?: string;
  total_introns_exported?: number;
  intron_directory_path?: string;
}

export interface ConcatenateConfig {
  input_paths: string[];
  output_file_prefix: string;
  output_format: AlignmentFormat;
  only_passing: boolean;
  write_raxml_partitions: boolean;
  write_nexus_partitions: boolean;
}

export interface ConcatenateResult {
  total_taxa: number;
  total_length: number;
  total_loci: number;
  supermatrix_path: string;
  raxml_partition_path?: string;
  nexus_partition_path?: string;
}

export interface GroupedConcatenateConfig {
  input_paths: string[];
  output_directory: string;
  gene_mapping_csv_path: string;
  output_format: AlignmentFormat;
  only_passing: boolean;
  write_raxml_partitions: boolean;
  write_nexus_partitions: boolean;
}

export interface GroupedConcatenateResult {
  total_genes: number;
  total_exons_processed: number;
  output_directory: string;
}
