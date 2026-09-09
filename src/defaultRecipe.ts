import { TrimmingRecipe } from './types';

/**
 * The default trimming recipe.
 *
 * This is the single definition used by the browser build. It must stay in step
 * with `TrimmingRecipe::default()` in `src-tauri/src/pipeline/recipe.rs`; the
 * `npm run parity` check fails when the two drift apart. The desktop build
 * replaces this with the presets returned by the `get_presets` command.
 */
export const DEFAULT_RECIPE: TrimmingRecipe = {
  name: 'AlignmentForge Default',
  description: 'Standard balanced phylogenomic filtering pipeline',

  replace_n_with_gap: true,
  ambiguity_strategy: 'keep',
  remove_gap_only_columns: true,

  trim_similarity: true,
  similarity_threshold: 0.4,

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
  stop_codon_action: 'removesample',
  macse_trim_terminal: true,
  macse_max_internal_sample: 3,
  macse_max_internal_locus: 10,
  max_stop_codons_sample: 2,
  max_stop_codons_locus: 5,
  genetic_code: 'standard',
  orf_search_mode: 'continuouscds',
  orf_min_shared_support_percent: 90.0,
  orf_min_segment_aa: 35,
  orf_min_coding_score: 40.0,
  exclude_uce: true,
  fail_if_no_orf: false,
  orf_use_references: false,
  orf_reference_sequences: {},

  trim_external: true,
  min_external_percent: 50.0,
  codon_preserving: false,

  trim_columns: false,
  min_column_gap_percent: 60.0,
  count_n_as_gap: true,

  enable_statistical_columns: false,
  stat_col_method: 'trimalsimilarity',
  stat_col_similarity_threshold: 0.35,
  stat_col_window_size: 3,
  stat_col_heuristic: 'custom',
  stat_col_min_block_length: 5,
  stat_col_max_nonconserved: 4,
  stat_col_gap_treatment: 'half',
  stat_col_entropy_threshold: 1.5,

  trim_coverage: true,
  min_coverage_bp: 60,
  min_coverage_percent: 50.0,
  relative_width: 'sample',
  min_sample_locus_occupancy_percent: 0.0,
  discarded_taxa: [],
  excluded_taxa: [],

  assess_alignment: true,
  min_taxa: 4,
  min_taxa_occupancy_percent: 50.0,
  min_length: 100,
  max_gap_percent: 50.0,
  min_pis_count: 0,
  min_pis_percent: 0.0,
  min_variable_count: 0,
  min_variable_percent: 0.0,
};
