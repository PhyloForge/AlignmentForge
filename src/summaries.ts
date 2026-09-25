import { AlignmentSummary, DatasetOverview, TrimmingRecipe } from './types';

/**
 * Aggregates summaries into the dataset overview.
 *
 * This only sums fields the engine already produced; it makes no filtering or
 * pass and fail decision of its own.
 */
export function buildDatasetOverviewFromSummaries(
  summaries: AlignmentSummary[],
  totalUniqueTaxa: number
): DatasetOverview {
  const totalLoci = summaries.length;
  let passed = 0;
  let sumTaxa = 0;
  let sumLength = 0;
  let sumGap = 0;
  let sumPis = 0;
  let totalBasepairs = 0;

  for (const summary of summaries) {
    if (summary.pass) passed++;
    sumTaxa += summary.num_taxa;
    sumLength += summary.length;
    sumGap += summary.gap_percent;
    sumPis += summary.pis_count;
    totalBasepairs += summary.total_basepairs;
  }

  return {
    total_alignments: totalLoci,
    passed_alignments: passed,
    discarded_alignments: totalLoci - passed,
    total_unique_taxa: totalUniqueTaxa,
    mean_taxa: totalLoci > 0 ? sumTaxa / totalLoci : 0,
    mean_length: totalLoci > 0 ? sumLength / totalLoci : 0,
    mean_gap_percent: totalLoci > 0 ? sumGap / totalLoci : 0,
    mean_pis: totalLoci > 0 ? sumPis / totalLoci : 0,
    total_matrix_basepairs: totalBasepairs,
  };
}

/** Tooltip for `mean_divergence`. The rule is in `compute_mean_divergence` (stats.rs). */
export const DIVERGENCE_NOTE =
  'Approximate value. For a locus with 40 or more samples, it uses 20 to 30 samples at equal intervals.';

/**
 * Recipe fields that decide only pass and fail. A change in them reuses the
 * measured loci. Keep in step with `GATING_FIELDS` in src-tauri/src/commands.rs.
 */
const GATING_FIELDS = [
  'name',
  'description',
  'assess_alignment',
  'min_taxa',
  'min_taxa_occupancy_percent',
  'min_length',
  'max_gap_percent',
  'min_pis_count',
  'min_pis_percent',
  'min_variable_count',
  'min_variable_percent',
] as const satisfies readonly (keyof TrimmingRecipe)[];

// A reference set is replaced, never edited in place, so each set object gets
// one number, and the key holds that number instead of every sequence.
const referenceSetIds = new WeakMap<object, number>();
let lastReferenceSetId = 0;

function referenceSetId(references: Record<string, string> | undefined): number {
  if (!references || Object.keys(references).length === 0) return 0;
  let id = referenceSetIds.get(references);
  if (id === undefined) {
    id = ++lastReferenceSetId;
    referenceSetIds.set(references, id);
  }
  return id;
}

/** Identifies the recipe settings that change the measured loci. */
export function processingRecipeKey(recipe: TrimmingRecipe): string {
  const fields: Partial<TrimmingRecipe> = { ...recipe };
  for (const field of GATING_FIELDS) delete fields[field];
  delete fields.orf_reference_sequences;
  return `${JSON.stringify(fields)}#${referenceSetId(recipe.orf_reference_sequences)}`;
}

/**
 * The ordinary alignment branch used by the catalog and the general export.
 * ORF extraction is a separate outcome and must not change what the catalog
 * measures or how it decides pass and fail.
 */
export function recipeWithoutOrfAnalysis(recipe: TrimmingRecipe): TrimmingRecipe {
  return {
    ...recipe,
    enable_orf: false,
    orf_use_references: false,
  };
}
