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
