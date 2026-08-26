import {
  Alignment,
  AlignmentSummary,
  AlignmentViewResponse,
  BatchExportConfig,
  BatchExportResult,
  CatalogUpdateResponse,
  ConcatenateConfig,
  ConcatenateResult,
  GroupedConcatenateConfig,
  GroupedConcatenateResult,
  DatasetOverview,
  ScanResponse,
  TaxonOccupancy,
  TrimmingRecipe,
} from './types';
import {
  buildDatasetOverviewFromSummaries,
  buildScanResponseFromAlignments,
  computeAlignmentSummary,
  executeClientTrimming,
  parseAlignmentText,
  recipeWithDatasetSampleFilter,
} from './parsers/clientParser';

// In-memory cache for client-side / browser loaded alignments
const clientAlignmentsCache = new Map<string, Alignment>();
let catalogJobCounter = 0;

export interface CatalogRecalculationProgress {
  current: number;
  total: number;
  percent: number;
  file_name: string;
}

// Detect if running inside Tauri desktop shell
export const isTauri =
  typeof window !== 'undefined' &&
  ('__TAURI_INTERNALS__' in window || '__TAURI__' in window);

async function invokeTauri<T>(cmd: string, args?: Record<string, unknown>): Promise<T> {
  if (isTauri) {
    const { invoke } = await import('@tauri-apps/api/core');
    return invoke<T>(cmd, args);
  }
  throw new Error('Not running inside Tauri');
}

export async function openDirectoryDialog(): Promise<string | null> {
  if (isTauri) {
    try {
      const { open } = await import('@tauri-apps/plugin-dialog');
      const selected = await open({
        directory: true,
        multiple: false,
        title: 'Select Alignments Directory',
      });
      return (selected as string) || null;
    } catch (err) {
      console.warn('Native dialog open error:', err);
      return null;
    }
  }
  return null;
}

export async function openFileDialog(): Promise<string | null> {
  if (isTauri) {
    try {
      const { open } = await import('@tauri-apps/plugin-dialog');
      const selected = await open({
        directory: false,
        multiple: false,
        title: 'Select Gene Mapping File (CSV/TSV/TXT)',
      });
      return (selected as string) || null;
    } catch (err) {
      console.warn('Native file dialog error:', err);
      return null;
    }
  }
  return null;
}
export async function openSaveDirectoryDialog(): Promise<string | null> {
  if (isTauri) {
    try {
      const { open } = await import('@tauri-apps/plugin-dialog');
      const selected = await open({
        directory: true,
        multiple: false,
        title: 'Select Export Destination Folder',
      });
      return (selected as string) || null;
    } catch (err) {
      console.warn('Native save dialog error:', err);
      return null;
    }
  }
  return 'output/trimmed_alignments';
}

export interface LoadedFilterConfig {
  filePath: string;
  recipe: TrimmingRecipe;
}

export async function exportFilterConfig(recipe: TrimmingRecipe): Promise<string | null> {
  if (!isTauri) {
    throw new Error('Filter configuration export is available in the desktop app.');
  }

  const { save } = await import('@tauri-apps/plugin-dialog');
  const safeName = recipe.name
    .trim()
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-+|-+$/g, '') || 'alignmentforge-filters';
  const filePath = await save({
    title: 'Export AlignmentForge Filter Configuration',
    defaultPath: `${safeName}.toml`,
    filters: [{ name: 'AlignmentForge Filter Config', extensions: ['toml'] }],
  });

  if (!filePath) return null;
  return invokeTauri<string>('save_filter_config', { filePath, recipe });
}

export async function loadFilterConfig(): Promise<LoadedFilterConfig | null> {
  if (!isTauri) {
    throw new Error('Filter configuration loading is available in the desktop app.');
  }

  const { open } = await import('@tauri-apps/plugin-dialog');
  const selected = await open({
    title: 'Load AlignmentForge Filter Configuration',
    directory: false,
    multiple: false,
    filters: [{ name: 'AlignmentForge Filter Config', extensions: ['toml'] }],
  });

  if (!selected || Array.isArray(selected)) return null;
  const recipe = await invokeTauri<TrimmingRecipe>('load_filter_config', { filePath: selected });
  return { filePath: selected, recipe };
}

export async function exportAlignmentStatsCsv(
  summaries: AlignmentSummary[]
): Promise<string | null> {
  if (!isTauri) {
    throw new Error('Alignment statistics CSV export is available in the desktop app.');
  }

  const { save } = await import('@tauri-apps/plugin-dialog');
  const filePath = await save({
    title: 'Export Alignment Statistics',
    defaultPath: 'alignmentforge-alignment-stats.csv',
    filters: [{ name: 'CSV Spreadsheet', extensions: ['csv'] }],
  });

  if (!filePath) return null;
  return invokeTauri<string>('save_alignment_stats_csv', { filePath, summaries });
}

export async function scanDirectory(dirPath: string): Promise<ScanResponse> {
  if (isTauri) {
    return invokeTauri<ScanResponse>('scan_directory', { dirPath });
  } else {
    if (clientAlignmentsCache.size > 0) {
      const aligns = Array.from(clientAlignmentsCache.values());
      return await buildScanResponseFromAlignments(aligns, defaultRecipe);
    }
    return {
      summaries: [],
      overview: {
        total_alignments: 0,
        passed_alignments: 0,
        discarded_alignments: 0,
        total_unique_taxa: 0,
        mean_taxa: 0,
        mean_length: 0,
        mean_gap_percent: 0,
        mean_pis: 0,
        total_matrix_basepairs: 0,
      },
      occupancy: [],
    };
  }
}

export async function loadDirectoryFromFiles(
  files: FileList | File[],
  recipe: TrimmingRecipe,
  onProgress?: (current: number, total: number, fileName: string) => void
): Promise<{ dirName: string; scanResponse: ScanResponse }> {
  clientAlignmentsCache.clear();
  const fileArray = Array.from(files).filter((file) => {
    const ext = file.name.split('.').pop()?.toLowerCase() || '';
    return ['fa', 'fasta', 'fna', 'faa', 'phy', 'phylip', 'nex', 'nexus'].includes(ext);
  });

  const total = fileArray.length;
  if (total === 0) {
    throw new Error('No supported alignment files (.fa, .fasta, .phy, .nex) found in folder.');
  }

  let dirName = 'Custom Alignments';
  const parsedAlignments: Alignment[] = [];
  const BATCH_SIZE = 30; // process in small non-blocking chunks

  for (let i = 0; i < total; i += BATCH_SIZE) {
    const chunk = fileArray.slice(i, i + BATCH_SIZE);

    for (const file of chunk) {
      if (file.webkitRelativePath) {
        const parts = file.webkitRelativePath.split('/');
        if (parts.length > 1) {
          dirName = parts[0];
        }
      }

      try {
        const text = await file.text();
        const align = parseAlignmentText(text, file.name, file.webkitRelativePath || file.name);
        parsedAlignments.push(align);
        clientAlignmentsCache.set(align.file_path, align);
      } catch (e) {
        console.warn('Failed to parse alignment:', file.name, e);
      }
    }

    if (onProgress) {
      const current = Math.min(total, i + BATCH_SIZE);
      const lastFile = chunk[chunk.length - 1]?.name || '';
      onProgress(current, total, lastFile);
    }

    // Yield control to browser renderer loop
    await new Promise((resolve) => setTimeout(resolve, 1));
  }

  const scanResponse = await buildScanResponseFromAlignments(parsedAlignments, recipe, 0, (pct) => {
    if (onProgress) onProgress(total * (pct / 100), total, 'Computing summaries...');
  });
  return { dirName, scanResponse };
}

export async function getAlignment(
  filePath: string,
  recipe: TrimmingRecipe,
  totalUniqueTaxa: number = 0
): Promise<AlignmentViewResponse> {
  if (isTauri) {
    return invokeTauri<AlignmentViewResponse>('get_alignment', { filePath, recipe, totalUniqueTaxa });
  } else {
    const cached = clientAlignmentsCache.get(filePath);
    if (cached) {
      const runtimeRecipe = recipeWithDatasetSampleFilter(
        recipe,
        Array.from(clientAlignmentsCache.values())
      );
      return executeClientTrimming(cached, runtimeRecipe, totalUniqueTaxa);
    }
    throw new Error(`Alignment not found: ${filePath}`);
  }
}

export async function recalculateCatalog(
  paths: string[],
  recipe: TrimmingRecipe,
  totalUniqueTaxa: number = 0,
  onProgress?: (progress: CatalogRecalculationProgress) => void
): Promise<CatalogUpdateResponse> {
  if (isTauri) {
    const jobId = `catalog-${Date.now()}-${++catalogJobCounter}`;
    const { listen } = await import('@tauri-apps/api/event');
    const unlisten = await listen<CatalogRecalculationProgress & { job_id: string }>(
      'catalog_recalculation_progress',
      (event) => {
        if (event.payload.job_id === jobId) {
          onProgress?.(event.payload);
        }
      }
    );

    onProgress?.({ current: 0, total: paths.length, percent: 0, file_name: '' });
    try {
      return await invokeTauri<CatalogUpdateResponse>('recalculate_catalog', {
        paths,
        recipe,
        totalUniqueTaxa,
        jobId,
      });
    } finally {
      unlisten();
    }
  } else {
    const aligns = paths
      .map((p) => clientAlignmentsCache.get(p))
      .filter((a): a is Alignment => a !== undefined);
    const summaries: AlignmentSummary[] = [];
    const batchSize = 10;
    const runtimeRecipe = recipeWithDatasetSampleFilter(recipe, aligns);

    onProgress?.({ current: 0, total: aligns.length, percent: 0, file_name: '' });
    for (let start = 0; start < aligns.length; start += batchSize) {
      const batch = aligns.slice(start, start + batchSize);
      for (const alignment of batch) {
        summaries.push(computeAlignmentSummary(alignment, runtimeRecipe, totalUniqueTaxa));
      }
      const current = Math.min(aligns.length, start + batch.length);
      onProgress?.({
        current,
        total: aligns.length,
        percent: aligns.length > 0 ? (current / aligns.length) * 100 : 100,
        file_name: batch[batch.length - 1]?.file_name || '',
      });
      await new Promise((resolve) => setTimeout(resolve, 0));
    }

    return {
      summaries,
      overview: buildDatasetOverviewFromSummaries(summaries, totalUniqueTaxa),
    };
  }
}

export async function runBatchExport(
  config: BatchExportConfig,
  recipe: TrimmingRecipe
): Promise<BatchExportResult> {
  if (isTauri) {
    return invokeTauri<BatchExportResult>('run_batch_export', { config, recipe });
  } else {
    return {
      total_processed: config.input_paths.length,
      total_exported: config.export_general_alignments ? config.input_paths.length : 0,
      total_discarded: 0,
      total_orfs_exported:
        config.export_orf_alignments && recipe.enable_orf ? config.input_paths.length : 0,
      alignment_directory_path: config.export_general_alignments
        ? `${config.output_directory}/${config.general_alignment_directory_name}`
        : undefined,
      orf_directory_path:
        config.export_orf_alignments && recipe.enable_orf
          ? `${config.output_directory}/${config.orf_alignment_directory_name}`
          : undefined,
      summary_csv_path: `${config.output_directory}/alignment-trimming_summary.csv`,
      recipe_json_path: `${config.output_directory}/recipe.json`,
      total_introns_exported: config.export_introns ? config.input_paths.length : 0,
      intron_directory_path: config.export_introns
        ? `${config.output_directory}/${config.intron_directory_name}`
        : undefined,
    };
  }
}

export async function runConcatenate(
  config: ConcatenateConfig,
  recipe: TrimmingRecipe
): Promise<ConcatenateResult> {
  if (isTauri) {
    return invokeTauri<ConcatenateResult>('run_concatenate', { config, recipe });
  } else {
    return {
      total_taxa: 20,
      total_length: 12000,
      total_loci: config.input_paths.length,
      supermatrix_path: `${config.output_file_prefix}.${config.output_format === 'phylip' ? 'phy' : 'fa'}`,
    };
  }
}

export async function runGroupedConcatenate(
  config: GroupedConcatenateConfig,
  recipe: TrimmingRecipe
): Promise<GroupedConcatenateResult> {
  if (isTauri) {
    return invokeTauri<GroupedConcatenateResult>('run_grouped_concatenate', { config, recipe });
  } else {
    return {
      total_genes: 5,
      total_exons_processed: config.input_paths.length,
      output_directory: config.output_directory,
    };
  }
}

export async function getPresets(): Promise<TrimmingRecipe[]> {
  if (isTauri) {
    return invokeTauri<TrimmingRecipe[]>('get_presets');
  } else {
    return [defaultRecipe];
  }
}

const defaultRecipe: TrimmingRecipe = {
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

export async function loadDirectoryFromUrl(
  urlPath: string,
  recipe: TrimmingRecipe,
  onProgress?: (current: number, total: number, fileName: string) => void
): Promise<{ dirName: string; scanResponse: ScanResponse }> {
  clientAlignmentsCache.clear();
  
  // Fetch manifest
  const manifestUrl = `${urlPath}/manifest.json`;
  const response = await fetch(manifestUrl);
  if (!response.ok) {
    throw new Error(`Failed to fetch manifest from ${manifestUrl}`);
  }
  
  const files: string[] = await response.json();
  if (!files || files.length === 0) {
    throw new Error('No files found in manifest.');
  }

  let dirName = urlPath.split('/').filter(Boolean).pop() || 'Example Alignments';
  const parsedAlignments: Alignment[] = [];
  const BATCH_SIZE = 10;
  const total = files.length;

  for (let i = 0; i < total; i += BATCH_SIZE) {
    const chunk = files.slice(i, i + BATCH_SIZE);

    const fetchPromises = chunk.map(async (fileName) => {
      try {
        const fileUrl = `${urlPath}/${fileName}`;
        const res = await fetch(fileUrl);
        if (!res.ok) throw new Error(`Failed to fetch ${fileName}`);
        const text = await res.text();
        const baseName = fileName.split('/').pop() || fileName;
        const align = parseAlignmentText(text, baseName, fileName);
        parsedAlignments.push(align);
        clientAlignmentsCache.set(align.file_path, align);
      } catch (e) {
        console.warn('Failed to parse alignment from URL:', fileName, e);
      }
    });
    
    await Promise.all(fetchPromises);

    if (onProgress) {
      const current = Math.min(total, i + BATCH_SIZE);
      const lastFile = chunk[chunk.length - 1] || '';
      onProgress(current, total, lastFile);
    }

    await new Promise((resolve) => setTimeout(resolve, 10));
  }

  const scanResponse = await buildScanResponseFromAlignments(parsedAlignments, recipe, 0, (pct) => {
    if (onProgress) onProgress(total * (pct / 100), total, 'Computing summaries...');
  });
  return { dirName, scanResponse };
}
