import {
  AlignmentSummary,
  AlignmentViewResponse,
  BatchExportConfig,
  BatchExportResult,
  CatalogUpdateResponse,
  ConcatenateConfig,
  ConcatenateResult,
  GroupedConcatenateConfig,
  GroupedConcatenateResult,
  ParseFailure,
  ScanResponse,
  TrimmingRecipe,
} from './types';
import { buildDatasetOverviewFromSummaries } from './summaries';
import { DEFAULT_RECIPE } from './defaultRecipe';
import { enginePool } from './engine/enginePool';
import type { AlignmentInput } from './engine/engineWorker';

/**
 * Extensions the desktop scanner accepts. The browser file picker must use the
 * same list, otherwise the two builds load different files from one folder.
 * Keep in step with `scan_alignment_directory` in `pipeline/catalog.rs`.
 */
export const SUPPORTED_ALIGNMENT_EXTENSIONS = [
  'fa',
  'fasta',
  'fna',
  'ffn',
  'phy',
  'phylip',
  'nex',
  'nexus',
  'aln',
  'txt',
] as const;

/** Folder levels that the desktop scan reads below the chosen folder. */
const MAX_SCAN_DEPTH = 4;
const SKIPPED_SCAN_NAMES = ['node_modules', 'target', '__MACOSX'];

/**
 * Applies the desktop scan rules to a picked or dropped file: at most
 * `MAX_SCAN_DEPTH` levels below the chosen folder, and no hidden or build
 * folder or file (for example `__MACOSX/._locus.phy`).
 */
function isScannedPath(relativePath: string): boolean {
  const parts = relativePath.split('/').filter(Boolean);
  // The first part is the chosen folder, which the desktop scan always reads.
  const below = parts.length > 1 ? parts.slice(1) : parts;
  return (
    below.length <= MAX_SCAN_DEPTH &&
    below.every((part) => !part.startsWith('.') && !SKIPPED_SCAN_NAMES.includes(part))
  );
}

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

/**
 * Writing files needs the desktop shell. The browser build must refuse instead
 * of reporting an export it cannot perform.
 */
export const DESKTOP_ONLY_EXPORT_MESSAGE =
  'Exporting files needs the AlignmentForge desktop app. The browser version cannot write to your disk. ' +
  'Download the desktop app to export alignments, supermatrices, and gene groups.';

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
  if (!isTauri) {
    throw new Error('Loading a folder by path is available in the desktop app.');
  }
  return invokeTauri<ScanResponse>('scan_directory', { dirPath });
}

export async function loadDirectoryFromFiles(
  files: FileList | File[],
  recipe: TrimmingRecipe,
  onProgress?: (current: number, total: number, fileName: string) => void
): Promise<{ dirName: string; scanResponse: ScanResponse }> {
  const fileArray = Array.from(files).filter((file) => {
    const ext = file.name.split('.').pop()?.toLowerCase() || '';
    return (
      (SUPPORTED_ALIGNMENT_EXTENSIONS as readonly string[]).includes(ext) &&
      isScannedPath(file.webkitRelativePath || file.name)
    );
  });

  const total = fileArray.length;
  if (total === 0) {
    throw new Error(
      `No supported alignment files found in folder. AlignmentForge reads: ${SUPPORTED_ALIGNMENT_EXTENSIONS.map((ext) => `.${ext}`).join(', ')}`
    );
  }

  let dirName = 'Custom Alignments';
  // The files themselves go to the workers. Each worker reads its own shard on
  // its own thread, so reading and parsing run in parallel and the contents are
  // never held twice in memory.
  const inputs: AlignmentInput[] = fileArray.map((file) => {
    if (file.webkitRelativePath) {
      const parts = file.webkitRelativePath.split('/');
      if (parts.length > 1) dirName = parts[0];
    }
    return {
      file,
      id: file.name.replace(/\.[^/.]+$/, ''),
      file_name: file.name,
      file_path: file.webkitRelativePath || file.name,
    };
  });

  const scanResponse = await loadIntoEngine(inputs, recipe, total, onProgress);
  return { dirName, scanResponse };
}

/**
 * Hands the alignments to the worker pool, then builds the scan response from
 * what the engine reports.
 */
async function loadIntoEngine(
  inputs: AlignmentInput[],
  recipe: TrimmingRecipe,
  total: number,
  onProgress?: (current: number, total: number, fileName: string) => void
): Promise<ScanResponse> {
  const loaded = await enginePool.load(inputs, (done, count) => {
    onProgress?.(
      total * ((done / Math.max(1, count)) * 0.75),
      total,
      'Reading and parsing alignments…'
    );
  });

  const [summaries, occupancy] = await Promise.all([
    enginePool.summarize(recipe, (done, count) => {
      onProgress?.(
        total * (0.75 + (done / Math.max(1, count)) * 0.25),
        total,
        'Computing summaries…'
      );
    }),
    enginePool.occupancy(),
  ]);

  return {
    summaries,
    overview: buildDatasetOverviewFromSummaries(summaries, loaded.totalUniqueTaxa),
    occupancy,
    parse_failures: loaded.parseFailures,
  };
}

export async function getAlignment(
  filePath: string,
  recipe: TrimmingRecipe,
  totalUniqueTaxa: number = 0
): Promise<AlignmentViewResponse> {
  if (isTauri) {
    return invokeTauri<AlignmentViewResponse>('get_alignment', { filePath, recipe, totalUniqueTaxa });
  } else {
    // The engine owns the alignments inside the workers, so the view is built
    // by the shard that holds this locus.
    return enginePool.view(filePath, recipe);
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
    // Every worker runs its shard at once, so the main thread stays free and
    // the work spreads across cores instead of blocking on one.
    onProgress?.({ current: 0, total: paths.length, percent: 0, file_name: '' });
    const summaries = await enginePool.summarize(recipe, (done, total) => {
      onProgress?.({
        current: done,
        total,
        percent: total > 0 ? (done / total) * 100 : 100,
        file_name: '',
      });
    });
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
  }
  throw new Error(DESKTOP_ONLY_EXPORT_MESSAGE);
}

export async function runConcatenate(
  config: ConcatenateConfig,
  recipe: TrimmingRecipe
): Promise<ConcatenateResult> {
  if (isTauri) {
    return invokeTauri<ConcatenateResult>('run_concatenate', { config, recipe });
  }
  throw new Error(DESKTOP_ONLY_EXPORT_MESSAGE);
}

export async function runGroupedConcatenate(
  config: GroupedConcatenateConfig,
  recipe: TrimmingRecipe
): Promise<GroupedConcatenateResult> {
  if (isTauri) {
    return invokeTauri<GroupedConcatenateResult>('run_grouped_concatenate', { config, recipe });
  }
  throw new Error(DESKTOP_ONLY_EXPORT_MESSAGE);
}

export async function getPresets(): Promise<TrimmingRecipe[]> {
  if (isTauri) {
    return invokeTauri<TrimmingRecipe[]>('get_presets');
  }
  try {
    // The engine's own presets, so the browser offers the desktop list.
    return await enginePool.presets();
  } catch {
    return [DEFAULT_RECIPE];
  }
}


/**
 * Restricts `?run=` to a folder shipped with the application.
 *
 * The value arrives from the page URL, so a crafted link could otherwise point
 * the app at any host. Only a relative path under `example_data/` is accepted.
 */
export function sanitizeExampleDataPath(rawPath: string): string | null {
  const trimmed = rawPath.trim().replace(/^\/+|\/+$/g, '');
  if (!trimmed) return null;
  // Reject absolute URLs, protocol-relative URLs, and parent traversal.
  if (/^[a-z][a-z0-9+.-]*:/i.test(trimmed)) return null;
  if (trimmed.startsWith('//')) return null;
  if (trimmed.split('/').some((segment) => segment === '..' || segment === '.')) return null;
  if (trimmed !== 'example_data' && !trimmed.startsWith('example_data/')) return null;
  return trimmed;
}

export async function loadDirectoryFromUrl(
  urlPath: string,
  recipe: TrimmingRecipe,
  onProgress?: (current: number, total: number, fileName: string) => void
): Promise<{ dirName: string; scanResponse: ScanResponse }> {
  const safePath = sanitizeExampleDataPath(urlPath);
  if (!safePath) {
    throw new Error(
      `Refused to load '${urlPath}'. The run parameter must name a folder inside example_data.`
    );
  }

  const manifestUrl = `${safePath}/manifest.json`;
  const response = await fetch(manifestUrl);
  if (!response.ok) {
    throw new Error(`Failed to fetch manifest from ${manifestUrl}`);
  }

  const manifest: unknown = await response.json();
  if (!Array.isArray(manifest) || manifest.length === 0) {
    throw new Error('No files found in manifest.');
  }
  // Manifest entries name files inside the same folder and nothing above it.
  const files: string[] = manifest.filter(
    (entry): entry is string =>
      typeof entry === 'string' &&
      entry.length > 0 &&
      !entry.startsWith('/') &&
      !/^[a-z][a-z0-9+.-]*:/i.test(entry) &&
      !entry.split('/').includes('..')
  );
  if (files.length === 0) {
    throw new Error('Manifest contained no usable file names.');
  }

  const dirName = safePath.split('/').filter(Boolean).pop() || 'Example Alignments';
  const total = files.length;
  // Results go into slots by manifest index, so the locus order does not
  // depend on which download finishes first.
  const inputSlots: (AlignmentInput | undefined)[] = new Array(total);
  const failureSlots: (ParseFailure | undefined)[] = new Array(total);
  const BATCH_SIZE = 10;

  for (let i = 0; i < total; i += BATCH_SIZE) {
    const chunk = files.slice(i, i + BATCH_SIZE);
    await Promise.all(
      chunk.map(async (fileName, offset) => {
        const baseName = fileName.split('/').pop() || fileName;
        try {
          const res = await fetch(`${safePath}/${fileName}`);
          if (!res.ok) throw new Error(`Failed to fetch ${fileName}`);
          inputSlots[i + offset] = {
            content: await res.text(),
            id: baseName.replace(/\.[^/.]+$/, ''),
            file_name: baseName,
            file_path: fileName,
          };
        } catch (e) {
          failureSlots[i + offset] = {
            file_name: baseName,
            file_path: fileName,
            error: e instanceof Error ? e.message : String(e),
          };
        }
      })
    );
    onProgress?.(Math.min(total, i + BATCH_SIZE) / 2, total, chunk[chunk.length - 1] || '');
  }

  const inputs = inputSlots.filter((input): input is AlignmentInput => input !== undefined);
  const fetchFailures = failureSlots.filter(
    (failure): failure is ParseFailure => failure !== undefined
  );
  if (inputs.length === 0) {
    throw new Error('None of the manifest files could be downloaded.');
  }

  // Downloads used the first half of the progress bar; parsing uses the second.
  const scanResponse = await loadIntoEngine(inputs, recipe, total, (current, count, fileName) =>
    onProgress?.((total + current) / 2, count, fileName)
  );
  return {
    dirName,
    scanResponse: {
      ...scanResponse,
      parse_failures: [...fetchFailures, ...(scanResponse.parse_failures ?? [])],
    },
  };
}
