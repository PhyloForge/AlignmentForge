/**
 * A pool of Web Workers, each running the Rust engine compiled to WebAssembly.
 *
 * The browser used to re-implement the trimming pipeline in TypeScript, which
 * had to be kept in step with the engine by hand. This runs the engine itself,
 * so the browser and the desktop app cannot disagree, and it spreads the loci
 * across cores instead of blocking the interface on one thread.
 */
import type {
  AlignmentSummary,
  AlignmentViewResponse,
  ParseFailure,
  TaxonOccupancy,
  TrimmingRecipe,
} from '../types';
import type { AlignmentInput } from './engineWorker';

interface Pending {
  resolve: (value: unknown) => void;
  reject: (reason: Error) => void;
}

class EngineWorker {
  private readonly worker: Worker;
  private readonly pending = new Map<number, Pending>();
  private nextRequestId = 1;
  /** File paths this worker holds, so a view request reaches the right shard. */
  readonly filePaths = new Set<string>();
  /** Reports partial completion of an in-flight request. */
  onProgress?: (requestId: number, done: number, total: number) => void;

  constructor() {
    this.worker = new Worker(new URL('./engineWorker.ts', import.meta.url), {
      type: 'module',
    });
    this.worker.onmessage = (event) => {
      const { requestId, ok, result, error, progressFor, done, total } = event.data;
      if (progressFor !== undefined) {
        this.onProgress?.(progressFor, done, total);
        return;
      }
      const entry = this.pending.get(requestId);
      if (!entry) return;
      this.pending.delete(requestId);
      if (ok) entry.resolve(result);
      else entry.reject(new Error(error));
    };
    this.worker.onerror = (event) => {
      const failure = new Error(event.message || 'Engine worker failed');
      for (const entry of this.pending.values()) entry.reject(failure);
      this.pending.clear();
    };
  }

  send<T>(request: Record<string, unknown>): Promise<T> {
    const requestId = this.nextRequestId++;
    return new Promise<T>((resolve, reject) => {
      this.pending.set(requestId, { resolve: resolve as (value: unknown) => void, reject });
      this.worker.postMessage({ ...request, requestId });
    });
  }

  terminate() {
    this.worker.terminate();
    this.pending.clear();
  }
}

export interface LoadedDataset {
  datasetTaxa: string[][];
  totalUniqueTaxa: number;
  parseFailures: ParseFailure[];
  /** Order in which loci are distributed, so summaries can be reassembled. */
  orderedPaths: string[];
}

export interface RetentionDetail {
  file_path: string;
  retained_taxa: string[];
  retained_taxon_basepairs: Record<string, number>;
  orf_retained_taxa: string[];
  orf_retained_taxon_basepairs: Record<string, number>;
}

/**
 * Identifies a recipe for the worker-side result cache.
 *
 * Only the fields that change processing are included. The pass thresholds are
 * left out because the engine caches the measured loci and applies gating on
 * top of them, so moving a threshold reuses the cached measurements.
 */
function recipeCacheKey(recipe: TrimmingRecipe): string {
  const fields: Record<string, unknown> = { ...recipe };
  for (const gatingField of [
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
  ]) {
    delete fields[gatingField];
  }
  return JSON.stringify(fields);
}

export class EnginePool {
  private workers: EngineWorker[] = [];
  private dataset: LoadedDataset | null = null;
  private lastCacheKey: string | null = null;

  get size(): number {
    return this.workers.length;
  }

  private ensureWorkers(count: number) {
    if (this.workers.length >= count) return;
    while (this.workers.length < count) this.workers.push(new EngineWorker());
  }

  /** Spreads alignments across workers and parses them there. */
  async load(
    alignments: AlignmentInput[],
    onProgress?: (loaded: number, total: number) => void
  ): Promise<LoadedDataset> {
    this.terminate();

    // One worker per core, capped by the number of loci so small datasets do
    // not pay for workers that would sit idle.
    const cores = Math.max(1, (navigator.hardwareConcurrency || 4) - 1);
    const workerCount = Math.max(1, Math.min(cores, alignments.length));
    this.ensureWorkers(workerCount);

    const shards: AlignmentInput[][] = Array.from({ length: workerCount }, () => []);
    alignments.forEach((alignment, index) => shards[index % workerCount].push(alignment));

    const doneByWorker = new Map<EngineWorker, number>();
    const reportLoaded = () => {
      let completed = 0;
      for (const value of doneByWorker.values()) completed += value;
      onProgress?.(Math.min(completed, alignments.length), alignments.length);
    };

    const perShard = await Promise.all(
      shards.map(async (shard, index) => {
        const worker = this.workers[index];
        worker.onProgress = (_requestId, done) => {
          doneByWorker.set(worker, done);
          reportLoaded();
        };
        const result = await worker
          .send<{ taxa: string[][]; failures: ParseFailure[] }>({
            kind: 'load',
            alignments: shard,
          })
          .finally(() => {
            worker.onProgress = undefined;
          });
        for (const alignment of shard) worker.filePaths.add(alignment.file_path);
        doneByWorker.set(worker, shard.length);
        reportLoaded();
        return { shard, ...result };
      })
    );

    const datasetTaxa: string[][] = [];
    const orderedPaths: string[] = [];
    const parseFailures: ParseFailure[] = [];
    for (const entry of perShard) {
      parseFailures.push(...entry.failures);
      const failedPaths = new Set(entry.failures.map((failure) => failure.file_path));
      const parsed = entry.shard.filter((alignment) => !failedPaths.has(alignment.file_path));
      parsed.forEach((alignment, index) => {
        orderedPaths.push(alignment.file_path);
        datasetTaxa.push(entry.taxa[index] ?? []);
      });
    }

    const totalUniqueTaxa = new Set(datasetTaxa.flat()).size;
    // Send the dataset-wide taxon lists once. Every later call reuses them
    // inside the worker instead of copying them across on each recipe change.
    await Promise.all(
      this.workers.map((worker) =>
        worker.send({ kind: 'setDataset', datasetTaxa, totalUniqueTaxa })
      )
    );

    this.dataset = { datasetTaxa, totalUniqueTaxa, parseFailures, orderedPaths };
    return this.dataset;
  }

  /**
   * Summarises every loaded locus, in parallel across the pool.
   *
   * Workers report partial completion as they go, and those counts are summed
   * so the caller sees one moving figure for the whole dataset.
   */
  async summarize(
    recipe: TrimmingRecipe,
    onProgress?: (done: number, total: number) => void
  ): Promise<AlignmentSummary[]> {
    if (!this.dataset) return [];
    const cacheKey = recipeCacheKey(recipe);
    this.lastCacheKey = cacheKey;
    const totalLoci = this.dataset.orderedPaths.length;
    const doneByWorker = new Map<EngineWorker, number>();

    const results = await Promise.all(
      this.workers.map((worker) => {
        worker.onProgress = (_requestId, done) => {
          doneByWorker.set(worker, done);
          let completed = 0;
          for (const value of doneByWorker.values()) completed += value;
          onProgress?.(Math.min(completed, totalLoci), totalLoci);
        };
        return worker
          .send<AlignmentSummary[]>({ kind: 'summarize', recipe, cacheKey })
          .finally(() => {
            worker.onProgress = undefined;
          });
      })
    );
    onProgress?.(totalLoci, totalLoci);
    return results.flat();
  }

  /** Full viewer payload for one locus, from the worker that holds it. */
  async view(filePath: string, recipe: TrimmingRecipe): Promise<AlignmentViewResponse> {
    if (!this.dataset) throw new Error('No dataset is loaded');
    const worker = this.workers.find((candidate) => candidate.filePaths.has(filePath));
    if (!worker) throw new Error(`Alignment not found: ${filePath}`);
    return worker.send<AlignmentViewResponse>({ kind: 'view', filePath, recipe });
  }

  /** Dataset-wide taxon occupancy, merged from every shard. */
  async occupancy(): Promise<TaxonOccupancy[]> {
    if (!this.dataset) return [];
    const totalLoci = this.dataset.orderedPaths.length;
    const shards = await Promise.all(
      this.workers.map((worker) =>
        worker.send<{ taxon: string; basepairs: number; gap_percent: number }[]>({
          kind: 'taxonStats',
        })
      )
    );

    const merged = new Map<string, { count: number; totalBp: number; sumGap: number }>();
    for (const stat of shards.flat()) {
      const entry = merged.get(stat.taxon) ?? { count: 0, totalBp: 0, sumGap: 0 };
      entry.count++;
      entry.totalBp += stat.basepairs;
      entry.sumGap += stat.gap_percent;
      merged.set(stat.taxon, entry);
    }

    const occupancy: TaxonOccupancy[] = Array.from(merged.entries()).map(([taxon_name, stats]) => ({
      taxon_name,
      present_loci_count: stats.count,
      present_loci_percent: (stats.count / Math.max(1, totalLoci)) * 100,
      mean_gap_percent: stats.sumGap / Math.max(1, stats.count),
      total_bp: stats.totalBp,
    }));
    occupancy.sort((a, b) => b.present_loci_count - a.present_loci_count);
    return occupancy;
  }

  /**
   * Per-sample retention for the QC occupancy chart.
   *
   * Kept out of the ordinary summaries because only that one chart needs it,
   * and it carries a taxon list for every locus.
   */
  async retentionDetails(): Promise<Map<string, RetentionDetail>> {
    const merged = new Map<string, RetentionDetail>();
    if (!this.dataset || !this.lastCacheKey) return merged;
    const shards = await Promise.all(
      this.workers.map((worker) =>
        worker
          .send<RetentionDetail[]>({ kind: 'retention', cacheKey: this.lastCacheKey })
          .catch(() => [] as RetentionDetail[])
      )
    );
    for (const detail of shards.flat()) merged.set(detail.file_path, detail);
    return merged;
  }

  async presets(): Promise<TrimmingRecipe[]> {
    this.ensureWorkers(1);
    return this.workers[0].send<TrimmingRecipe[]>({ kind: 'presets' });
  }

  hasDataset(): boolean {
    return this.dataset !== null && this.dataset.orderedPaths.length > 0;
  }

  terminate() {
    for (const worker of this.workers) worker.terminate();
    this.workers = [];
    this.dataset = null;
    this.lastCacheKey = null;
  }
}

export const enginePool = new EnginePool();
