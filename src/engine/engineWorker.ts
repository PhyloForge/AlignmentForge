/// <reference lib="webworker" />
/**
 * One shard of the dataset, processed by the Rust engine compiled to WebAssembly.
 *
 * Each worker owns an `EngineSession` holding its own alignments, so no state is
 * shared between workers and the main thread never runs pipeline work.
 */
import init, { EngineSession, enginePresets } from '../engine-wasm/engine.js';

export interface AlignmentInput {
  /** Either the text itself, or a file the worker reads on its own thread. */
  content?: string;
  file?: File;
  id: string;
  file_name: string;
  file_path: string;
}

type Request =
  | { kind: 'init'; module: WebAssembly.Module }
  | { kind: 'load'; requestId: number; alignments: AlignmentInput[] }
  | { kind: 'setDataset'; requestId: number; datasetTaxa: string[][]; totalUniqueTaxa: number }
  | { kind: 'summarize'; requestId: number; recipe: unknown; cacheKey: string }
  | { kind: 'retention'; requestId: number; cacheKey: string }
  | { kind: 'view'; requestId: number; filePath: string; recipe: unknown; includeRaw: boolean }
  | { kind: 'taxonStats'; requestId: number }
  | { kind: 'presets'; requestId: number };

let ready: Promise<void> | null = null;
let session: EngineSession | null = null;

// The pool compiles the engine once and sends the module to every worker.
let receiveEngine: (module: WebAssembly.Module) => void = () => {};
const engineModule = new Promise<WebAssembly.Module>((resolve) => {
  receiveEngine = resolve;
});

function ensureReady(): Promise<void> {
  if (!ready) {
    ready = engineModule
      .then((module) => init({ module_or_path: module }))
      .then(() => {
        session = new EngineSession();
      });
  }
  return ready;
}

/** The newest summarize request. An older run stops at its next chunk. */
let latestSummarize = 0;

// One channel for every yield. Held here, so it stays alive between messages.
const yieldChannel = new MessageChannel();
const waitingYields: (() => void)[] = [];
yieldChannel.port1.onmessage = () => waitingYields.shift()?.();

/** Lets queued messages run, so a newer recipe can reach this worker. */
function yieldToMessages(): Promise<void> {
  return new Promise((resolve) => {
    waitingYields.push(resolve);
    yieldChannel.port2.postMessage(null);
  });
}

self.onmessage = async (event: MessageEvent<Request>) => {
  const request = event.data;
  if (request.kind === 'init') {
    receiveEngine(request.module);
    return;
  }
  // Counted before any await, so a run in progress sees it at its next check.
  const summarizeRun = request.kind === 'summarize' ? ++latestSummarize : 0;
  try {
    await ensureReady();
    if (!session) throw new Error('Engine session was not created');

    switch (request.kind) {
      case 'load': {
        // Parse failures are reported per file so the caller can list them
        // instead of losing the locus silently.
        let loadedInShard = 0;
        const taxa: string[][] = [];
        const failures: { file_name: string; file_path: string; error: string }[] = [];
        for (const alignment of request.alignments) {
          try {
            // Reading happens here so it runs on this thread, in parallel with
            // every other worker, rather than serially on the main thread.
            const content =
              alignment.content ?? (alignment.file ? await alignment.file.text() : '');
            taxa.push(
              session.addAlignment(
                content,
                alignment.id,
                alignment.file_name,
                alignment.file_path
              ) as string[]
            );
          } catch (error) {
            failures.push({
              file_name: alignment.file_name,
              file_path: alignment.file_path,
              error: error instanceof Error ? error.message : String(error),
            });
          }
          loadedInShard++;
          if (loadedInShard % 8 === 0 || loadedInShard === request.alignments.length) {
            self.postMessage({
              progressFor: request.requestId,
              done: loadedInShard,
              total: request.alignments.length,
            });
          }
        }
        self.postMessage({ requestId: request.requestId, ok: true, result: { taxa, failures } });
        break;
      }
      case 'summarize': {
        const total = session.alignmentCount();
        const cached = session.summarizeBegin(request.recipe, request.cacheKey);
        if (!cached) {
          // Walk the shard in ranges so progress moves while the work runs.
          // Nothing is serialised until the end, so ranges stay cheap.
          const chunk = Math.max(1, Math.min(16, Math.ceil(total / 8)));
          for (let start = 0; start < total; start += chunk) {
            session.summarizeChunk(start, chunk);
            self.postMessage({
              progressFor: request.requestId,
              done: Math.min(total, start + chunk),
              total,
            });
            // The session holds one run, so a replaced run must stop before it
            // touches the session again.
            await yieldToMessages();
            if (summarizeRun !== latestSummarize) {
              self.postMessage({
                requestId: request.requestId,
                ok: false,
                error: 'Summaries superseded by a newer recipe',
              });
              return;
            }
          }
        } else {
          self.postMessage({ progressFor: request.requestId, done: total, total });
        }
        self.postMessage({
          requestId: request.requestId,
          ok: true,
          result: session.summarizeFinish(request.cacheKey),
        });
        break;
      }
      case 'retention': {
        self.postMessage({
          requestId: request.requestId,
          ok: true,
          result: session.retentionDetails(request.cacheKey),
        });
        break;
      }
      case 'view': {
        const result = session.viewAlignment(request.filePath, request.recipe, request.includeRaw);
        self.postMessage({ requestId: request.requestId, ok: true, result });
        break;
      }
      case 'setDataset': {
        session.setDataset(request.datasetTaxa, request.totalUniqueTaxa);
        self.postMessage({ requestId: request.requestId, ok: true, result: null });
        break;
      }
      case 'taxonStats': {
        self.postMessage({ requestId: request.requestId, ok: true, result: session.taxonStats() });
        break;
      }
      case 'presets': {
        self.postMessage({ requestId: request.requestId, ok: true, result: enginePresets() });
        break;
      }
    }
  } catch (error) {
    self.postMessage({
      requestId: request.requestId,
      ok: false,
      error: error instanceof Error ? error.message : String(error),
    });
  }
};
