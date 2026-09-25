/**
 * Compares the browser engine (the WebAssembly build) with the desktop engine
 * on the example datasets.
 *
 * The two builds share the Rust source, but the browser reads text, splits each
 * run across several calls, and converts values through wasm-bindgen. This
 * check runs both builds on the same files and recipes, and compares every
 * field of every locus summary. Run `npm run build:wasm` first.
 */
import { execFileSync } from 'node:child_process';
import { existsSync, readFileSync } from 'node:fs';
import { basename, dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const DATA_DIR = join(ROOT, 'public/example_data');
const DATASETS = ['exons', 'uces', 'all_markers'];
const WASM_PATH = join(ROOT, 'src/engine-wasm/engine_bg.wasm');

if (!existsSync(WASM_PATH)) {
  console.error('No WebAssembly build found. Run `npm run build:wasm` first.');
  process.exit(1);
}
const { EngineSession, enginePresets, initSync } = await import('../src/engine-wasm/engine.js');
initSync({ module: readFileSync(WASM_PATH) });

const presets = enginePresets();
const recipes = [
  ...presets,
  { ...presets[0], name: 'Default with ContinuousCds ORF', enable_orf: true },
  {
    ...presets[0],
    name: 'Default with BestSharedSegment ORF',
    enable_orf: true,
    orf_search_mode: 'bestsharedsegment',
  },
];

/** Runs the recipes through the calls that a browser worker makes. */
function browserSummaries(files) {
  const session = new EngineSession();
  const datasetTaxa = files.map((file) => {
    const name = basename(file);
    // The locus name rule of src/tauriClient.ts.
    const id = name.replace(/\.[^/.]+$/, '');
    return session.addAlignment(readFileSync(file, 'utf8'), id, name, file);
  });
  session.setDataset(datasetTaxa, new Set(datasetTaxa.flat()).size);

  return recipes.map((recipe, index) => {
    const cacheKey = `recipe-${index}`;
    session.summarizeBegin(recipe, cacheKey);
    session.summarizeChunk(0, files.length);
    const summaries = session.summarizeFinish(cacheKey);
    // The browser fetches the per-sample lists with a separate call.
    const retention = new Map(
      session.retentionDetails(cacheKey).map((detail) => [detail.file_path, detail])
    );
    return summaries.map((summary) => ({ ...summary, ...retention.get(summary.file_path) }));
  });
}

/** Runs the recipes through the desktop engine (src-tauri/examples/engine_summaries.rs). */
function desktopSummaries(files) {
  const output = execFileSync(
    'cargo',
    [
      'run',
      '--quiet',
      '--manifest-path',
      join(ROOT, 'src-tauri/Cargo.toml'),
      '--example',
      'engine_summaries',
    ],
    {
      input: JSON.stringify({ files, recipes }),
      encoding: 'utf8',
      maxBuffer: 1024 * 1024 * 1024,
      stdio: ['pipe', 'pipe', 'inherit'],
    }
  );
  return JSON.parse(output);
}

// JSON has no NaN, so a NaN on the desktop arrives as null. Show both as they are.
const show = (value) => (typeof value === 'number' ? String(value) : JSON.stringify(value));

function collectDifferences(desktop, browser, where, differences) {
  if (Array.isArray(desktop) && Array.isArray(browser) && desktop.length === browser.length) {
    desktop.forEach((value, index) =>
      collectDifferences(value, browser[index], `${where}[${index}]`, differences)
    );
  } else if (
    desktop !== null &&
    browser !== null &&
    typeof desktop === 'object' &&
    typeof browser === 'object' &&
    !Array.isArray(desktop) &&
    !Array.isArray(browser)
  ) {
    for (const key of new Set([...Object.keys(desktop), ...Object.keys(browser)])) {
      collectDifferences(desktop[key], browser[key], `${where}.${key}`, differences);
    }
  } else if (!Object.is(desktop, browser)) {
    differences.push(`${where}: desktop ${show(desktop)}, browser ${show(browser)}`);
  }
}

const differences = [];
let comparedSummaries = 0;
for (const dataset of DATASETS) {
  const directory = join(DATA_DIR, dataset);
  const files = JSON.parse(readFileSync(join(directory, 'manifest.json'), 'utf8')).map((name) =>
    join(directory, name)
  );
  const desktop = desktopSummaries(files);
  const browser = browserSummaries(files);

  recipes.forEach((recipe, index) => {
    const label = `${dataset} / ${recipe.name}`;
    if (desktop[index].length !== browser[index].length) {
      differences.push(
        `${label}: desktop has ${desktop[index].length} loci, browser has ${browser[index].length}`
      );
      return;
    }
    desktop[index].forEach((summary, locus) => {
      collectDifferences(summary, browser[index][locus], `${label} / ${summary.id}`, differences);
    });
    comparedSummaries += desktop[index].length;
  });
}

if (differences.length > 0) {
  console.error(`The browser and desktop engines differ in ${differences.length} values:`);
  for (const line of differences.slice(0, 20)) console.error(`  ${line}`);
  process.exit(1);
}
console.log(
  `Engine parity: ${comparedSummaries} locus summaries match (${DATASETS.length} datasets, ${recipes.length} recipes).`
);
