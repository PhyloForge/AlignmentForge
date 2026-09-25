/**
 * Rebuilds the example-data manifests.
 *
 * The browser build cannot list a directory over HTTP, so each example folder
 * ships a manifest.json naming its files. Run this after you add or remove an
 * example alignment: `npm run manifests`.
 */
import { readdirSync, statSync, writeFileSync, readFileSync, existsSync } from 'node:fs';
import { join, resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const DATA_DIR = join(ROOT, 'public/example_data');
const EXTENSIONS = ['.fa', '.fasta', '.fna', '.ffn', '.phy', '.phylip', '.nex', '.nexus'];

const isAlignment = (name) => EXTENSIONS.some((ext) => name.toLowerCase().endsWith(ext));

const checkOnly = process.argv.includes('--check');
let changed = false;

function writeManifest(path, entries) {
  const next = `${JSON.stringify(entries)}\n`;
  const current = existsSync(path) ? readFileSync(path, 'utf8') : null;
  if (current === next) return;
  changed = true;
  if (checkOnly) {
    console.error(`Out of date: ${path.replace(`${ROOT}/`, '')}`);
    return;
  }
  writeFileSync(path, next);
  console.log(`Wrote ${path.replace(`${ROOT}/`, '')} (${entries.length} files)`);
}

const subdirectories = readdirSync(DATA_DIR)
  .filter((entry) => statSync(join(DATA_DIR, entry)).isDirectory())
  .sort();

const combined = [];
for (const directory of subdirectories) {
  const files = readdirSync(join(DATA_DIR, directory)).filter(isAlignment).sort();
  writeManifest(join(DATA_DIR, directory, 'manifest.json'), files);
  combined.push(...files.map((file) => `${directory}/${file}`));
}

writeManifest(join(DATA_DIR, 'manifest.json'), combined);

if (checkOnly && changed) {
  console.error('Example manifests are out of date. Run: npm run manifests');
  process.exit(1);
}
if (!changed) console.log('Example manifests are up to date.');
