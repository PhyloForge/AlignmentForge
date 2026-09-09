/**
 * Regenerates the desktop application icons from the logo.
 *
 * The icons under `src-tauri/icons` are committed artifacts, so this only needs
 * running when the logo changes. It rasterises `public/logo.svg` and hands the
 * result to Tauri's icon generator, which produces every size the bundlers
 * need, including the macOS `.icns` and Windows `.ico` containers.
 */
import { execFileSync } from 'node:child_process';
import { mkdtempSync, existsSync, rmSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';

const SOURCE = 'public/logo.svg';
const SIZE = 1024;

if (!existsSync(SOURCE)) {
  console.error(`Cannot find ${SOURCE}`);
  process.exit(1);
}

if (process.platform !== 'darwin') {
  console.error(
    'This script rasterises the logo with qlmanage, which is macOS only.\n' +
      'On another system, export public/logo.svg to a 1024x1024 PNG with any\n' +
      'tool, then run: npx tauri icon <that file> --output src-tauri/icons'
  );
  process.exit(1);
}

const work = mkdtempSync(join(tmpdir(), 'af-icons-'));
try {
  execFileSync('qlmanage', ['-t', '-s', String(SIZE), '-o', work, SOURCE], {
    stdio: 'ignore',
  });
  const raster = join(work, 'logo.svg.png');
  if (!existsSync(raster)) {
    throw new Error('qlmanage did not produce a PNG');
  }
  execFileSync('npx', ['tauri', 'icon', raster, '--output', 'src-tauri/icons'], {
    stdio: 'inherit',
  });
  console.log('\nDesktop icons regenerated from', SOURCE);
} finally {
  rmSync(work, { recursive: true, force: true });
}
