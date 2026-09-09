/**
 * Guards the viewer's display helpers against the engine's rules.
 *
 * Alignment processing runs in the engine, but the canvas translates codons
 * itself for speed. An earlier edit left the case and uracil normalisation out,
 * which silently stopped every stop codon from being recognised, so these
 * assertions run as part of `npm run check`.
 */
import { shouldSkipOrfLocus, translateClientCodon } from '../src/sequenceDisplay';

const failures: string[] = [];
const check = (label: string, actual: unknown, expected: unknown) => {
  if (actual !== expected) failures.push(`${label}: expected ${expected}, got ${actual}`);
};

// Lowercase input must translate exactly like uppercase input.
for (const [codon, expected] of [
  ['taa', '*'], ['tag', '*'], ['tga', '*'],
  ['atg', 'M'], ['ggg', 'G'], ['ttt', 'F'],
] as const) {
  check(`lowercase '${codon}'`, translateClientCodon(codon, 'standard'), expected);
  check(`uppercase '${codon}'`, translateClientCodon(codon.toUpperCase(), 'standard'), expected);
}

// RNA spelling must behave like DNA spelling in every code.
for (const code of ['standard', 'vertebratemitochondrial', 'invertebratemitochondrial'] as const) {
  check(`AUG under ${code}`, translateClientCodon('AUG', code), 'M');
  check(`UAA under ${code}`, translateClientCodon('UAA', code), '*');
  check(`aug under ${code}`, translateClientCodon('aug', code), 'M');
}

// Genetic code differences the viewer relies on.
check('TGA standard', translateClientCodon('TGA', 'standard'), '*');
check('TGA vertebrate mito', translateClientCodon('TGA', 'vertebratemitochondrial'), 'W');
check('AGA vertebrate mito', translateClientCodon('AGA', 'vertebratemitochondrial'), '*');
check('AGA invertebrate mito', translateClientCodon('AGA', 'invertebratemitochondrial'), 'S');

// Missing data must never read as a stop.
check('gap codon', translateClientCodon('ta-', 'standard'), '-');
check('ambiguous codon', translateClientCodon('tan', 'standard'), 'X');

// Locus skipping must match `should_skip_orf_locus` in the engine.
check('uce locus skipped', shouldSkipOrfLocus('uce-0044', 'continuouscds'), true);
check('intron locus skipped', shouldSkipOrfLocus('gene_intron_1', 'continuouscds'), true);
check('exon locus kept', shouldSkipOrfLocus('anura-00037_ppan', 'continuouscds'), false);
check('supercontig skipped for continuous', shouldSkipOrfLocus('x_supercontig', 'continuouscds'), true);
check('supercontig kept for shared segment', shouldSkipOrfLocus('x_supercontig', 'bestsharedsegment'), false);

if (failures.length > 0) {
  console.error('Display helper check FAILED:');
  for (const failure of failures) console.error(`  - ${failure}`);
  process.exit(1);
}
console.log('Display helper check OK.');
