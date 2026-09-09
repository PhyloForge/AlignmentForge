import { GeneticCode, TrimmingRecipe } from './types';

/**
 * Display helpers for the viewer.
 *
 * Alignment processing happens in the engine, compiled to WebAssembly. These
 * two functions stay in TypeScript because the canvas needs them per frame,
 * where a round trip to a worker for every codon would not be practical.
 * `npm run parity:codons` checks the translator against the engine's rules.
 */

const STANDARD_AMINO_ACIDS: Record<string, string> = {
  TTT: 'F', TTC: 'F', TTA: 'L', TTG: 'L', TCT: 'S', TCC: 'S', TCA: 'S', TCG: 'S',
  TAT: 'Y', TAC: 'Y', TAA: '*', TAG: '*', TGT: 'C', TGC: 'C', TGA: '*', TGG: 'W',
  CTT: 'L', CTC: 'L', CTA: 'L', CTG: 'L', CCT: 'P', CCC: 'P', CCA: 'P', CCG: 'P',
  CAT: 'H', CAC: 'H', CAA: 'Q', CAG: 'Q', CGT: 'R', CGC: 'R', CGA: 'R', CGG: 'R',
  ATT: 'I', ATC: 'I', ATA: 'I', ATG: 'M', ACT: 'T', ACC: 'T', ACA: 'T', ACG: 'T',
  AAT: 'N', AAC: 'N', AAA: 'K', AAG: 'K', AGT: 'S', AGC: 'S', AGA: 'R', AGG: 'R',
  GTT: 'V', GTC: 'V', GTA: 'V', GTG: 'V', GCT: 'A', GCC: 'A', GCA: 'A', GCG: 'A',
  GAT: 'D', GAC: 'D', GAA: 'E', GAG: 'E', GGT: 'G', GGC: 'G', GGA: 'G', GGG: 'G',
};

/** Uppercases a codon and rewrites uracil as thymine. */
function toDnaCodon(codon: string): string {
  return codon.toUpperCase().replace(/U/g, 'T');
}

export function translateClientCodon(codon: string, geneticCode: GeneticCode): string {
  const normalized = toDnaCodon(codon);
  if (normalized.length !== 3) return 'X';
  if (normalized.includes('-')) return '-';
  if (!/^[ACGT]{3}$/.test(normalized)) return 'X';
  if (geneticCode === 'vertebratemitochondrial') {
    if (normalized === 'AGA' || normalized === 'AGG') return '*';
    if (normalized === 'ATA') return 'M';
    if (normalized === 'TGA') return 'W';
  } else if (geneticCode === 'invertebratemitochondrial') {
    if (normalized === 'AGA' || normalized === 'AGG') return 'S';
    if (normalized === 'ATA') return 'M';
    if (normalized === 'TGA') return 'W';
  }
  return STANDARD_AMINO_ACIDS[normalized] ?? 'X';
}

export function shouldSkipOrfLocus(
  locusId: string,
  searchMode: TrimmingRecipe['orf_search_mode'] = 'continuouscds'
): boolean {
  const lower = locusId.toLowerCase();
  const alwaysSkip =
    lower.startsWith('uce-') ||
    lower.startsWith('uce_') ||
    lower.startsWith('uce') ||
    lower.includes('noncoding') ||
    lower.includes('non-coding') ||
    lower.includes('intergenic') ||
    lower.includes('intron');
  const continuousCdsOnly =
    lower.includes('supercontig') || lower.includes('flanking');
  return alwaysSkip || (searchMode === 'continuouscds' && continuousCdsOnly);
}
