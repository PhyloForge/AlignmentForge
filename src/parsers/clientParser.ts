import {
  Alignment,
  AlignmentFormat,
  AlignmentSummary,
  AlignmentViewResponse,
  DatasetOverview,
  GeneticCode,
  ScanResponse,
  StopCodonPos,
  TaxonOccupancy,
  TrimmingDiff,
  TrimmingRecipe,
} from '../types';

function clientStopCodonsForCode(geneticCode: GeneticCode): Set<string> {
  if (geneticCode === 'vertebratemitochondrial') {
    return new Set(['TAA', 'TAG', 'AGA', 'AGG']);
  }
  if (geneticCode === 'invertebratemitochondrial') {
    return new Set(['TAA', 'TAG']);
  }
  return new Set(['TAA', 'TAG', 'TGA', 'UAA', 'UAG', 'UGA']);
}

function reverseComplementClient(sequence: string): string {
  const complements: Record<string, string> = {
    A: 'T', T: 'A', U: 'A', C: 'G', G: 'C',
    R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K',
    B: 'V', D: 'H', H: 'D', V: 'B', N: 'N', '-': '-', '?': '?',
  };
  return sequence
    .toUpperCase()
    .split('')
    .reverse()
    .map((base) => complements[base] ?? base)
    .join('');
}

interface ClientReadingFrame {
  isReverse: boolean;
  offset: number;
  cleanTaxaCount: number;
  score: number;
}

function evaluateClientReadingFrame(
  sequences: string[],
  geneticCode: GeneticCode,
  isReverse: boolean,
  offset: number
): ClientReadingFrame {
  const stops = clientStopCodonsForCode(geneticCode);
  const oriented = isReverse ? sequences.map(reverseComplementClient) : sequences;
  let cleanCodons = 0;
  let stopTaxa = 0;
  let cleanTaxaCount = 0;
  let startCodons = 0;
  let terminalStops = 0;
  for (const sequence of oriented) {
    const codonCount = Math.floor((sequence.length - offset) / 3);
    if (codonCount <= 0) continue;
    let hasInternalStop = false;
    for (let codonIndex = 0; codonIndex < codonCount; codonIndex++) {
      const start = offset + codonIndex * 3;
      const codon = sequence.slice(start, start + 3).toUpperCase();
      if (codonIndex === 0 && codon === 'ATG') startCodons++;
      if (stops.has(codon)) {
        if (codonIndex + 1 === codonCount) terminalStops++;
        else hasInternalStop = true;
      }
    }
    if (hasInternalStop) stopTaxa++;
    else {
      cleanTaxaCount++;
      cleanCodons += codonCount;
    }
  }
  return {
    isReverse,
    offset,
    cleanTaxaCount,
    score: cleanCodons * 100 + startCodons * 250 + terminalStops * 150 - stopTaxa * 80,
  };
}

function findClientReadingFrame(
  sequences: string[],
  geneticCode: GeneticCode
): ClientReadingFrame {
  let best = evaluateClientReadingFrame(sequences, geneticCode, false, 0);
  for (const isReverse of [false, true]) {
    for (let offset = 0; offset < 3; offset++) {
      const candidate = evaluateClientReadingFrame(sequences, geneticCode, isReverse, offset);
      if (candidate.score > best.score) best = candidate;
    }
  }
  return best;
}

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

export function translateClientCodon(codon: string, geneticCode: GeneticCode): string {
  const normalized = codon.toUpperCase().replace(/U/g, 'T');
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

interface ClientSharedOrfSegment {
  frame: number;
  isReverse: boolean;
  start: number;
  end: number;
  supportCount: number;
  lengthCodons: number;
  informativeCodons: number;
  codingScore: number;
  aminoAcidConservation: number;
  frameContrast: number;
}

function clientProteinProfileConservation(
  sequences: string[],
  offset: number,
  geneticCode: GeneticCode
): number {
  const length = Math.min(...sequences.map((sequence) => sequence.length));
  const codonCount = Math.floor((length - offset) / 3);
  if (!Number.isFinite(length) || codonCount <= 0) return 0;
  let sum = 0;
  let informativeColumns = 0;
  for (let codonIndex = 0; codonIndex < codonCount; codonIndex++) {
    const counts = new Map<string, number>();
    const start = offset + codonIndex * 3;
    for (const sequence of sequences) {
      const aminoAcid = translateClientCodon(sequence.slice(start, start + 3), geneticCode);
      if (!['*', 'X', '-'].includes(aminoAcid)) {
        counts.set(aminoAcid, (counts.get(aminoAcid) ?? 0) + 1);
      }
    }
    const total = Array.from(counts.values()).reduce((sum, count) => sum + count, 0);
    if (total >= 2) {
      sum += Math.max(...counts.values()) / total;
      informativeColumns++;
    }
  }
  return informativeColumns > 0 ? sum / informativeColumns : 0;
}

function clientSynonymousFraction(sequences: string[], geneticCode: GeneticCode): number {
  const length = Math.min(...sequences.map((sequence) => sequence.length));
  const codonCount = Math.floor(length / 3);
  let substitutions = 0;
  let synonymous = 0;
  for (let codonIndex = 0; codonIndex < codonCount; codonIndex++) {
    const start = codonIndex * 3;
    const counts = new Map<string, number>();
    for (const sequence of sequences) {
      const codon = sequence.slice(start, start + 3).toUpperCase().replace(/U/g, 'T');
      const aminoAcid = translateClientCodon(codon, geneticCode);
      if (!['*', 'X', '-'].includes(aminoAcid)) {
        counts.set(codon, (counts.get(codon) ?? 0) + 1);
      }
    }
    const consensus = Array.from(counts.entries()).sort((a, b) => b[1] - a[1])[0]?.[0];
    if (!consensus) continue;
    const consensusAminoAcid = translateClientCodon(consensus, geneticCode);
    for (const [codon, count] of counts) {
      if (codon !== consensus) {
        substitutions += count;
        if (translateClientCodon(codon, geneticCode) === consensusAminoAcid) synonymous += count;
      }
    }
  }
  return substitutions >= 5 ? synonymous / substitutions : 0;
}

function clientCodonPositionVariability(sequences: string[]): [number, number, number] {
  const length = Math.min(...sequences.map((sequence) => sequence.length));
  const codonCount = Math.floor(length / 3);
  const sums = [0, 0, 0];
  const columns = [0, 0, 0];
  for (let codonIndex = 0; codonIndex < codonCount; codonIndex++) {
    for (let position = 0; position < 3; position++) {
      const counts = new Map<string, number>();
      const column = codonIndex * 3 + position;
      for (const sequence of sequences) {
        const state = resolvedNucleotideState(sequence[column]);
        if (state) counts.set(state, (counts.get(state) ?? 0) + 1);
      }
      const total = Array.from(counts.values()).reduce((sum, count) => sum + count, 0);
      if (total >= 2) {
        sums[position] += 1 - Math.max(...counts.values()) / total;
        columns[position]++;
      }
    }
  }
  return sums.map((sum, position) => columns[position] > 0 ? sum / columns[position] : 0) as [number, number, number];
}

function clientCodingEvidence(
  orientedSequences: string[],
  start: number,
  end: number,
  geneticCode: GeneticCode
): { score: number; conservation: number; contrast: number } {
  const segments = orientedSequences
    .filter((sequence) => sequence.length >= end)
    .map((sequence) => sequence.slice(start, end));
  if (segments.length < 2 || end <= start) return { score: 0, conservation: 0, contrast: 0 };
  const conservation = clientProteinProfileConservation(segments, 0, geneticCode);
  const reverseSegments = segments.map(reverseComplementClient);
  const alternative = Math.max(
    clientProteinProfileConservation(segments, 1, geneticCode),
    clientProteinProfileConservation(segments, 2, geneticCode),
    clientProteinProfileConservation(reverseSegments, 0, geneticCode),
    clientProteinProfileConservation(reverseSegments, 1, geneticCode),
    clientProteinProfileConservation(reverseSegments, 2, geneticCode)
  );
  const contrast = Math.max(0, conservation - alternative);
  const synonymous = clientSynonymousFraction(segments, geneticCode);
  const variability = clientCodonPositionVariability(segments);
  const periodicity = Math.max(0, variability[2] - (variability[0] + variability[1]) / 2);
  const startFraction = segments.filter((segment) => segment.startsWith('ATG')).length / segments.length;
  const boundaryStopFraction = orientedSequences.filter((sequence) =>
    translateClientCodon(sequence.slice(end, end + 3), geneticCode) === '*'
  ).length / orientedSequences.length;
  const clamp = (value: number) => Math.min(1, Math.max(0, value));
  const score =
    clamp((conservation - 0.75) / 0.20) * 25 +
    clamp(contrast / 0.08) * 30 +
    clamp((synonymous - 0.25) / 0.50) * 25 +
    clamp(periodicity / 0.10) * 10 +
    startFraction * 5 +
    boundaryStopFraction * 5;
  return { score: Math.min(100, Math.max(0, score)), conservation, contrast };
}

/**
 * Finds the aligned stop-free interval that maximizes sample support × length
 * across all six reading frames. A sample supports an interval only when it has
 * no stop codon there and at least 25% of its codons are resolved nucleotides.
 */
function findClientBestSharedOrfSegment(
  sequences: string[],
  geneticCode: GeneticCode,
  minSharedSupportPercent: number,
  minSegmentAa: number,
  minCodingScore: number
): ClientSharedOrfSegment | null {
  if (sequences.length === 0) return null;

  const stops = clientStopCodonsForCode(geneticCode);
  const minimumCodons = Math.max(1, Math.floor(minSegmentAa));
  const minSupport = Math.min(
    sequences.length,
    Math.max(1, Math.ceil(sequences.length * Math.min(100, Math.max(0, minSharedSupportPercent)) / 100))
  );
  const candidates: ClientSharedOrfSegment[] = [];

  for (const isReverse of [false, true]) {
    const oriented = isReverse ? sequences.map(reverseComplementClient) : sequences;
    const alignmentLength = Math.min(...oriented.map((sequence) => sequence.length));

    for (let offset = 0; offset < 3; offset++) {
      if (alignmentLength < offset + 3) continue;
      const codonCount = Math.floor((alignmentLength - offset) / 3);
      if (codonCount < minimumCodons) continue;

      const stopPrefixes: number[][] = [];
      const resolvedPrefixes: number[][] = [];
      const candidateIntervals = new Set<string>();

      for (const sequence of oriented) {
        const stopPrefix = new Array<number>(codonCount + 1).fill(0);
        const resolvedPrefix = new Array<number>(codonCount + 1).fill(0);
        let runStart = 0;

        for (let codonIndex = 0; codonIndex < codonCount; codonIndex++) {
          const start = offset + codonIndex * 3;
          const codon = sequence.slice(start, start + 3).toUpperCase().replace(/U/g, 'T');
          const isStop = stops.has(codon);
          const isResolved = !isStop && /^[ACGT]{3}$/.test(codon);
          stopPrefix[codonIndex + 1] = stopPrefix[codonIndex] + (isStop ? 1 : 0);
          resolvedPrefix[codonIndex + 1] = resolvedPrefix[codonIndex] + (isResolved ? 1 : 0);

          if (isStop) {
            if (codonIndex - runStart >= minimumCodons) {
              candidateIntervals.add(`${runStart}:${codonIndex}`);
            }
            runStart = codonIndex + 1;
          }
        }

        if (codonCount - runStart >= minimumCodons) {
          candidateIntervals.add(`${runStart}:${codonCount}`);
        }
        stopPrefixes.push(stopPrefix);
        resolvedPrefixes.push(resolvedPrefix);
      }

      for (const interval of candidateIntervals) {
        const [startCodon, endCodon] = interval.split(':').map(Number);
        const lengthCodons = endCodon - startCodon;
        const minimumResolved = Math.max(1, Math.ceil(lengthCodons * 0.25));
        let supportCount = 0;
        let informativeCodons = 0;

        for (let sampleIndex = 0; sampleIndex < oriented.length; sampleIndex++) {
          const stopCount =
            stopPrefixes[sampleIndex][endCodon] - stopPrefixes[sampleIndex][startCodon];
          const resolvedCount =
            resolvedPrefixes[sampleIndex][endCodon] - resolvedPrefixes[sampleIndex][startCodon];
          if (stopCount === 0 && resolvedCount >= minimumResolved) {
            supportCount++;
            informativeCodons += resolvedCount;
          }
        }

        if (supportCount < minSupport) continue;
        const candidate: ClientSharedOrfSegment = {
          frame: isReverse ? -(offset + 1) : offset + 1,
          isReverse,
          start: offset + startCodon * 3,
          end: offset + endCodon * 3,
          supportCount,
          lengthCodons,
          informativeCodons,
          codingScore: 0,
          aminoAcidConservation: 0,
          frameContrast: 0,
        };
        candidates.push(candidate);
        candidates.sort((left, right) =>
          (right.supportCount * right.lengthCodons - left.supportCount * left.lengthCodons) ||
          (right.supportCount - left.supportCount) ||
          (right.lengthCodons - left.lengthCodons) ||
          (right.informativeCodons - left.informativeCodons)
        );
        if (candidates.length > 64) candidates.length = 64;
      }
    }
  }

  for (const candidate of candidates) {
    const oriented = candidate.isReverse ? sequences.map(reverseComplementClient) : sequences;
    const evidence = clientCodingEvidence(oriented, candidate.start, candidate.end, geneticCode);
    candidate.codingScore = evidence.score;
    candidate.aminoAcidConservation = evidence.conservation;
    candidate.frameContrast = evidence.contrast;
  }
  return candidates.find((candidate) => candidate.codingScore >= minCodingScore) ??
    candidates.sort((left, right) =>
      (right.codingScore - left.codingScore) ||
      (right.supportCount * right.lengthCodons - left.supportCount * left.lengthCodons)
    )[0] ?? null;
}

function detectClientStopCodons(
  taxa: string[],
  sequences: string[],
  geneticCode: GeneticCode
): StopCodonPos[] {
  const stops = clientStopCodonsForCode(geneticCode);
  const result: StopCodonPos[] = [];
  taxa.forEach((taxon, taxonIndex) => {
    const sequence = sequences[taxonIndex] || '';
    const codonCount = Math.floor(sequence.length / 3);
    for (let codonIndex = 0; codonIndex < codonCount; codonIndex++) {
      const start = codonIndex * 3;
      const codon = sequence.slice(start, start + 3).toUpperCase();
      if (stops.has(codon)) {
        result.push({
          taxon,
          start,
          end: start + 3,
          codon,
          is_terminal: codonIndex + 1 === codonCount,
        });
      }
    }
  });
  return result;
}

function detectClientRawStopsInFrame(
  taxa: string[],
  sequences: string[],
  frame: number,
  geneticCode: GeneticCode,
  rawColMap: number[],
  selectedRegionRawMap: number[]
): StopCodonPos[] {
  if (frame === 0) return [];
  const isReverse = frame < 0;
  const frameOffset = Math.max(0, Math.abs(frame) - 1);
  const orientedSequences = isReverse ? sequences.map(reverseComplementClient) : sequences;
  const orientedMap = isReverse ? [...rawColMap].reverse() : rawColMap;
  const selectedOriginColumn = isReverse
    ? selectedRegionRawMap[selectedRegionRawMap.length - 1]
    : selectedRegionRawMap[0];
  const selectedOrigin = Math.max(0, orientedMap.indexOf(selectedOriginColumn));
  const offset = (selectedOrigin + frameOffset) % 3;
  const result: StopCodonPos[] = [];
  taxa.forEach((taxon, taxonIndex) => {
    const sequence = orientedSequences[taxonIndex] || '';
    const codonCount = Math.floor((sequence.length - offset) / 3);
    for (let codonIndex = 0; codonIndex < codonCount; codonIndex++) {
      const localStart = offset + codonIndex * 3;
      const codon = sequence.slice(localStart, localStart + 3).toUpperCase();
      if (translateClientCodon(codon, geneticCode) !== '*') continue;
      const mapped = orientedMap.slice(localStart, localStart + 3);
      if (mapped.length !== 3) continue;
      const rawStart = Math.min(...mapped);
      const rawEnd = Math.max(...mapped) + 1;
      if (rawEnd - rawStart === 3) {
        result.push({
          taxon,
          start: rawStart,
          end: rawEnd,
          codon,
          is_terminal: codonIndex + 1 === codonCount,
        });
      }
    }
  });
  return result;
}

interface ClientReferenceMatch {
  start: number;
  end: number;
  identity: number;
  coverage: number;
  isReverse: boolean;
}

function normalizeClientReference(sequence: string): string {
  return sequence.toUpperCase().replace(/U/g, 'T').replace(/[^ACGT]/g, '');
}

function clientConsensusWithMap(sequences: string[]): { consensus: string; columns: number[] } {
  const length = Math.min(...sequences.map((sequence) => sequence.length));
  const consensus: string[] = [];
  const columns: number[] = [];
  for (let column = 0; column < length; column++) {
    const counts = new Map<string, number>();
    for (const sequence of sequences) {
      const state = resolvedNucleotideState(sequence[column]);
      if (state) counts.set(state, (counts.get(state) ?? 0) + 1);
    }
    const selected = Array.from(counts.entries()).sort((left, right) => right[1] - left[1])[0];
    if (selected) {
      consensus.push(selected[0]);
      columns.push(column);
    }
  }
  return { consensus: consensus.join(''), columns };
}

function clientLocalReferenceMatch(reference: string, target: string):
  | { start: number; end: number; identity: number; coverage: number; score: number }
  | null {
  if (reference.length < 12 || target.length < 12) return null;
  const width = target.length + 1;
  let previousScore = new Int32Array(width);
  let currentScore = new Int32Array(width);
  let previousStart = new Int32Array(width);
  let currentStart = new Int32Array(width);
  let previousMatches = new Int32Array(width);
  let currentMatches = new Int32Array(width);
  let previousReferenceBases = new Int32Array(width);
  let currentReferenceBases = new Int32Array(width);
  let previousTargetBases = new Int32Array(width);
  let currentTargetBases = new Int32Array(width);
  let bestScore = 0;
  let bestStart = 0;
  let bestEnd = 0;
  let bestMatches = 0;
  let bestReferenceBases = 0;
  let bestTargetBases = 0;

  for (let referenceIndex = 1; referenceIndex <= reference.length; referenceIndex++) {
    currentScore.fill(0);
    currentStart.fill(0);
    currentMatches.fill(0);
    currentReferenceBases.fill(0);
    currentTargetBases.fill(0);
    for (let targetIndex = 1; targetIndex <= target.length; targetIndex++) {
      const matched = reference[referenceIndex - 1] === target[targetIndex - 1];
      const diagonalScore = previousScore[targetIndex - 1] + (matched ? 3 : -2);
      const upScore = previousScore[targetIndex] - 3;
      const leftScore = currentScore[targetIndex - 1] - 3;
      const score = Math.max(0, diagonalScore, upScore, leftScore);
      if (score === 0) continue;

      if (score === diagonalScore) {
        currentStart[targetIndex] = previousScore[targetIndex - 1] > 0
          ? previousStart[targetIndex - 1]
          : targetIndex - 1;
        currentMatches[targetIndex] = previousMatches[targetIndex - 1] + (matched ? 1 : 0);
        currentReferenceBases[targetIndex] = previousReferenceBases[targetIndex - 1] + 1;
        currentTargetBases[targetIndex] = previousTargetBases[targetIndex - 1] + 1;
      } else if (score === upScore) {
        currentStart[targetIndex] = previousStart[targetIndex];
        currentMatches[targetIndex] = previousMatches[targetIndex];
        currentReferenceBases[targetIndex] = previousReferenceBases[targetIndex] + 1;
        currentTargetBases[targetIndex] = previousTargetBases[targetIndex];
      } else {
        currentStart[targetIndex] = currentStart[targetIndex - 1];
        currentMatches[targetIndex] = currentMatches[targetIndex - 1];
        currentReferenceBases[targetIndex] = currentReferenceBases[targetIndex - 1];
        currentTargetBases[targetIndex] = currentTargetBases[targetIndex - 1] + 1;
      }
      currentScore[targetIndex] = score;
      if (score > bestScore) {
        bestScore = score;
        bestStart = currentStart[targetIndex];
        bestEnd = targetIndex;
        bestMatches = currentMatches[targetIndex];
        bestReferenceBases = currentReferenceBases[targetIndex];
        bestTargetBases = currentTargetBases[targetIndex];
      }
    }
    [previousScore, currentScore] = [currentScore, previousScore];
    [previousStart, currentStart] = [currentStart, previousStart];
    [previousMatches, currentMatches] = [currentMatches, previousMatches];
    [previousReferenceBases, currentReferenceBases] = [currentReferenceBases, previousReferenceBases];
    [previousTargetBases, currentTargetBases] = [currentTargetBases, previousTargetBases];
  }

  if (bestScore <= 0 || bestEnd <= bestStart || bestReferenceBases === 0) return null;
  return {
    start: bestStart,
    end: bestEnd,
    identity: bestMatches / Math.max(1, bestReferenceBases, bestTargetBases) * 100,
    coverage: bestReferenceBases / reference.length * 100,
    score: bestScore,
  };
}

function matchClientReferenceToAlignment(
  sequences: string[],
  referenceSequence: string
): ClientReferenceMatch | null {
  const { consensus, columns } = clientConsensusWithMap(sequences);
  const forward = normalizeClientReference(referenceSequence);
  const reverse = normalizeClientReference(reverseComplementClient(referenceSequence));
  const candidates = [[forward, false], [reverse, true]] as const;
  const matches = candidates
    .map(([reference, isReverse]) => {
      const match = clientLocalReferenceMatch(reference, consensus);
      return match ? { ...match, isReverse } : null;
    })
    .filter((match): match is NonNullable<typeof match> => match !== null)
    .filter((match) => match.identity >= 70 && match.coverage >= 70)
    .sort((left, right) => right.score - left.score);
  const selected = matches[0];
  if (!selected) return null;
  const start = columns[selected.start];
  const endColumn = columns[selected.end - 1];
  if (start === undefined || endColumn === undefined || endColumn < start) return null;
  return {
    start,
    end: endColumn + 1,
    identity: selected.identity,
    coverage: selected.coverage,
    isReverse: selected.isReverse,
  };
}

interface ClientSiteStatistics {
  variableCount: number;
  variablePercent: number;
  pisCount: number;
  pisPercent: number;
  pisMask: boolean[];
}

function resolvedNucleotideState(base: string | undefined): string | null {
  const upper = base?.toUpperCase();
  if (upper === 'U') return 'T';
  return upper === 'A' || upper === 'C' || upper === 'G' || upper === 'T' ? upper : null;
}

function convertClientAmbiguities(
  sequences: string[],
  strategy: TrimmingRecipe['ambiguity_strategy']
): string[] {
  if (strategy === 'keep' || sequences.length === 0) return sequences;

  const ambiguity = new Set(['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V']);
  if (strategy === 'converttogap') {
    return sequences.map((sequence) =>
      sequence
        .split('')
        .map((base) => (ambiguity.has(base.toUpperCase()) ? '-' : base))
        .join('')
    );
  }

  if (strategy === 'fixedstandard') {
    const fixed: Record<string, string> = {
      R: 'A', Y: 'T', S: 'G', W: 'A', K: 'T', M: 'A',
      B: 'T', D: 'T', H: 'T', V: 'A',
    };
    return sequences.map((sequence) =>
      sequence
        .split('')
        .map((base) => fixed[base.toUpperCase()] ?? base)
        .join('')
    );
  }

  const length = sequences[0]?.length || 0;
  const majority = new Array<string>(length).fill('A');
  for (let column = 0; column < length; column++) {
    const counts: Record<string, number> = { A: 0, C: 0, G: 0, T: 0 };
    for (const sequence of sequences) {
      const state = resolvedNucleotideState(sequence[column]);
      if (state) counts[state]++;
    }
    majority[column] = Object.entries(counts).reduce(
      (best, entry) => (entry[1] > best[1] ? entry : best),
      ['A', -1] as [string, number]
    )[0];
  }
  return sequences.map((sequence) =>
    sequence
      .split('')
      .map((base, column) =>
        ambiguity.has(base.toUpperCase()) ? majority[column] : base
      )
      .join('')
  );
}

function calculateClientSiteStatistics(sequences: string[]): ClientSiteStatistics {
  const length = sequences[0]?.length || 0;
  let variableCount = 0;
  let pisCount = 0;
  const pisMask = new Array<boolean>(length).fill(false);
  const numSeqs = sequences.length;

  for (let column = 0; column < length; column++) {
    let countA = 0, countC = 0, countG = 0, countT = 0;
    for (let i = 0; i < numSeqs; i++) {
      const char = sequences[i][column];
      if (char === 'A' || char === 'a') countA++;
      else if (char === 'C' || char === 'c') countC++;
      else if (char === 'G' || char === 'g') countG++;
      else if (char === 'T' || char === 't' || char === 'U' || char === 'u') countT++;
    }
    
    let diffStates = 0;
    let repeatedStates = 0;
    
    if (countA > 0) { diffStates++; if (countA >= 2) repeatedStates++; }
    if (countC > 0) { diffStates++; if (countC >= 2) repeatedStates++; }
    if (countG > 0) { diffStates++; if (countG >= 2) repeatedStates++; }
    if (countT > 0) { diffStates++; if (countT >= 2) repeatedStates++; }
    
    if (diffStates >= 2) variableCount++;
    if (repeatedStates >= 2) {
      pisCount++;
      pisMask[column] = true;
    }
  }

  return {
    variableCount,
    variablePercent: length > 0 ? (variableCount / length) * 100 : 0,
    pisCount,
    pisPercent: length > 0 ? (pisCount / length) * 100 : 0,
    pisMask,
  };
}

function calculateClientGapPercent(sequences: string[]): number {
  let total = 0;
  let gaps = 0;
  for (const sequence of sequences) {
    for (const base of sequence) {
      total++;
      if (base === '-' || base === '?' || base.toUpperCase() === 'N') gaps++;
    }
  }
  return total > 0 ? (gaps / total) * 100 : 0;
}

function calculateClientMajorityConsensus(sequences: string[]): string {
  const length = sequences[0]?.length || 0;
  const consensus: string[] = [];
  for (let column = 0; column < length; column++) {
    const counts: Record<string, number> = {};
    for (const sequence of sequences) {
      const state = resolvedNucleotideState(sequence[column]);
      if (state) counts[state] = (counts[state] || 0) + 1;
    }
    let majority = '-';
    let majorityCount = 0;
    for (const [state, count] of Object.entries(counts)) {
      if (count > majorityCount) {
        majority = state;
        majorityCount = count;
      }
    }
    consensus.push(majority);
  }
  return consensus.join('');
}

export function parseFastaText(text: string, fileName: string, filePath: string): Alignment {
  const lines = text.split(/\r?\n/);
  const taxa: string[] = [];
  const sequences: string[] = [];
  let currentHeader = '';
  let currentSeq = '';

  for (const rawLine of lines) {
    const line = rawLine.trim();
    if (!line) continue;

    if (line.startsWith('>') || line.startsWith(';')) {
      if (currentHeader) {
        taxa.push(currentHeader);
        sequences.push(currentSeq);
        currentSeq = '';
      }
      currentHeader = line.slice(1).trim();
    } else {
      currentSeq += line;
    }
  }

  if (currentHeader) {
    taxa.push(currentHeader);
    sequences.push(currentSeq);
  }

  const id = fileName.replace(/\.[^/.]+$/, '');
  return {
    id,
    file_name: fileName,
    file_path: filePath,
    format: 'fasta',
    taxa,
    sequences,
    length: sequences[0]?.length || 0,
    num_taxa: taxa.length,
  };
}

export function parsePhylipText(text: string, fileName: string, filePath: string): Alignment {
  const lines = text
    .split(/\r?\n/)
    .map((l) => l.trim())
    .filter((l) => l.length > 0);

  if (lines.length === 0) throw new Error('Empty PHYLIP file');

  const headerParts = lines[0].split(/\s+/);
  const expectedTaxa = parseInt(headerParts[0], 10);
  const expectedLength = parseInt(headerParts[1], 10);

  const taxa: string[] = [];
  const sequences: string[] = [];
  const remaining = lines.slice(1);

  if (remaining.length === expectedTaxa) {
    for (const line of remaining) {
      const tokens = line.split(/\s+/);
      if (tokens.length >= 2) {
        taxa.push(tokens[0]);
        sequences.push(tokens.slice(1).join(''));
      } else if (line.length > 10) {
        taxa.push(line.slice(0, 10).trim());
        sequences.push(line.slice(10).replace(/\s+/g, ''));
      }
    }
  } else {
    // Interleaved
    let curIdx = 0;
    for (const line of remaining) {
      const tokens = line.split(/\s+/);
      if (tokens.length === 0) continue;

      if (taxa.length < expectedTaxa) {
        if (tokens.length >= 2) {
          taxa.push(tokens[0]);
          sequences.push(tokens.slice(1).join(''));
        } else if (line.length > 10) {
          taxa.push(line.slice(0, 10).trim());
          sequences.push(line.slice(10).replace(/\s+/g, ''));
        }
      } else {
        if (curIdx >= expectedTaxa) curIdx = 0;
        const chunk = tokens.length >= 2 && tokens[0] === taxa[curIdx]
          ? tokens.slice(1).join('')
          : tokens.join('');
        sequences[curIdx] += chunk;
        curIdx++;
      }
    }
  }

  const id = fileName.replace(/\.[^/.]+$/, '');
  return {
    id,
    file_name: fileName,
    file_path: filePath,
    format: 'phylip',
    taxa,
    sequences,
    length: sequences[0]?.length || expectedLength || 0,
    num_taxa: taxa.length,
  };
}

export function parseNexusText(text: string, fileName: string, filePath: string): Alignment {
  const lines = text.split(/\r?\n/);
  let inMatrix = false;
  const taxa: string[] = [];
  const sequences: string[] = [];
  const map = new Map<string, number>();

  for (const raw of lines) {
    const trimmed = raw.trim();
    if (!trimmed || trimmed.startsWith('[')) continue;

    const upper = trimmed.toUpperCase();
    if (upper.startsWith('MATRIX')) {
      inMatrix = true;
      continue;
    }

    if (inMatrix) {
      if (trimmed === ';' || upper.startsWith('END;') || upper === 'END') {
        inMatrix = false;
        break;
      }
      const tokens = trimmed.split(/\s+/);
      if (tokens.length >= 2) {
        const name = tokens[0].replace(/['"]/g, '');
        const seq = tokens.slice(1).join('');
        if (map.has(name)) {
          const idx = map.get(name)!;
          sequences[idx] += seq;
        } else {
          map.set(name, taxa.length);
          taxa.push(name);
          sequences.push(seq);
        }
      }
    }
  }

  const id = fileName.replace(/\.[^/.]+$/, '');
  return {
    id,
    file_name: fileName,
    file_path: filePath,
    format: 'nexus',
    taxa,
    sequences,
    length: sequences[0]?.length || 0,
    num_taxa: taxa.length,
  };
}

export function parseAlignmentText(text: string, fileName: string, filePath: string): Alignment {
  const ext = fileName.split('.').pop()?.toLowerCase() || '';
  if (['fa', 'fasta', 'fna', 'faa'].includes(ext)) {
    return parseFastaText(text, fileName, filePath);
  } else if (['phy', 'phylip'].includes(ext)) {
    return parsePhylipText(text, fileName, filePath);
  } else if (['nex', 'nexus'].includes(ext)) {
    return parseNexusText(text, fileName, filePath);
  } else {
    // Auto-detect by content
    if (text.trim().startsWith('>') || text.trim().startsWith(';')) {
      return parseFastaText(text, fileName, filePath);
    }
    if (text.toUpperCase().includes('#NEXUS')) {
      return parseNexusText(text, fileName, filePath);
    }
    return parsePhylipText(text, fileName, filePath);
  }
}

/** Ordinary alignment processing used by Catalog QC and general export. */
export function recipeWithoutOrfAnalysis(recipe: TrimmingRecipe): TrimmingRecipe {
  return {
    ...recipe,
    enable_orf: false,
    orf_use_references: false,
    fail_if_no_orf: false,
  };
}

export function computeAlignmentSummary(
  rawAlign: Alignment,
  recipe: TrimmingRecipe,
  totalDatasetTaxa: number = 0
): AlignmentSummary {
  const rawNumTaxa = rawAlign.taxa.length;
  const rawLength = rawAlign.length;

  let rawTotalChars = rawNumTaxa * rawLength;
  let rawGapCount = 0;
  for (const seq of rawAlign.sequences) {
    const gaps = (seq.match(/[-?Nn]/g) || []).length;
    rawGapCount += gaps;
  }
  const rawGapPercent = rawTotalChars > 0 ? (rawGapCount / rawTotalChars) * 100 : 0;

  // Catalog alignment QC and ORF analysis are independent branches. Ordinary
  // metrics and pass/fail always come from the ORF-disabled branch.
  const trimmedView = executeClientTrimming(
    rawAlign,
    recipeWithoutOrfAnalysis(recipe),
    totalDatasetTaxa
  );
  const orfView = recipe.enable_orf
    ? executeClientTrimming(rawAlign, recipe, totalDatasetTaxa)
    : trimmedView;
  const orfDiff = orfView.diff;
  const align = trimmedView.trimmed_alignment;

  const numTaxa = align.taxa.length;
  const length = align.length;
  let totalChars = 0;
  let gapCount = 0;
  let gcCount = 0;
  let nonGapCount = 0;
  const retainedTaxonBasepairs: Record<string, number> = {};

  for (let sequenceIndex = 0; sequenceIndex < align.sequences.length; sequenceIndex++) {
    const seq = align.sequences[sequenceIndex];
    let sampleBasepairs = 0;
    for (let i = 0; i < seq.length; i++) {
      totalChars++;
      const c = seq[i].toUpperCase();
      if (c === '-' || c === '?' || c === 'N') {
        gapCount++;
      } else {
        sampleBasepairs++;
        nonGapCount++;
        if (c === 'G' || c === 'C') gcCount++;
      }
    }
    const taxon = align.taxa[sequenceIndex];
    if (taxon) {
      retainedTaxonBasepairs[taxon] =
        (retainedTaxonBasepairs[taxon] ?? 0) + sampleBasepairs;
    }
  }

  const gapPercent = totalChars > 0 ? (gapCount / totalChars) * 100 : 0;
  const gcPercent = nonGapCount > 0 ? (gcCount / nonGapCount) * 100 : 0;

  const siteStats = calculateClientSiteStatistics(align.sequences);
  const variableCount = siteStats.variableCount;
  const variablePercent = siteStats.variablePercent;
  const pisCount = siteStats.pisCount;
  const pisPercent = siteStats.pisPercent;

  // Majority consensus
  const consensusChars: string[] = [];
  for (let col = 0; col < length; col++) {
    const counts: Record<string, number> = {};
    for (const seq of align.sequences) {
      const state = resolvedNucleotideState(seq[col]);
      if (state) {
        counts[state] = (counts[state] || 0) + 1;
      }
    }
    let maxChar = '-';
    let maxCnt = 0;
    for (const [ch, cnt] of Object.entries(counts)) {
      if (cnt > maxCnt) {
        maxCnt = cnt;
        maxChar = ch;
      }
    }
    consensusChars.push(maxChar);
  }
  const consensus = consensusChars.join('');

  // Mean divergence to consensus
  let sumDist = 0;
  for (const seq of align.sequences) {
    let mismatches = 0;
    let overlap = 0;
    for (let col = 0; col < length; col++) {
      const c = seq[col]?.toUpperCase();
      const cons = consensus[col];
      if (c && cons && c !== '-' && c !== '?' && c !== 'N' && cons !== '-') {
        overlap++;
        if (c !== cons) mismatches++;
      }
    }
    sumDist += overlap > 0 ? mismatches / overlap : 0;
  }
  const meanDivergence = numTaxa > 0 ? sumDist / numTaxa : 0;

  // Gating against post-trimmed metrics
  const failReasons: string[] = [...(trimmedView.diff.fail_reasons || [])];

  if (numTaxa === 0) {
    const r0 = '0 surviving taxa (all samples pruned)';
    if (!failReasons.includes(r0)) failReasons.push(r0);
  }

  if (recipe.assess_alignment) {
    if (recipe.min_taxa > 0 && numTaxa < recipe.min_taxa) {
      const r = `Taxa (${numTaxa} < min ${recipe.min_taxa})`;
      if (!failReasons.includes(r)) failReasons.push(r);
    }
    if (totalDatasetTaxa > 0 && recipe.min_taxa_occupancy_percent > 0) {
      const occPct = (numTaxa / totalDatasetTaxa) * 100;
      if (occPct < recipe.min_taxa_occupancy_percent) {
        const r = `Taxon occupancy (${occPct.toFixed(1)}% < min ${recipe.min_taxa_occupancy_percent}%, ${numTaxa}/${totalDatasetTaxa} taxa)`;
        if (!failReasons.includes(r)) failReasons.push(r);
      }
    }
    if (recipe.min_length > 0 && length < recipe.min_length) {
      const r = `Length (${length} bp < min ${recipe.min_length} bp)`;
      if (!failReasons.includes(r)) failReasons.push(r);
    }
    if (recipe.max_gap_percent > 0 && gapPercent > recipe.max_gap_percent) {
      const r = `Gap % (${gapPercent.toFixed(1)}% > max ${recipe.max_gap_percent}%)`;
      if (!failReasons.includes(r)) failReasons.push(r);
    }
    if ((recipe.min_variable_count ?? 0) > 0 && variableCount < recipe.min_variable_count) {
      const r = `Variable sites (${variableCount} < min ${recipe.min_variable_count})`;
      if (!failReasons.includes(r)) failReasons.push(r);
    }
    if (
      (recipe.min_variable_percent ?? 0) > 0 &&
      variablePercent < recipe.min_variable_percent
    ) {
      const r = `Variable-site proportion (${variablePercent.toFixed(1)}% < min ${recipe.min_variable_percent.toFixed(1)}%)`;
      if (!failReasons.includes(r)) failReasons.push(r);
    }
    if ((recipe.min_pis_count ?? 0) > 0 && pisCount < recipe.min_pis_count) {
      const r = `Parsimony-informative sites (${pisCount} < min ${recipe.min_pis_count})`;
      if (!failReasons.includes(r)) failReasons.push(r);
    }
    if ((recipe.min_pis_percent ?? 0) > 0 && pisPercent < recipe.min_pis_percent) {
      const r = `Parsimony-informative proportion (${pisPercent.toFixed(1)}% < min ${recipe.min_pis_percent.toFixed(1)}%)`;
      if (!failReasons.includes(r)) failReasons.push(r);
    }
  }

  const pass = trimmedView.diff.pass && failReasons.length === 0 && numTaxa > 0;

  return {
    id: rawAlign.id,
    file_name: rawAlign.file_name,
    file_path: rawAlign.file_path,
    format: rawAlign.format,
    num_taxa: numTaxa,
    length,
    total_basepairs: totalChars - gapCount,
    gap_count: gapCount,
    gap_percent: Number(gapPercent.toFixed(1)),
    variable_count: variableCount,
    variable_percent: Number(variablePercent.toFixed(1)),
    pis_count: pisCount,
    pis_percent: Number(pisPercent.toFixed(1)),
    mean_divergence: Number(meanDivergence.toFixed(3)),
    gc_percent: Number(gcPercent.toFixed(1)),
    pass,
    fail_reasons: failReasons,
    orf_valid: orfDiff.found_valid_orf ?? true,
    orf_evaluated: orfDiff.orf_evaluated ?? false,
    orf_candidate_found: orfDiff.orf_candidate_found ?? false,
    orf_frame: orfDiff.orf_frame,
    orf_start: orfDiff.orf_start,
    orf_end: orfDiff.orf_end,
    orf_support_count: orfDiff.orf_support_count ?? 0,
    orf_support_percent: orfDiff.orf_support_percent ?? 0,
    orf_retained_samples: orfDiff.orf_retained_samples ?? 0,
    orf_candidate_length_aa: orfDiff.orf_candidate_length_aa ?? 0,
    orf_coding_score: orfDiff.orf_coding_score ?? 0,
    orf_amino_acid_conservation: orfDiff.orf_amino_acid_conservation ?? 0,
    orf_frame_contrast: orfDiff.orf_frame_contrast ?? 0,
    orf_reference_evaluated: orfDiff.orf_reference_evaluated ?? false,
    orf_reference_matched: orfDiff.orf_reference_matched ?? false,
    orf_reference_identity: orfDiff.orf_reference_identity ?? 0,
    orf_reference_coverage: orfDiff.orf_reference_coverage ?? 0,
    orf_intron_length: orfDiff.orf_intron_length ?? 0,
    raw_num_taxa: rawNumTaxa,
    raw_length: rawLength,
    raw_gap_percent: Number(rawGapPercent.toFixed(1)),
    retained_taxa: [...align.taxa],
    retained_taxon_basepairs: retainedTaxonBasepairs,
    orf_retained_taxa: [...orfView.trimmed_alignment.taxa],
    orf_retained_taxon_basepairs: Object.fromEntries(
      orfView.trimmed_alignment.taxa.map((taxon, index) => {
        const sequence = orfView.trimmed_alignment.sequences[index] ?? '';
        const basepairs = Array.from(sequence).filter(
          (state) => !['-', '?', 'N'].includes(state.toUpperCase())
        ).length;
        return [taxon, basepairs];
      })
    ),
  };
}

type OrfAuditMetadata = {
  orf_candidate_found?: boolean;
  orf_coding_score?: number;
  orf_retained_samples?: number;
  orf_frame?: number;
  orf_reference_evaluated?: boolean;
  orf_reference_matched?: boolean;
};

function candidateOrfFailureReason(
  metadata: OrfAuditMetadata,
  recipe: TrimmingRecipe
): string {
  if (
    recipe.orf_search_mode === 'referenceguided' &&
    metadata.orf_reference_evaluated &&
    !metadata.orf_reference_matched
  ) {
    return 'Reference-guided exon failed: no matching reference sequence could be anchored for this locus';
  }
  if (!metadata.orf_candidate_found) {
    if (
      recipe.orf_search_mode === 'continuouscds' ||
      recipe.orf_search_mode === 'referenceguided'
    ) {
      return 'Candidate ORF failed: no usable continuous reading frame was found';
    }
    return `Candidate ORF failed: no segment met ${recipe.orf_min_shared_support_percent.toFixed(0)}% support and ${recipe.orf_min_segment_aa} aa minimums`;
  }
  if ((metadata.orf_frame ?? 0) < 0 && !recipe.auto_flip_reverse) {
    return 'Candidate ORF failed: the best candidate is reverse-strand but Auto-Flip Reverse Strand is disabled';
  }
  if (
    (recipe.orf_search_mode === 'bestsharedsegment' ||
      recipe.orf_search_mode === 'referencecandidateorf') &&
    (metadata.orf_coding_score ?? 0) < recipe.orf_min_coding_score
  ) {
    return `Candidate ORF failed: coding evidence ${(metadata.orf_coding_score ?? 0).toFixed(1)} < minimum ${recipe.orf_min_coding_score.toFixed(1)}`;
  }
  if ((metadata.orf_retained_samples ?? 0) === 0) {
    return 'Candidate ORF failed: no samples remain after stop-codon screening';
  }
  return 'Candidate ORF failed: premature stop codons remain in retained samples';
}

/** Adds a derived exclusion list for taxa below a dataset-wide locus occupancy threshold. */
export function recipeWithDatasetSampleFilter(
  recipe: TrimmingRecipe,
  alignments: Alignment[]
): TrimmingRecipe {
  const threshold = recipe.min_sample_locus_occupancy_percent ?? 0;
  if (threshold <= 0 || alignments.length === 0) {
    return { ...recipe, excluded_taxa: [] };
  }

  const presenceCounts = new Map<string, number>();
  for (const alignment of alignments) {
    for (const taxon of new Set(alignment.taxa)) {
      presenceCounts.set(taxon, (presenceCounts.get(taxon) ?? 0) + 1);
    }
  }

  const excludedTaxa = Array.from(presenceCounts.entries())
    .filter(([, count]) => (count / alignments.length) * 100 < threshold)
    .map(([taxon]) => taxon)
    .sort();
  return { ...recipe, excluded_taxa: excludedTaxa };
}

export async function buildScanResponseFromAlignments(
  alignments: Alignment[],
  recipe: TrimmingRecipe,
  totalDatasetTaxa: number = 0,
  onProgress?: (percent: number) => void
): Promise<ScanResponse> {
  const uniqueTaxa = new Set(alignments.flatMap((alignment) => alignment.taxa)).size;
  const assessmentTaxa = totalDatasetTaxa > 0 ? totalDatasetTaxa : uniqueTaxa;
  const runtimeRecipe = recipeWithDatasetSampleFilter(recipe, alignments);
  
  const summaries: AlignmentSummary[] = [];
  for (let i = 0; i < alignments.length; i++) {
    summaries.push(computeAlignmentSummary(alignments[i], runtimeRecipe, assessmentTaxa));
    if (i % 20 === 0) {
      if (onProgress) onProgress(90 + (i / alignments.length) * 9.9);
      await new Promise(r => setTimeout(r, 0));
    }
  }
  const totalLoci = summaries.length;

  // Taxa occupancy map
  const taxaMap = new Map<string, { count: number; totalBp: number; sumGap: number }>();
  for (const align of alignments) {
    for (let i = 0; i < align.taxa.length; i++) {
      const taxon = align.taxa[i];
      const seq = align.sequences[i] || '';
      const gaps = (seq.match(/[-?Nn]/g) || []).length;
      const bp = seq.length - gaps;
      const gapPct = seq.length > 0 ? (gaps / seq.length) * 100 : 100;

      const entry = taxaMap.get(taxon) || { count: 0, totalBp: 0, sumGap: 0 };
      entry.count++;
      entry.totalBp += bp;
      entry.sumGap += gapPct;
      taxaMap.set(taxon, entry);
    }
  }

  const occupancy: TaxonOccupancy[] = Array.from(taxaMap.entries()).map(([taxon_name, stats]) => ({
    taxon_name,
    present_loci_count: stats.count,
    present_loci_percent: Number(((stats.count / Math.max(1, totalLoci)) * 100).toFixed(1)),
    mean_gap_percent: Number((stats.sumGap / Math.max(1, stats.count)).toFixed(1)),
    total_bp: stats.totalBp,
  }));

  occupancy.sort((a, b) => b.present_loci_count - a.present_loci_count);

  const overview = buildDatasetOverviewFromSummaries(summaries, assessmentTaxa);

  return { summaries, overview, occupancy };
}

export function buildDatasetOverviewFromSummaries(
  summaries: AlignmentSummary[],
  totalUniqueTaxa: number
): DatasetOverview {
  const totalLoci = summaries.length;
  const passedCount = summaries.filter((summary) => summary.pass).length;
  const sumTaxa = summaries.reduce((acc, summary) => acc + summary.num_taxa, 0);
  const sumLength = summaries.reduce((acc, summary) => acc + summary.length, 0);
  const sumGap = summaries.reduce((acc, summary) => acc + summary.gap_percent, 0);
  const sumPis = summaries.reduce((acc, summary) => acc + summary.pis_count, 0);
  const totalBasepairs = summaries.reduce((acc, summary) => acc + summary.total_basepairs, 0);

  return {
    total_alignments: totalLoci,
    passed_alignments: passedCount,
    discarded_alignments: totalLoci - passedCount,
    total_unique_taxa: totalUniqueTaxa,
    mean_taxa: totalLoci > 0 ? Number((sumTaxa / totalLoci).toFixed(1)) : 0,
    mean_length: totalLoci > 0 ? Number((sumLength / totalLoci).toFixed(0)) : 0,
    mean_gap_percent: totalLoci > 0 ? Number((sumGap / totalLoci).toFixed(1)) : 0,
    mean_pis: totalLoci > 0 ? Number((sumPis / totalLoci).toFixed(0)) : 0,
    total_matrix_basepairs: totalBasepairs,
  };
}

/**
 * Re-evaluates quality pass/fail criteria across all summaries instantaneously in memory (0ms).
 */
export function reevaluateSummaries(
  summaries: AlignmentSummary[],
  recipe: TrimmingRecipe,
  totalDatasetTaxa: number
): { summaries: AlignmentSummary[]; overview: DatasetOverview } {
  if (summaries.length === 0) {
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
    };
  }

  const updated: AlignmentSummary[] = summaries.map((s) => {
    const failReasons: string[] = [];

    if (s.num_taxa === 0) {
      failReasons.push('0 surviving taxa (all samples pruned)');
    }

    if (recipe.assess_alignment) {
      if (recipe.min_taxa > 0 && s.num_taxa < recipe.min_taxa) {
        failReasons.push(`Taxa (${s.num_taxa} < min ${recipe.min_taxa})`);
      }
      if (totalDatasetTaxa > 0 && recipe.min_taxa_occupancy_percent > 0) {
        const occPct = (s.num_taxa / totalDatasetTaxa) * 100;
        if (occPct < recipe.min_taxa_occupancy_percent) {
          failReasons.push(
            `Taxon occupancy (${occPct.toFixed(1)}% < min ${recipe.min_taxa_occupancy_percent}%, ${s.num_taxa}/${totalDatasetTaxa} taxa)`
          );
        }
      }
      if (recipe.min_length > 0 && s.length < recipe.min_length) {
        failReasons.push(`Length (${s.length} bp < min ${recipe.min_length} bp)`);
      }
      if (recipe.max_gap_percent > 0 && s.gap_percent > recipe.max_gap_percent) {
        failReasons.push(`Gap % (${s.gap_percent.toFixed(1)}% > max ${recipe.max_gap_percent}%)`);
      }
      if (
        (recipe.min_variable_count ?? 0) > 0 &&
        (s.variable_count ?? 0) < recipe.min_variable_count
      ) {
        failReasons.push(
          `Variable sites (${s.variable_count ?? 0} < min ${recipe.min_variable_count})`
        );
      }
      if (
        (recipe.min_variable_percent ?? 0) > 0 &&
        (s.variable_percent ?? 0) < recipe.min_variable_percent
      ) {
        failReasons.push(
          `Variable-site proportion (${(s.variable_percent ?? 0).toFixed(1)}% < min ${recipe.min_variable_percent.toFixed(1)}%)`
        );
      }
      if ((recipe.min_pis_count ?? 0) > 0 && s.pis_count < recipe.min_pis_count) {
        failReasons.push(
          `Parsimony-informative sites (${s.pis_count} < min ${recipe.min_pis_count})`
        );
      }
      if ((recipe.min_pis_percent ?? 0) > 0 && s.pis_percent < recipe.min_pis_percent) {
        failReasons.push(
          `Parsimony-informative proportion (${s.pis_percent.toFixed(1)}% < min ${recipe.min_pis_percent.toFixed(1)}%)`
        );
      }
    }

    return {
      ...s,
      pass: failReasons.length === 0,
      fail_reasons: failReasons,
    };
  });

  const passedCount = updated.filter((s) => s.pass).length;
  const discardedCount = updated.length - passedCount;
  const totalLoci = updated.length;

  const sumTaxa = updated.reduce((acc, s) => acc + s.num_taxa, 0);
  const sumLen = updated.reduce((acc, s) => acc + s.length, 0);
  const sumGap = updated.reduce((acc, s) => acc + s.gap_percent, 0);
  const sumPis = updated.reduce((acc, s) => acc + s.pis_count, 0);
  const totalBp = updated.reduce((acc, s) => acc + s.total_basepairs, 0);

  const overview: DatasetOverview = {
    total_alignments: totalLoci,
    passed_alignments: passedCount,
    discarded_alignments: discardedCount,
    total_unique_taxa: totalDatasetTaxa,
    mean_taxa: totalLoci > 0 ? Number((sumTaxa / totalLoci).toFixed(1)) : 0,
    mean_length: totalLoci > 0 ? Number((sumLen / totalLoci).toFixed(0)) : 0,
    mean_gap_percent: totalLoci > 0 ? Number((sumGap / totalLoci).toFixed(1)) : 0,
    mean_pis: totalLoci > 0 ? Number((sumPis / totalLoci).toFixed(0)) : 0,
    total_matrix_basepairs: totalBp,
  };

  return { summaries: updated, overview };
}

/** Re-gates an already processed alignment preview without rerunning trimming. */
export function reevaluateAlignmentView(
  viewData: AlignmentViewResponse,
  recipe: TrimmingRecipe,
  totalDatasetTaxa: number
): AlignmentViewResponse {
  const { diff } = viewData;
  const failReasons: string[] = [];
  const isCodingWithOrf =
    recipe.enable_orf &&
    !(recipe.exclude_uce && shouldSkipOrfLocus(diff.id, recipe.orf_search_mode));
  const hadOrfFailure = diff.fail_reasons.some(
    (reason) =>
      reason.startsWith('ORF check failed') ||
      reason.startsWith('Candidate ORF failed') ||
      reason.startsWith('Reference-guided exon failed')
  );
  const foundValidOrf = diff.found_valid_orf ?? !hadOrfFailure;

  if (diff.new_taxa_count === 0) {
    failReasons.push('0 surviving taxa (all samples pruned)');
  }

  if (
    isCodingWithOrf &&
    (!foundValidOrf ||
      diff.new_taxa_count === 0 ||
      hadOrfFailure ||
      diff.final_stop_codons?.some((stop) => !stop.is_terminal))
  ) {
    failReasons.push(candidateOrfFailureReason(diff, recipe));
  }

  if (recipe.assess_alignment) {
    if (recipe.min_taxa > 0 && diff.new_taxa_count < recipe.min_taxa) {
      failReasons.push(`Taxa (${diff.new_taxa_count} < min ${recipe.min_taxa})`);
    }
    if (totalDatasetTaxa > 0 && recipe.min_taxa_occupancy_percent > 0) {
      const occupancyPercent = (diff.new_taxa_count / totalDatasetTaxa) * 100;
      if (occupancyPercent < recipe.min_taxa_occupancy_percent) {
        failReasons.push(
          `Taxon occupancy (${occupancyPercent.toFixed(1)}% < min ${recipe.min_taxa_occupancy_percent}%, ${diff.new_taxa_count}/${totalDatasetTaxa} taxa)`
        );
      }
    }
    if (recipe.min_length > 0 && diff.new_length < recipe.min_length) {
      failReasons.push(`Length (${diff.new_length} bp < min ${recipe.min_length} bp)`);
    }
    if (recipe.max_gap_percent > 0 && diff.new_gap_percent > recipe.max_gap_percent) {
      failReasons.push(
        `Gap % (${diff.new_gap_percent.toFixed(1)}% > max ${recipe.max_gap_percent}%)`
      );
    }
    const variablePercent = diff.new_length > 0 ? (diff.new_variable / diff.new_length) * 100 : 0;
    if (
      (recipe.min_variable_count ?? 0) > 0 &&
      diff.new_variable < recipe.min_variable_count
    ) {
      failReasons.push(
        `Variable sites (${diff.new_variable} < min ${recipe.min_variable_count})`
      );
    }
    if (
      (recipe.min_variable_percent ?? 0) > 0 &&
      variablePercent < recipe.min_variable_percent
    ) {
      failReasons.push(
        `Variable-site proportion (${variablePercent.toFixed(1)}% < min ${recipe.min_variable_percent.toFixed(1)}%)`
      );
    }
    const pisPercent = diff.new_length > 0 ? (diff.new_pis / diff.new_length) * 100 : 0;
    if ((recipe.min_pis_count ?? 0) > 0 && diff.new_pis < recipe.min_pis_count) {
      failReasons.push(
        `Parsimony-informative sites (${diff.new_pis} < min ${recipe.min_pis_count})`
      );
    }
    if ((recipe.min_pis_percent ?? 0) > 0 && pisPercent < recipe.min_pis_percent) {
      failReasons.push(
        `Parsimony-informative proportion (${pisPercent.toFixed(1)}% < min ${recipe.min_pis_percent.toFixed(1)}%)`
      );
    }
  }

  return {
    ...viewData,
    diff: {
      ...diff,
      pass: failReasons.length === 0,
      fail_reasons: failReasons,
    },
  };
}

export function executeClientTrimming(
  alignment: Alignment,
  recipe: TrimmingRecipe,
  totalDatasetTaxa: number = 0
): AlignmentViewResponse {
  let curTaxa = [...alignment.taxa];
  let curSeqs = [...alignment.sequences];

  const oldTaxaCount = curTaxa.length;
  const oldLength = alignment.length;

  const droppedTaxa: string[] = [];
  const droppedTaxaReasons: Record<string, string> = {};
  let rawColMap: number[] = Array.from({ length: oldLength }, (_, i) => i);
  let columnReasons: Record<number, string> = {};
  const rawSiteStats = calculateClientSiteStatistics(alignment.sequences);
  const pisMask = rawSiteStats.pisMask;

  // Step 0: Drop samples that fall below the dataset-wide locus occupancy threshold.
  const excludedTaxa = new Set(recipe.excluded_taxa ?? []);
  if (excludedTaxa.size > 0) {
    const keptTaxa: string[] = [];
    const keptSeqs: string[] = [];
    for (let i = 0; i < curTaxa.length; i++) {
      const taxon = curTaxa[i];
      if (excludedTaxa.has(taxon)) {
        droppedTaxa.push(taxon);
        droppedTaxaReasons[taxon] =
          `Dataset-wide sample occupancy (< ${recipe.min_sample_locus_occupancy_percent.toFixed(0)}% of loci)`;
      } else {
        keptTaxa.push(taxon);
        keptSeqs.push(curSeqs[i]);
      }
    }
    curTaxa = keptTaxa;
    curSeqs = keptSeqs;
  }

  // Step 1: Character sanitation. Do not remove columns until the early
  // sample filters have removed every sequence that should not contribute.
  if (recipe.replace_n_with_gap) {
    curSeqs = curSeqs.map((s) => s.replace(/[Nn?]/g, '-'));
  }
  curSeqs = convertClientAmbiguities(curSeqs, recipe.ambiguity_strategy);

  // Step 2: Sample-level filtering before all column-level decisions.
  const coverageAlignmentLength = curSeqs[0]?.length ?? 0;
  if (
    recipe.trim_coverage &&
    curSeqs.length > 2 &&
    !(recipe.min_coverage_bp >= coverageAlignmentLength && coverageAlignmentLength > 0)
  ) {
    const baseCounts = curSeqs.map((sequence) =>
      [...sequence].filter((base) => {
        const upper = base.toUpperCase();
        return upper !== '-' && upper !== '?' && upper !== 'N';
      }).length
    );
    const maxSampleBasePairs = Math.max(1, ...baseCounts);
    const referenceWidth = recipe.relative_width === 'sample'
      ? maxSampleBasePairs
      : Math.max(1, coverageAlignmentLength);
    const keptTaxa: string[] = [];
    const keptSeqs: string[] = [];
    for (let index = 0; index < curTaxa.length; index++) {
      const basePairs = baseCounts[index];
      const coveragePercent = (basePairs / referenceWidth) * 100;
      if (
        basePairs >= recipe.min_coverage_bp &&
        coveragePercent >= recipe.min_coverage_percent
      ) {
        keptTaxa.push(curTaxa[index]);
        keptSeqs.push(curSeqs[index]);
      } else {
        const taxon = curTaxa[index];
        droppedTaxa.push(taxon);
        droppedTaxaReasons[taxon] =
          `Low coverage (< ${recipe.min_coverage_bp} bp / ${recipe.min_coverage_percent.toFixed(0)}%)`;
      }
    }
    curTaxa = keptTaxa;
    curSeqs = keptSeqs;
  }

  if (recipe.trim_similarity && curSeqs.length > 2) {
    const similarityConsensus = calculateClientMajorityConsensus(curSeqs);
    const comparisonLength = curSeqs[0]?.length ?? 0;
    const keptTaxa: string[] = [];
    const keptSeqs: string[] = [];
    for (let index = 0; index < curTaxa.length; index++) {
      const sequence = curSeqs[index];
      let mismatches = 0;
      let overlap = 0;
      for (let column = 0; column < comparisonLength; column++) {
        const state = sequence[column]?.toUpperCase();
        const consensusState = similarityConsensus[column];
        if (state && consensusState && state !== '-' && state !== '?' && state !== 'N' && consensusState !== '-') {
          overlap++;
          if (state !== consensusState) mismatches++;
        }
      }
      const divergence = overlap > 0 ? mismatches / overlap : 0;
      if (divergence < recipe.similarity_threshold) {
        keptTaxa.push(curTaxa[index]);
        keptSeqs.push(sequence);
      } else {
        const taxon = curTaxa[index];
        droppedTaxa.push(taxon);
        droppedTaxaReasons[taxon] =
          `Divergent outlier (> ${(recipe.similarity_threshold * 100).toFixed(0)}% divergence)`;
      }
    }
    curTaxa = keptTaxa;
    curSeqs = keptSeqs;
  }

  // Step 3: Remove missing-only columns using only surviving samples.
  if ((recipe.remove_gap_only_columns ?? true) && curSeqs.length > 0) {
    const sanitizedLength = Math.min(...curSeqs.map((sequence) => sequence.length));
    const keptColumns: number[] = [];
    for (let column = 0; column < sanitizedLength; column++) {
      const hasObservedState = curSeqs.some((sequence) => {
        const base = sequence[column]?.toUpperCase();
        return base !== undefined && base !== '-' && base !== 'N' && base !== '?';
      });
      if (hasObservedState) {
        keptColumns.push(column);
      } else {
        columnReasons[rawColMap[column]] = 'Gap-Only / Missing-Only Sanitation';
      }
    }
    if (keptColumns.length !== sanitizedLength) {
      rawColMap = keptColumns.map((column) => rawColMap[column]);
      curSeqs = curSeqs.map((sequence) =>
        keptColumns.map((column) => sequence[column]).join('')
      );
    }
  }

  // Keep all original taxa on the sanitized column coordinate system so the
  // raw viewer can display stops in trimmed regions and ORF-pruned samples.
  const rawOrfTaxa = [...alignment.taxa];
  const rawOrfSequences = alignment.sequences.map((sequence) =>
    rawColMap.map((column) => sequence[column] ?? '-').join('')
  );
  const rawOrfColMap = [...rawColMap];

  // Viewer consensus uses the same surviving-sample set and sanitized columns.
  const majorityConsensus = calculateClientMajorityConsensus(curSeqs);

  const skipOrf = shouldSkipOrfLocus(alignment.id, recipe.orf_search_mode);
  const isCodingWithOrf = recipe.enable_orf && !(recipe.exclude_uce && skipOrf);
  const usesReferenceMode =
    recipe.orf_search_mode === 'referenceguided' ||
    recipe.orf_search_mode === 'referencecandidateorf';
  const allowsReferenceFallback = recipe.orf_search_mode === 'referencecandidateorf';
  const preReferenceSeqs = [...curSeqs];
  const preReferenceRawColMap = [...rawColMap];
  const preReferenceColumnReasons = { ...columnReasons };
  let orfReferenceEvaluated = false;
  let orfReferenceMatched = false;
  let orfReferenceIdentity = 0;
  let orfReferenceCoverage = 0;
  let orfIntronLength = 0;
  let referenceFrameHint: Pick<ClientReadingFrame, 'isReverse' | 'offset'> | null = null;

  if (isCodingWithOrf && usesReferenceMode && curSeqs.length > 0) {
    orfReferenceEvaluated = true;
    const reference = recipe.orf_reference_sequences?.[alignment.id];
    if (reference) {
      const referenceMatch = matchClientReferenceToAlignment(curSeqs, reference);
      if (referenceMatch) {
        const currentLength = curSeqs[0]?.length ?? 0;
        const start = Math.min(currentLength, referenceMatch.start);
        const end = Math.max(start, Math.min(currentLength, referenceMatch.end));
        if (end > start) {
          orfReferenceMatched = true;
          orfReferenceIdentity = referenceMatch.identity;
          orfReferenceCoverage = referenceMatch.coverage;
          orfIntronLength = currentLength - (end - start);
          const referenceFrame = findClientReadingFrame(
            [normalizeClientReference(reference)],
            recipe.genetic_code
          );
          referenceFrameHint = {
            isReverse: referenceFrame.isReverse !== referenceMatch.isReverse,
            offset: referenceFrame.offset,
          };
          for (let column = 0; column < start; column++) {
            columnReasons[rawColMap[column]] = 'Reference-Anchored Intron';
          }
          for (let column = end; column < currentLength; column++) {
            columnReasons[rawColMap[column]] = 'Reference-Anchored Intron';
          }
          rawColMap = rawColMap.slice(start, end);
          curSeqs = curSeqs.map((sequence) => sequence.slice(start, end));
        }
      }
    }
  }

  // Step 3a: Profile HMM segment cleaner (TAPIR-Style)
  if (recipe.trim_hmm && curSeqs.length > 2) {
    const { cleanedSeqs } = cleanWithProfileHMM(
      curTaxa,
      curSeqs,
      recipe.hmm_min_posterior,
      recipe.hmm_min_segment_length
    );
    curSeqs = cleanedSeqs;
  }

  let foundValidOrf = true;
  let orfEvaluated = false;
  let orfCandidateFound = false;
  let orfFrame: number | undefined;
  let orfStart: number | undefined;
  let orfEnd: number | undefined;
  let orfSupportCount = 0;
  let orfSupportPercent = 0;
  let orfRetainedSamples = 0;
  let orfCandidateLengthAa = 0;
  let orfCodingScore = 0;
  let orfAminoAcidConservation = 0;
  let orfFrameContrast = 0;
  let rawOrfStopCodons: StopCodonPos[] | null = null;

  // Step 3b: ORF optimization & stop codon sample pruning
  const referenceBlocksOrf =
    isCodingWithOrf && recipe.orf_search_mode === 'referenceguided' && !orfReferenceMatched;
  if (referenceBlocksOrf) foundValidOrf = false;
  if (isCodingWithOrf && !referenceBlocksOrf && curSeqs.length > 0) {
    let selectedRegionRawMap = [...rawColMap];
    orfEvaluated = true;
    const useContinuousCandidate =
      referenceFrameHint !== null ||
      recipe.orf_search_mode === 'continuouscds' ||
      recipe.orf_search_mode === 'referenceguided';
    const bestShared = !useContinuousCandidate
      ? findClientBestSharedOrfSegment(
          curSeqs,
          recipe.genetic_code,
          recipe.orf_min_shared_support_percent,
          recipe.orf_min_segment_aa,
          recipe.orf_min_coding_score
        )
      : null;
    let frame = useContinuousCandidate
      ? (() => {
          const continuousFrame = referenceFrameHint
            ? evaluateClientReadingFrame(
                curSeqs,
                recipe.genetic_code,
                referenceFrameHint.isReverse,
                referenceFrameHint.offset
              )
            : findClientReadingFrame(curSeqs, recipe.genetic_code);
          const currentLength = curSeqs[0]?.length || 0;
          const trailingRemainder = Math.max(0, currentLength - continuousFrame.offset) % 3;
          const end = currentLength - trailingRemainder;
          const oriented = continuousFrame.isReverse
            ? curSeqs.map(reverseComplementClient)
            : curSeqs;
          const evidence = clientCodingEvidence(
            oriented,
            continuousFrame.offset,
            end,
            recipe.genetic_code
          );
          return {
            frame: continuousFrame.isReverse
              ? -(continuousFrame.offset + 1)
              : continuousFrame.offset + 1,
            isReverse: continuousFrame.isReverse,
            start: continuousFrame.offset,
            end,
            candidateFound: continuousFrame.cleanTaxaCount > 0,
            found: continuousFrame.cleanTaxaCount > 0,
            supportCount: continuousFrame.cleanTaxaCount,
            lengthAa: Math.floor((end - continuousFrame.offset) / 3),
            codingScore: evidence.score,
            conservation: evidence.conservation,
            contrast: evidence.contrast,
          };
        })()
      : {
          frame: bestShared?.frame,
          isReverse: bestShared?.isReverse ?? false,
          start: bestShared?.start ?? 0,
          end: bestShared?.end ?? (curSeqs[0]?.length || 0),
          candidateFound: bestShared !== null,
          found:
            bestShared !== null &&
            bestShared.codingScore >= recipe.orf_min_coding_score,
          supportCount: bestShared?.supportCount ?? 0,
          lengthAa: bestShared?.lengthCodons ?? 0,
          codingScore: bestShared?.codingScore ?? 0,
          conservation: bestShared?.aminoAcidConservation ?? 0,
          contrast: bestShared?.frameContrast ?? 0,
        };

    let usedReferenceFrame = referenceFrameHint !== null;
    const guidedAttemptFailed =
      usedReferenceFrame &&
      (!frame.found ||
        (frame.isReverse && !recipe.auto_flip_reverse) ||
        (recipe.stop_codon_action === 'keep' && frame.supportCount < curTaxa.length));
    if (allowsReferenceFallback && orfReferenceMatched && guidedAttemptFailed) {
      curSeqs = [...preReferenceSeqs];
      rawColMap = [...preReferenceRawColMap];
      columnReasons = { ...preReferenceColumnReasons };

      if (recipe.trim_hmm && curSeqs.length > 2) {
        const { cleanedSeqs } = cleanWithProfileHMM(
          curTaxa,
          curSeqs,
          recipe.hmm_min_posterior,
          recipe.hmm_min_segment_length
        );
        curSeqs = cleanedSeqs;
      }

      selectedRegionRawMap = [...rawColMap];
      const fallback = findClientBestSharedOrfSegment(
        curSeqs,
        recipe.genetic_code,
        recipe.orf_min_shared_support_percent,
        recipe.orf_min_segment_aa,
        recipe.orf_min_coding_score
      );
      frame = {
        frame: fallback?.frame,
        isReverse: fallback?.isReverse ?? false,
        start: fallback?.start ?? 0,
        end: fallback?.end ?? (curSeqs[0]?.length || 0),
        candidateFound: fallback !== null,
        found: fallback !== null && fallback.codingScore >= recipe.orf_min_coding_score,
        supportCount: fallback?.supportCount ?? 0,
        lengthAa: fallback?.lengthCodons ?? 0,
        codingScore: fallback?.codingScore ?? 0,
        conservation: fallback?.aminoAcidConservation ?? 0,
        contrast: fallback?.frameContrast ?? 0,
      };
      usedReferenceFrame = false;
    }
    orfCandidateFound = frame.candidateFound;
    orfFrame = frame.frame;
    orfSupportCount = frame.supportCount;
    orfSupportPercent = curTaxa.length > 0 ? frame.supportCount / curTaxa.length * 100 : 0;
    orfCandidateLengthAa = frame.lengthAa;
    orfCodingScore = frame.codingScore;
    orfAminoAcidConservation = frame.conservation;
    orfFrameContrast = frame.contrast;
    foundValidOrf = frame.found;
    if (frame.isReverse && !recipe.auto_flip_reverse) {
      foundValidOrf = false;
    }

    const orientedMetadataMap = frame.isReverse ? [...rawColMap].reverse() : rawColMap;
    const mappedCandidate = orientedMetadataMap.slice(frame.start, frame.end);
    if (frame.candidateFound && mappedCandidate.length > 0) {
      orfStart = Math.min(...mappedCandidate);
      orfEnd = Math.max(...mappedCandidate) + 1;
    }
    if (orfFrame !== undefined) {
      rawOrfStopCodons = detectClientRawStopsInFrame(
        rawOrfTaxa,
        rawOrfSequences,
        orfFrame,
        recipe.genetic_code,
        rawOrfColMap,
        selectedRegionRawMap
      );
    }

    if (frame.isReverse && recipe.auto_flip_reverse && foundValidOrf) {
      curSeqs = curSeqs.map(reverseComplementClient);
      rawColMap.reverse();
    }

    const shouldTrim =
      usedReferenceFrame ||
      recipe.auto_shift_frame ||
      recipe.orf_search_mode === 'bestsharedsegment' ||
      recipe.orf_search_mode === 'referencecandidateorf';
    if (shouldTrim && foundValidOrf) {
      const currentLength = curSeqs[0]?.length || 0;
      const start = Math.min(frame.start, currentLength);
      const end = Math.max(start, Math.min(frame.end, currentLength));
      for (let column = 0; column < start; column++) {
        columnReasons[rawColMap[column]] = 'ORF Frame Shift / Codon Boundary';
      }
      for (let column = end; column < currentLength; column++) {
        columnReasons[rawColMap[column]] = 'ORF Frame Shift / Codon Boundary';
      }
      rawColMap = rawColMap.slice(start, end);
      curSeqs = curSeqs.map((sequence) => sequence.slice(start, end));
    }

    const stops = clientStopCodonsForCode(recipe.genetic_code);
    if (foundValidOrf && recipe.stop_codon_action === 'removesample') {
      const keptTaxa: string[] = [];
      const keptSeqs: string[] = [];
      for (let i = 0; i < curTaxa.length; i++) {
        const taxon = curTaxa[i];
        const seq = curSeqs[i];
        let hasInternalStop = false;
        const codonCount = Math.floor(seq.length / 3);
        for (let codonIndex = 0; codonIndex < codonCount - 1; codonIndex++) {
          const start = codonIndex * 3;
          if (stops.has(seq.slice(start, start + 3).toUpperCase())) {
            hasInternalStop = true;
            break;
          }
        }
        if (hasInternalStop) {
          droppedTaxa.push(taxon);
          droppedTaxaReasons[taxon] = 'Premature internal stop codon in exon';
        } else {
          keptTaxa.push(taxon);
          keptSeqs.push(seq);
        }
      }
      curTaxa = keptTaxa;
      curSeqs = keptSeqs;
    } else if (foundValidOrf && recipe.stop_codon_action === 'maskcodon') {
      curSeqs = curSeqs.map((sequence) => {
        const chars = sequence.split('');
        const codonCount = Math.floor(sequence.length / 3);
        for (let codonIndex = 0; codonIndex < codonCount - 1; codonIndex++) {
          const start = codonIndex * 3;
          if (stops.has(sequence.slice(start, start + 3).toUpperCase())) {
            chars[start] = '-';
            chars[start + 1] = '-';
            chars[start + 2] = '-';
          }
        }
        return chars.join('');
      });
    }
    orfRetainedSamples = curTaxa.length;
  }

  // Step 4: External ragged edge trimming
  if (recipe.trim_external && curSeqs.length > 0) {
    const minTaxa = Math.ceil(curSeqs.length * (recipe.min_external_percent / 100));
    let start = 0;
    let end = rawColMap.length;

    for (let col = 0; col < rawColMap.length; col++) {
      let count = 0;
      for (const seq of curSeqs) {
        const c = seq[col];
        if (c && c !== '-') count++;
      }
      if (count >= minTaxa) {
        start = col;
        break;
      }
    }

    for (let col = rawColMap.length - 1; col >= 0; col--) {
      let count = 0;
      for (const seq of curSeqs) {
        const c = seq[col];
        if (c && c !== '-') count++;
      }
      if (count >= minTaxa) {
        end = col + 1;
        break;
      }
    }

    if (recipe.codon_preserving || isCodingWithOrf) {
      const rem = start % 3;
      if (rem === 1) start += 2;
      else if (rem === 2) start += 1;
      end -= Math.max(0, end - start) % 3;
    }

    for (let col = 0; col < start; col++) {
      const rawC = rawColMap[col];
      columnReasons[rawC] = `Ragged 5' End (< ${recipe.min_external_percent.toFixed(0)}% taxa coverage)`;
    }
    for (let col = end; col < rawColMap.length; col++) {
      const rawC = rawColMap[col];
      columnReasons[rawC] = `Ragged 3' End (< ${recipe.min_external_percent.toFixed(0)}% taxa coverage)`;
    }

    rawColMap = rawColMap.slice(start, end);
    curSeqs = curSeqs.map((s) => (start < end ? s.slice(start, end) : ''));
  }

  // Step 5: Column gap filter (Bypassed on coding loci when ORF is enabled to preserve reading frames)
  if (recipe.trim_columns && !isCodingWithOrf && curSeqs.length > 0) {
    const curLen = curSeqs[0]?.length || 0;
    const keptIdxs: number[] = [];
    for (let col = 0; col < curLen; col++) {
      let gaps = 0;
      for (const seq of curSeqs) {
        if (seq[col] === '-' || seq[col] === '?' || (recipe.count_n_as_gap && seq[col]?.toUpperCase() === 'N')) {
          gaps++;
        }
      }
      const gapPct = (gaps / curSeqs.length) * 100;
      if (gapPct < recipe.min_column_gap_percent) {
        keptIdxs.push(col);
      } else {
        const rawC = rawColMap[col];
        columnReasons[rawC] = `High Gap Column (${gapPct.toFixed(1)}% gaps > max ${recipe.min_column_gap_percent.toFixed(1)}%)`;
      }
    }

    rawColMap = keptIdxs.map((idx) => rawColMap[idx]);
    curSeqs = curSeqs.map((s) => keptIdxs.map((idx) => s[idx]).join(''));
  }

  // Step 6: Statistical Column Trimming (trimAl & Gblocks) (Bypassed on coding loci when ORF is enabled)
  if (recipe.enable_statistical_columns && !isCodingWithOrf && curSeqs.length > 0 && curSeqs[0]?.length > 0) {
    const curLen = curSeqs[0].length;
    const keptIdxs: number[] = [];
    const droppedLocalCols: number[] = [];

    const method = recipe.stat_col_method || 'trimalsimilarity';
    const windowSize = recipe.stat_col_window_size || 3;
    const heuristic = recipe.stat_col_heuristic || 'custom';

    if (method === 'trimalsimilarity') {
      const rawScores = Array.from({ length: curLen }, (_, col) =>
        computeClientColumnSimilarity(curSeqs, col)
      );
      const smoothed = applyClientSlidingWindow(rawScores, windowSize);
      const avg = smoothed.reduce((a, b) => a + b, 0) / Math.max(1, curLen);

      let effectiveThresh = recipe.stat_col_similarity_threshold ?? 0.35;
      if (heuristic === 'gappyout') {
        effectiveThresh = Math.max(0.15, Math.min(0.60, avg * 0.75));
      } else if (heuristic === 'strict') {
        effectiveThresh = Math.max(0.25, Math.min(0.70, avg * 0.90));
      } else if (heuristic === 'strictplus') {
        effectiveThresh = Math.max(0.35, Math.min(0.85, avg * 1.10));
      }

      for (let col = 0; col < curLen; col++) {
        const score = smoothed[col];
        const gapFrac = computeClientColumnGapFraction(curSeqs, col);

        let isDropped = false;
        if (heuristic === 'gappyout') {
          isDropped = gapFrac > 0.60 || score < effectiveThresh;
        } else if (heuristic === 'strict' || heuristic === 'strictplus') {
          isDropped = gapFrac > 0.40 || score < effectiveThresh;
        } else {
          isDropped = score < effectiveThresh;
        }

        if (isDropped) {
          droppedLocalCols.push(col);
          const rawC = rawColMap[col];
          const label =
            heuristic === 'gappyout'
              ? `trimAl Gappyout (Score ${score.toFixed(2)}, Gap ${(gapFrac * 100).toFixed(0)}%)`
              : heuristic === 'strict'
              ? `trimAl Strict (Score ${score.toFixed(2)}, Gap ${(gapFrac * 100).toFixed(0)}%)`
              : heuristic === 'strictplus'
              ? `trimAl StrictPlus (Score ${score.toFixed(2)}, Gap ${(gapFrac * 100).toFixed(0)}%)`
              : `trimAl Similarity (${score.toFixed(2)} < min ${effectiveThresh.toFixed(2)}, win=${windowSize})`;
          columnReasons[rawC] = label;
        } else {
          keptIdxs.push(col);
        }
      }
    } else if (method === 'gblocksblocks') {
      const rawScores = Array.from({ length: curLen }, (_, col) =>
        computeClientColumnSimilarity(curSeqs, col)
      );
      const gapFracs = Array.from({ length: curLen }, (_, col) =>
        computeClientColumnGapFraction(curSeqs, col)
      );
      const minBlockLen = recipe.stat_col_min_block_length || 5;
      const gapTreatment = recipe.stat_col_gap_treatment || 'half';

      const isConserved = Array.from({ length: curLen }, (_, col) => {
        const sim = rawScores[col];
        const gap = gapFracs[col];
        const gapOk =
          gapTreatment === 'none' ? gap === 0 : gapTreatment === 'half' ? gap <= 0.5 : true;
        return sim >= 0.5 && gapOk;
      });

      const blockMembership = new Array(curLen).fill(false);
      let col = 0;
      while (col < curLen) {
        if (isConserved[col]) {
          const blockStart = col;
          let blockEnd = col;
          while (blockEnd < curLen && isConserved[blockEnd]) {
            blockEnd++;
          }
          const blockLen = blockEnd - blockStart;
          if (blockLen >= minBlockLen) {
            for (let b = blockStart; b < blockEnd; b++) {
              blockMembership[b] = true;
            }
          }
          col = blockEnd;
        } else {
          col++;
        }
      }

      for (let c = 0; c < curLen; c++) {
        if (blockMembership[c]) {
          keptIdxs.push(c);
        } else {
          droppedLocalCols.push(c);
          const rawC = rawColMap[c];
          columnReasons[rawC] = !isConserved[c]
            ? `Gblocks Non-Conserved / Gap (sim ${rawScores[c].toFixed(2)})`
            : `Gblocks Fragment (< ${minBlockLen} bp conserved block)`;
        }
      }
    } else if (method === 'entropy') {
      const rawEntropy = Array.from({ length: curLen }, (_, col) =>
        computeClientColumnEntropy(curSeqs, col)
      );
      const smoothed = applyClientSlidingWindow(rawEntropy, windowSize);
      const maxEntropy = recipe.stat_col_entropy_threshold ?? 1.5;

      for (let col = 0; col < curLen; col++) {
        const ent = smoothed[col];
        if (ent > maxEntropy) {
          droppedLocalCols.push(col);
          const rawC = rawColMap[col];
          columnReasons[rawC] = `High Entropy / Noise (${ent.toFixed(2)} > max ${maxEntropy.toFixed(2)} bits, win=${windowSize})`;
        } else {
          keptIdxs.push(col);
        }
      }
    } else {
      for (let col = 0; col < curLen; col++) keptIdxs.push(col);
    }

    rawColMap = keptIdxs.map((idx) => rawColMap[idx]);
    curSeqs = curSeqs.map((s) => keptIdxs.map((idx) => s[idx]).join(''));
  }

  const keptSet = new Set(rawColMap);
  const finalTrimmedCols: number[] = [];
  for (let c = 0; c < oldLength; c++) {
    if (!keptSet.has(c)) {
      finalTrimmedCols.push(c);
    }
  }

  const newTaxaCount = curTaxa.length;
  const newLength = curSeqs[0]?.length || 0;
  const oldVariable = rawSiteStats.variableCount;
  const oldPis = rawSiteStats.pisCount;
  const newSiteStats = calculateClientSiteStatistics(curSeqs);
  const newVariable = newSiteStats.variableCount;
  const newPis = newSiteStats.pisCount;
  const oldGapPercent = calculateClientGapPercent(alignment.sequences);
  const newGapPercent = calculateClientGapPercent(curSeqs);
  const finalStopCodons = detectClientStopCodons(curTaxa, curSeqs, recipe.genetic_code);
  const mappedFinalStopCodons = finalStopCodons.flatMap((stop) => {
    const mapped = rawColMap.slice(stop.start, stop.end);
    if (mapped.length !== 3) return [];
    const rawStart = Math.min(...mapped);
    const rawEnd = Math.max(...mapped) + 1;
    return rawEnd - rawStart === 3 ? [{ ...stop, start: rawStart, end: rawEnd }] : [];
  });
  const stopCodons = rawOrfStopCodons ?? mappedFinalStopCodons;

  const trimmed_alignment: Alignment = {
    id: alignment.id,
    file_name: alignment.file_name,
    file_path: alignment.file_path,
    format: alignment.format,
    taxa: curTaxa,
    sequences: curSeqs,
    length: newLength,
    num_taxa: newTaxaCount,
  };

  const diff: TrimmingDiff = {
    id: alignment.id,
    old_taxa_count: oldTaxaCount,
    new_taxa_count: newTaxaCount,
    dropped_taxa: droppedTaxa,
    old_length: oldLength,
    new_length: newLength,
    trimmed_columns: finalTrimmedCols,
    masked_segments: [],
    column_reasons: columnReasons,
    dropped_taxa_reasons: droppedTaxaReasons,
    stop_codons: stopCodons,
    final_stop_codons: finalStopCodons,
    old_gap_percent: oldGapPercent,
    new_gap_percent: newGapPercent,
    old_variable: oldVariable,
    new_variable: newVariable,
    old_pis: oldPis,
    new_pis: newPis,
    found_valid_orf: foundValidOrf,
    orf_evaluated: orfEvaluated,
    orf_candidate_found: orfCandidateFound,
    orf_frame: orfFrame,
    orf_start: orfStart,
    orf_end: orfEnd,
    orf_support_count: orfSupportCount,
    orf_support_percent: orfSupportPercent,
    orf_retained_samples: orfRetainedSamples,
    orf_candidate_length_aa: orfCandidateLengthAa,
    orf_coding_score: orfCodingScore,
    orf_amino_acid_conservation: orfAminoAcidConservation,
    orf_frame_contrast: orfFrameContrast,
    orf_reference_evaluated: orfReferenceEvaluated,
    orf_reference_matched: orfReferenceMatched,
    orf_reference_identity: orfReferenceIdentity,
    orf_reference_coverage: orfReferenceCoverage,
    orf_intron_length: orfIntronLength,
    pass: (() => {
      if (newTaxaCount === 0) return false;
      const retainedInternalStop = finalStopCodons.some((stop) => !stop.is_terminal);
      if (
        isCodingWithOrf &&
        (!foundValidOrf || newTaxaCount === 0 || retainedInternalStop)
      ) return false;
      if (recipe.assess_alignment) {
        if (recipe.min_taxa > 0 && newTaxaCount < recipe.min_taxa) return false;
        const occupancyDenominator = totalDatasetTaxa > 0 ? totalDatasetTaxa : oldTaxaCount;
        if (occupancyDenominator > 0 && recipe.min_taxa_occupancy_percent > 0) {
          const occPct = (newTaxaCount / occupancyDenominator) * 100;
          if (occPct < recipe.min_taxa_occupancy_percent) return false;
        }
        if (recipe.min_length > 0 && newLength < recipe.min_length) return false;
        if (recipe.max_gap_percent > 0 && newGapPercent > recipe.max_gap_percent) return false;
        if ((recipe.min_variable_count ?? 0) > 0 && newVariable < recipe.min_variable_count) return false;
        if ((recipe.min_variable_percent ?? 0) > 0 && newSiteStats.variablePercent < recipe.min_variable_percent) return false;
        if ((recipe.min_pis_count ?? 0) > 0 && newPis < recipe.min_pis_count) return false;
        if ((recipe.min_pis_percent ?? 0) > 0 && newSiteStats.pisPercent < recipe.min_pis_percent) return false;
      }
      return true;
    })(),
    fail_reasons: (() => {
      const reasons: string[] = [];
      if (newTaxaCount === 0) {
        reasons.push('0 surviving taxa (all samples pruned)');
      }
      const retainedInternalStop = finalStopCodons.some((stop) => !stop.is_terminal);
      if (
        isCodingWithOrf &&
        (!foundValidOrf || newTaxaCount === 0 || retainedInternalStop)
      ) {
        if (referenceBlocksOrf) {
          reasons.push(
            'Reference-guided exon failed: no matching reference sequence could be anchored for this locus'
          );
        } else if (orfEvaluated && !orfCandidateFound) {
          reasons.push(
            recipe.orf_search_mode === 'continuouscds' ||
            recipe.orf_search_mode === 'referenceguided'
              ? 'Candidate ORF failed: no usable continuous reading frame was found'
              : `Candidate ORF failed: no segment met ${recipe.orf_min_shared_support_percent.toFixed(0)}% support and ${recipe.orf_min_segment_aa} aa minimums`
          );
        } else if ((orfFrame ?? 0) < 0 && !recipe.auto_flip_reverse) {
          reasons.push(
            'Candidate ORF failed: the best candidate is reverse-strand but Auto-Flip Reverse Strand is disabled'
          );
        } else if (
          orfEvaluated &&
          orfCandidateFound &&
          (recipe.orf_search_mode === 'bestsharedsegment' ||
            recipe.orf_search_mode === 'referencecandidateorf') &&
          orfCodingScore < recipe.orf_min_coding_score
        ) {
          reasons.push(
            `Candidate ORF failed: coding evidence ${orfCodingScore.toFixed(1)} < minimum ${recipe.orf_min_coding_score.toFixed(1)}`
          );
        } else {
          reasons.push(
            'Candidate ORF failed: no retained samples or premature stop codons remain'
          );
        }
      }
      if (recipe.assess_alignment) {
        if (recipe.min_taxa > 0 && newTaxaCount < recipe.min_taxa) {
          reasons.push(`Taxa (${newTaxaCount} < min ${recipe.min_taxa})`);
        }
        const occupancyDenominator = totalDatasetTaxa > 0 ? totalDatasetTaxa : oldTaxaCount;
        if (occupancyDenominator > 0 && recipe.min_taxa_occupancy_percent > 0) {
          const occPct = (newTaxaCount / occupancyDenominator) * 100;
          if (occPct < recipe.min_taxa_occupancy_percent) {
            reasons.push(
              `Taxon occupancy (${occPct.toFixed(1)}% < min ${recipe.min_taxa_occupancy_percent}%, ${newTaxaCount}/${occupancyDenominator} taxa)`
            );
          }
        }
        if (recipe.min_length > 0 && newLength < recipe.min_length) {
          reasons.push(`Length (${newLength} bp < min ${recipe.min_length} bp)`);
        }
        if (recipe.max_gap_percent > 0 && newGapPercent > recipe.max_gap_percent) {
          reasons.push(`Gap % (${newGapPercent.toFixed(1)}% > max ${recipe.max_gap_percent}%)`);
        }
        if ((recipe.min_variable_count ?? 0) > 0 && newVariable < recipe.min_variable_count) {
          reasons.push(
            `Variable sites (${newVariable} < min ${recipe.min_variable_count})`
          );
        }
        if (
          (recipe.min_variable_percent ?? 0) > 0 &&
          newSiteStats.variablePercent < recipe.min_variable_percent
        ) {
          reasons.push(
            `Variable-site proportion (${newSiteStats.variablePercent.toFixed(1)}% < min ${recipe.min_variable_percent.toFixed(1)}%)`
          );
        }
        if ((recipe.min_pis_count ?? 0) > 0 && newPis < recipe.min_pis_count) {
          reasons.push(
            `Parsimony-informative sites (${newPis} < min ${recipe.min_pis_count})`
          );
        }
        if ((recipe.min_pis_percent ?? 0) > 0 && newSiteStats.pisPercent < recipe.min_pis_percent) {
          reasons.push(
            `Parsimony-informative proportion (${newSiteStats.pisPercent.toFixed(1)}% < min ${recipe.min_pis_percent.toFixed(1)}%)`
          );
        }
      }
      return reasons;
    })(),
  };

  return {
    raw_alignment: alignment,
    trimmed_alignment,
    diff,
    pis_mask: pisMask,
    majority_consensus: majorityConsensus,
  };
}

export function cleanWithProfileHMM(
  taxa: string[],
  sequences: string[],
  minPosterior: number,
  minSegLen: number
): { cleanedSeqs: string[]; maskedSegments: { taxon: string; start: number; end: number }[] } {
  if (sequences.length < 3 || sequences[0].length === 0) {
    return { cleanedSeqs: sequences, maskedSegments: [] };
  }

  const length = sequences[0].length;
  const alpha = 0.5;
  const nullProb = 0.25;

  const colCounts: number[][] = Array.from({ length }, () => [0, 0, 0, 0]);
  const colTotals: number[] = new Array(length).fill(0);
  const baseMap: Record<string, number> = { A: 0, C: 1, G: 2, T: 3, U: 3 };

  for (let col = 0; col < length; col++) {
    for (const seq of sequences) {
      const b = seq[col]?.toUpperCase();
      const idx = baseMap[b];
      if (idx !== undefined) {
        colCounts[col][idx]++;
        colTotals[col]++;
      }
    }
  }

  const cleanedSeqs: string[] = [];
  const maskedSegments: { taxon: string; start: number; end: number }[] = [];

  for (let t = 0; t < sequences.length; t++) {
    const seq = sequences[t];
    const chars = seq.split('');
    const confidences: number[] = new Array(length).fill(1.0);

    for (let col = 0; col < length; col++) {
      const b = seq[col]?.toUpperCase();
      const idx = baseMap[b];
      if (idx !== undefined) {
        const countWithoutI = Math.max(0, colCounts[col][idx] - 1);
        const totalWithoutI = Math.max(0, colTotals[col] - 1);
        const denom = totalWithoutI + alpha;
        const pEmit = (countWithoutI + alpha * nullProb) / denom;
        confidences[col] = pEmit / (pEmit + nullProb);
      }
    }

    let segStart: number | null = null;
    for (let col = 0; col < length; col++) {
      const isLow = confidences[col] < minPosterior;
      if (isLow) {
        if (segStart === null) segStart = col;
      } else {
        if (segStart !== null) {
          if (col - segStart >= minSegLen) {
            for (let c = segStart; c < col; c++) chars[c] = '-';
            maskedSegments.push({ taxon: taxa[t], start: segStart, end: col });
          }
          segStart = null;
        }
      }
    }

    if (segStart !== null && length - segStart >= minSegLen) {
      for (let c = segStart; c < length; c++) chars[c] = '-';
      maskedSegments.push({ taxon: taxa[t], start: segStart, end: length });
    }

    cleanedSeqs.push(chars.join(''));
  }

  return { cleanedSeqs, maskedSegments };
}

export function computeClientColumnSimilarity(seqs: string[], col: number): number {
  const chars: string[] = [];
  for (const s of seqs) {
    const c = s[col]?.toUpperCase();
    if (c && c !== '-' && c !== '?' && c !== 'N') {
      chars.push(c);
    }
  }
  const n = chars.length;
  if (n <= 1) return 0.0;

  const counts: Record<string, number> = {};
  for (const c of chars) {
    counts[c] = (counts[c] || 0) + 1;
  }

  let matches = 0;
  for (const cnt of Object.values(counts)) {
    if (cnt >= 2) matches += (cnt * (cnt - 1)) / 2;
  }
  const totalPairs = (n * (n - 1)) / 2;
  return matches / totalPairs;
}

export function applyClientSlidingWindow(scores: number[], windowSize: number): number[] {
  const len = scores.length;
  if (len === 0 || windowSize <= 1) return [...scores];

  const half = Math.floor(windowSize / 2);
  const smoothed: number[] = [];
  for (let i = 0; i < len; i++) {
    const start = Math.max(0, i - half);
    const end = Math.min(len, i + half + 1);
    const slice = scores.slice(start, end);
    const sum = slice.reduce((a, b) => a + b, 0);
    smoothed.push(sum / slice.length);
  }
  return smoothed;
}

export function computeClientColumnEntropy(seqs: string[], col: number): number {
  const counts: Record<string, number> = {};
  let total = 0;
  for (const s of seqs) {
    const u = s[col]?.toUpperCase();
    if (u === 'A' || u === 'C' || u === 'G' || u === 'T' || u === 'U') {
      counts[u] = (counts[u] || 0) + 1;
      total++;
    }
  }
  if (total <= 1) return 0.0;

  let entropy = 0.0;
  for (const cnt of Object.values(counts)) {
    const p = cnt / total;
    entropy -= p * Math.log2(p);
  }
  return entropy;
}

export function computeClientColumnGapFraction(seqs: string[], col: number): number {
  if (seqs.length === 0) return 1.0;
  let gaps = 0;
  for (const s of seqs) {
    const u = s[col]?.toUpperCase();
    if (!u || u === '-' || u === '?' || u === 'N') gaps++;
  }
  return gaps / seqs.length;
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
