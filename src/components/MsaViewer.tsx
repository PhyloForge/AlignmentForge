import React, { useEffect, useRef, useState, useMemo, useCallback } from 'react';
import {
  ZoomIn,
  ZoomOut,
  Palette,
  Eye,
  Sliders,
  Sparkles,
  Settings2,
  Check,
  MoveHorizontal,
  Info,
  CheckCircle,
  XCircle,
  Languages,
} from 'lucide-react';
import {
  AlignmentViewResponse,
  AminoAcidColorScheme,
  AminoAcidViewerSettings,
  ColorScheme,
  GeneticCode,
  StopCodonPos,
} from '../types';
import { translateClientCodon } from '../parsers/clientParser';

interface MsaViewerProps {
  viewData: AlignmentViewResponse;
  colorScheme: ColorScheme;
  onChangeColorScheme: (scheme: ColorScheme) => void;
  showDiffOverlay: boolean;
  onToggleDiffOverlay: () => void;
  geneticCode: GeneticCode;
  aminoAcidViewerSettings: AminoAcidViewerSettings;
  onChangeAminoAcidViewerSettings: (settings: AminoAcidViewerSettings) => void;
}

export type OverlayColorOption = 'grey' | 'slate' | 'amber' | 'red' | 'dark';

// Nucleotide and Gap Palette Definitions
const NUC_COLORS: Record<string, string> = {
  A: '#4caf50', // Green
  C: '#2196f3', // Blue
  G: '#ff9800', // Orange
  T: '#f44336', // Red
  U: '#f44336',
  '-': '#2b3545', // Muted gap
  '?': '#374151',
  N: '#4b5563',
};

// Okabe-Ito-inspired colors avoid relying on red/green discrimination.
const COLORBLIND_NUC_COLORS: Record<string, string> = {
  A: '#0072b2', // Blue
  C: '#e69f00', // Orange
  G: '#009e73', // Bluish green
  T: '#cc79a7', // Reddish purple
  U: '#cc79a7',
  '-': '#2b3545',
  '?': '#374151',
  N: '#4b5563',
};

const HIGH_CONTRAST_NUC_COLORS: Record<string, string> = {
  A: '#1d4ed8',
  C: '#b45309',
  G: '#047857',
  T: '#9d174d',
  U: '#9d174d',
  '-': '#111827',
  '?': '#374151',
  N: '#4b5563',
};

const PASTEL_NUC_COLORS: Record<string, string> = {
  A: '#a7f3d0',
  C: '#bfdbfe',
  G: '#fde68a',
  T: '#fbcfe8',
  U: '#fbcfe8',
  '-': '#d1d5db',
  '?': '#cbd5e1',
  N: '#cbd5e1',
};

const OCEAN_NUC_COLORS: Record<string, string> = {
  A: '#0891b2',
  C: '#2563eb',
  G: '#7c3aed',
  T: '#db2777',
  U: '#db2777',
  '-': '#1e293b',
  '?': '#475569',
  N: '#64748b',
};

const SUNSET_NUC_COLORS: Record<string, string> = {
  A: '#f59e0b',
  C: '#f97316',
  G: '#e11d48',
  T: '#9333ea',
  U: '#9333ea',
  '-': '#292524',
  '?': '#57534e',
  N: '#78716c',
};

const FOREST_NUC_COLORS: Record<string, string> = {
  A: '#65a30d',
  C: '#059669',
  G: '#0d9488',
  T: '#a16207',
  U: '#a16207',
  '-': '#1c2e27',
  '?': '#3f5148',
  N: '#607269',
};

const MONOCHROME_NUC_COLORS: Record<string, string> = {
  A: '#263141',
  C: '#263141',
  G: '#263141',
  T: '#263141',
  U: '#263141',
  '-': '#171b22',
  '?': '#303946',
  N: '#3b4554',
};

const NUCLEOTIDE_PALETTES: Partial<
  Record<ColorScheme, { colors: Record<string, string>; textColor: string }>
> = {
  nucleotide: { colors: NUC_COLORS, textColor: '#ffffff' },
  'nucleotide-colorblind': { colors: COLORBLIND_NUC_COLORS, textColor: '#0b1220' },
  'nucleotide-high-contrast': { colors: HIGH_CONTRAST_NUC_COLORS, textColor: '#ffffff' },
  'nucleotide-pastel': { colors: PASTEL_NUC_COLORS, textColor: '#111827' },
  'nucleotide-ocean': { colors: OCEAN_NUC_COLORS, textColor: '#ffffff' },
  'nucleotide-sunset': { colors: SUNSET_NUC_COLORS, textColor: '#ffffff' },
  'nucleotide-forest': { colors: FOREST_NUC_COLORS, textColor: '#ffffff' },
  monochrome: { colors: MONOCHROME_NUC_COLORS, textColor: '#d1d5db' },
};

const CLUSTAL_COLORS: Record<string, string> = {
  A: '#80a0f0',
  C: '#f08080',
  G: '#90ee90',
  T: '#e066ff',
  '-': '#232833',
  '?': '#374151',
  N: '#4b5563',
};

const CHEMISTRY_AMINO_ACID_COLORS: Record<string, string> = {
  A: '#6b8eec', V: '#6b8eec', I: '#6b8eec', L: '#6b8eec', M: '#6b8eec',
  F: '#d97706', W: '#d97706', Y: '#d97706',
  K: '#2563eb', R: '#2563eb', H: '#2563eb',
  D: '#dc2626', E: '#dc2626',
  S: '#16a34a', T: '#16a34a', N: '#16a34a', Q: '#16a34a',
  G: '#9333ea', P: '#9333ea',
  C: '#ca8a04',
  '*': '#000000',
  X: '#4b5563', '?': '#4b5563', '-': '#232833',
};

const COLORBLIND_AMINO_ACID_COLORS: Record<string, string> = {
  A: '#0072b2', V: '#0072b2', I: '#0072b2', L: '#0072b2', M: '#0072b2',
  F: '#e69f00', W: '#e69f00', Y: '#e69f00',
  K: '#56b4e9', R: '#56b4e9', H: '#56b4e9',
  D: '#d55e00', E: '#d55e00',
  S: '#009e73', T: '#009e73', N: '#009e73', Q: '#009e73',
  G: '#cc79a7', P: '#cc79a7', C: '#f0e442',
  '*': '#000000', X: '#4b5563', '?': '#4b5563', '-': '#232833',
};

const ZAPPO_AMINO_ACID_COLORS: Record<string, string> = {
  A: '#f59e8b', V: '#f59e8b', I: '#f59e8b', L: '#f59e8b', M: '#f59e8b',
  F: '#f59e0b', W: '#f59e0b', Y: '#f59e0b',
  K: '#3b82f6', R: '#3b82f6', H: '#3b82f6',
  D: '#ef4444', E: '#ef4444',
  S: '#22c55e', T: '#22c55e', N: '#22c55e', Q: '#22c55e',
  G: '#e879f9', P: '#38bdf8', C: '#eab308',
  '*': '#000000', X: '#4b5563', '?': '#4b5563', '-': '#232833',
};

const TAYLOR_AMINO_ACID_COLORS: Record<string, string> = {
  A: '#84cc16', C: '#eab308', D: '#ef4444', E: '#f43f5e', F: '#10b981',
  G: '#f97316', H: '#3b82f6', I: '#22c55e', K: '#8b5cf6', L: '#16a34a',
  M: '#059669', N: '#c026d3', P: '#f59e0b', Q: '#db2777', R: '#4f46e5',
  S: '#fb7185', T: '#fb923c', V: '#65a30d', W: '#06b6d4', Y: '#14b8a6',
  '*': '#000000', X: '#4b5563', '?': '#4b5563', '-': '#232833',
};

const HYDROPHOBICITY_AMINO_ACID_COLORS: Record<string, string> = {
  A: '#2563eb', V: '#1d4ed8', I: '#1e40af', L: '#1e3a8a', M: '#3b82f6',
  F: '#4338ca', W: '#3730a3', Y: '#4f46e5',
  K: '#f59e0b', R: '#d97706', H: '#fbbf24',
  D: '#dc2626', E: '#ef4444',
  S: '#10b981', T: '#14b8a6', N: '#22c55e', Q: '#059669',
  G: '#8b5cf6', P: '#a855f7', C: '#ca8a04',
  '*': '#000000', X: '#4b5563', '?': '#4b5563', '-': '#232833',
};

const MONOCHROME_AMINO_ACID_COLORS: Record<string, string> = Object.fromEntries(
  'ACDEFGHIKLMNPQRSTVWY'.split('').map((residue) => [residue, '#263141'])
);
Object.assign(MONOCHROME_AMINO_ACID_COLORS, {
  '*': '#000000', X: '#3b4554', '?': '#3b4554', '-': '#171b22',
});

const AMINO_ACID_PALETTES: Record<
  AminoAcidColorScheme,
  { colors: Record<string, string>; textColor: string }
> = {
  chemistry: { colors: CHEMISTRY_AMINO_ACID_COLORS, textColor: '#ffffff' },
  colorblind: { colors: COLORBLIND_AMINO_ACID_COLORS, textColor: '#ffffff' },
  zappo: { colors: ZAPPO_AMINO_ACID_COLORS, textColor: '#ffffff' },
  taylor: { colors: TAYLOR_AMINO_ACID_COLORS, textColor: '#ffffff' },
  hydrophobicity: { colors: HYDROPHOBICITY_AMINO_ACID_COLORS, textColor: '#ffffff' },
  monochrome: { colors: MONOCHROME_AMINO_ACID_COLORS, textColor: '#d1d5db' },
};

function translateViewerSequence(
  sequence: string,
  geneticCode: GeneticCode,
  hideTerminalStops: boolean
): string {
  const residues: string[] = [];
  for (let start = 0; start + 2 < sequence.length; start += 3) {
    const codon = sequence.slice(start, start + 3);
    residues.push(
      codon === '---'
        ? '-'
        : codon.includes('-')
          ? 'X'
          : translateClientCodon(codon, geneticCode)
    );
  }
  if (hideTerminalStops) {
    for (let index = residues.length - 1; index >= 0; index--) {
      if (residues[index] === '-') continue;
      if (residues[index] === '*') residues[index] = '-';
      break;
    }
  }
  return residues.join('');
}

function calculateViewerConsensus(sequences: string[]): string {
  const length = sequences.reduce((maximum, sequence) => Math.max(maximum, sequence.length), 0);
  let consensus = '';
  for (let column = 0; column < length; column++) {
    const counts = new Map<string, number>();
    for (const sequence of sequences) {
      const state = sequence[column]?.toUpperCase();
      if (!state || state === '-' || state === '?' || state === 'X') continue;
      counts.set(state, (counts.get(state) ?? 0) + 1);
    }
    let bestState = '-';
    let bestCount = 0;
    for (const [state, count] of counts) {
      if (count > bestCount) {
        bestState = state;
        bestCount = count;
      }
    }
    consensus += bestState;
  }
  return consensus;
}

const OVERLAY_COLORS: Record<OverlayColorOption, { fill: string; stroke: string; label: string }> = {
  grey: {
    fill: 'rgba(71, 85, 105, 0.78)', // Muted Slate Grey (Default)
    stroke: '#64748b',
    label: 'Muted Grey (Default)',
  },
  slate: {
    fill: 'rgba(30, 41, 59, 0.88)', // Dark Slate
    stroke: '#475569',
    label: 'Dark Slate',
  },
  amber: {
    fill: 'rgba(217, 119, 6, 0.45)', // Amber
    stroke: '#d97706',
    label: 'Subtle Amber',
  },
  red: {
    fill: 'rgba(239, 68, 68, 0.4)', // Red
    stroke: '#ef4444',
    label: 'Classic Red',
  },
  dark: {
    fill: 'rgba(15, 23, 42, 0.94)', // Dark Dim
    stroke: '#1e293b',
    label: 'Dark Dim',
  },
};

export const MsaViewer: React.FC<MsaViewerProps> = ({
  viewData,
  colorScheme,
  onChangeColorScheme,
  showDiffOverlay,
  onToggleDiffOverlay,
  geneticCode,
  aminoAcidViewerSettings,
  onChangeAminoAcidViewerSettings,
}) => {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const containerRef = useRef<HTMLDivElement | null>(null);

  const [charWidth, setCharWidth] = useState<number>(12); // px per character
  const [charHeight, setCharHeight] = useState<number>(18); // px per row
  const [scrollX, setScrollX] = useState<number>(0);
  const [scrollY, setScrollY] = useState<number>(0);

  // Mouse Drag Panning State
  const [isDragging, setIsDragging] = useState<boolean>(false);
  const lastMousePosRef = useRef<{ x: number; y: number } | null>(null);

  // User-customizable Taxa Column Width
  const [customTaxaWidth, setCustomTaxaWidth] = useState<number | null>(null);
  const [isResizingDivider, setIsResizingDivider] = useState<boolean>(false);
  const dragStartXRef = useRef<number>(0);
  const dragStartWidthRef = useRef<number>(0);

  // View Mode: 'overlays' (raw + ghost diffs) vs 'final' (actual trimmed alignment only)
  const [viewMode, setViewMode] = useState<'overlays' | 'final'>('overlays');
  const [overlayColor, setOverlayColor] = useState<OverlayColorOption>('grey');
  const [showColorMenu, setShowColorMenu] = useState<boolean>(false);
  const [overlayMenuPosition, setOverlayMenuPosition] = useState({ left: 0, top: 0 });

  const [hoverInfo, setHoverInfo] = useState<{
    taxon: string;
    col: number;
    base: string;
    consensus: string;
    isPis: boolean;
    isTrimmed: boolean;
    trimReason: string | null;
    isDroppedTaxon: boolean;
    droppedReason: string | null;
    isStopCodon: boolean;
    stopCodonSeq: string | null;
    screenX: number;
    screenY: number;
  } | null>(null);

  const { raw_alignment, trimmed_alignment, diff, pis_mask, majority_consensus } = viewData;

  const aminoAcidAvailable = Boolean(
    diff.orf_evaluated &&
      diff.found_valid_orf &&
      trimmed_alignment.sequences.length > 0 &&
      trimmed_alignment.length >= 3
  );
  const aminoAcidMode = aminoAcidViewerSettings.enabled && aminoAcidAvailable;
  const effectiveViewMode = aminoAcidMode ? 'final' : viewMode;
  const activeAminoAcidPalette =
    AMINO_ACID_PALETTES[aminoAcidViewerSettings.colorScheme ?? 'chemistry'];

  // Protein translation always uses the final processed alignment because its
  // strand, frame, and codon boundaries have already been resolved by ORF QC.
  const activeTaxa =
    effectiveViewMode === 'final' ? trimmed_alignment.taxa : raw_alignment.taxa;
  const nucleotideActiveSeqs =
    effectiveViewMode === 'final' ? trimmed_alignment.sequences : raw_alignment.sequences;
  const activeSeqs = useMemo(
    () =>
      aminoAcidMode
        ? nucleotideActiveSeqs.map((sequence) =>
            translateViewerSequence(
              sequence,
              geneticCode,
              aminoAcidViewerSettings.hideTerminalStops
            )
          )
        : nucleotideActiveSeqs,
    [
      aminoAcidMode,
      nucleotideActiveSeqs,
      geneticCode,
      aminoAcidViewerSettings.hideTerminalStops,
    ]
  );
  const numTaxa = activeTaxa.length;
  const length = activeSeqs[0]?.length ?? 0;
  const displayConsensus = useMemo(
    () =>
      aminoAcidMode || effectiveViewMode === 'final'
        ? calculateViewerConsensus(activeSeqs)
        : majority_consensus,
    [aminoAcidMode, effectiveViewMode, activeSeqs, majority_consensus]
  );

  useEffect(() => {
    setScrollX(0);
    setScrollY(0);
    setHoverInfo(null);
  }, [aminoAcidMode, effectiveViewMode, raw_alignment.id]);

  useEffect(() => {
    if (effectiveViewMode !== 'overlays') setShowColorMenu(false);
  }, [effectiveViewMode]);

  // Stop codons list
  const stopCodonsList = useMemo(() => {
    if (aminoAcidMode) return [];
    const annotatedStops =
      effectiveViewMode === 'final' ? diff.final_stop_codons : diff.stop_codons;
    if (annotatedStops !== undefined) {
      return annotatedStops;
    }
    const list: StopCodonPos[] = [];
    for (let r = 0; r < activeTaxa.length; r++) {
      const taxon = activeTaxa[r];
      const seq = activeSeqs[r] || '';
      let col = 0;
      while (col + 2 < seq.length) {
        const triplet = (seq[col] + seq[col + 1] + seq[col + 2]).toUpperCase();
        if (['TAA', 'TAG', 'TGA', 'UAA', 'UAG', 'UGA'].includes(triplet)) {
          const isTerminal =
            col + 3 >= seq.length ||
            seq
              .slice(col + 3)
              .split('')
              .every((c) => c === '-' || c === '?');
          list.push({ taxon, start: col, end: col + 3, codon: triplet, is_terminal: isTerminal });
        }
        col += 3;
      }
    }
    return list;
  }, [
    aminoAcidMode,
    effectiveViewMode,
    diff.stop_codons,
    diff.final_stop_codons,
    activeTaxa,
    activeSeqs,
  ]);

  // Auto-calculated Taxon Label Width based on longest taxon name
  const autoTaxaWidth = useMemo(() => {
    const maxLen = Math.max(...activeTaxa.map((t) => t.length), 10);
    return Math.min(550, Math.max(220, Math.ceil(maxLen * 8.2 + 38)));
  }, [activeTaxa]);

  const taxaNameWidth = customTaxaWidth ?? autoTaxaWidth;
  const topHeaderHeight = 44; // Ruler + PIS + Consensus row

  // Handle Drag Resizing of Taxon Column Divider
  const handleDividerMouseDown = (e: React.MouseEvent) => {
    e.preventDefault();
    setIsResizingDivider(true);
    dragStartXRef.current = e.clientX;
    dragStartWidthRef.current = taxaNameWidth;
  };

  useEffect(() => {
    const handleMouseMove = (e: MouseEvent) => {
      if (!isResizingDivider) return;
      const deltaX = e.clientX - dragStartXRef.current;
      const newWidth = Math.max(220, Math.min(800, dragStartWidthRef.current + deltaX));
      setCustomTaxaWidth(newWidth);
    };

    const handleMouseUp = () => {
      setIsResizingDivider(false);
    };

    if (isResizingDivider) {
      window.addEventListener('mousemove', handleMouseMove);
      window.addEventListener('mouseup', handleMouseUp);
    }
    return () => {
      window.removeEventListener('mousemove', handleMouseMove);
      window.removeEventListener('mouseup', handleMouseUp);
    };
  }, [isResizingDivider]);

  // Canvas Render Loop
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const dpr = window.devicePixelRatio || 1;
    const width = canvas.clientWidth;
    const height = canvas.clientHeight;

    canvas.width = width * dpr;
    canvas.height = height * dpr;
    ctx.scale(dpr, dpr);

    // Clear background
    ctx.fillStyle = '#0e1014';
    ctx.fillRect(0, 0, width, height);

    if (numTaxa === 0 || length === 0) return;

    // Viewport calculations
    const visibleCols = Math.ceil((width - taxaNameWidth) / charWidth) + 2;
    const startCol = Math.max(0, Math.floor(scrollX / charWidth));
    const endCol = Math.min(length, startCol + visibleCols);

    const visibleRows = Math.ceil((height - topHeaderHeight) / charHeight) + 2;
    const startRow = Math.max(0, Math.floor(scrollY / charHeight));
    const endRow = Math.min(numTaxa, startRow + visibleRows);

    const trimmedColsSet = new Set(diff.trimmed_columns);
    const droppedTaxaSet = new Set(diff.dropped_taxa);
    const activeOverlay = OVERLAY_COLORS[overlayColor];

    // 1. Draw Sequences & Matrix Cells
    for (let r = startRow; r < endRow; r++) {
      const taxonName = activeTaxa[r];
      const seq = activeSeqs[r];
      const isDropped = effectiveViewMode === 'overlays' && droppedTaxaSet.has(taxonName);
      const rowY = topHeaderHeight + (r * charHeight - scrollY);

      for (let c = startCol; c < endCol; c++) {
        const colX = taxaNameWidth + (c * charWidth - scrollX);
        const char = seq ? seq[c] || '-' : '-';
        const uChar = char.toUpperCase();
        const consChar = displayConsensus[c] || '-';
        const isTrimmedCol =
          effectiveViewMode === 'overlays' && showDiffOverlay && trimmedColsSet.has(c);

        // Check if part of a stop codon
        const isStopCodon = stopCodonsList.some(
          (sc) => sc.taxon === taxonName && c >= sc.start && c < sc.end
        );

        // Color cell background or text
        let fillColor = '#171b22';
        let textColor = '#c9d1d9';

        const nucleotidePalette = NUCLEOTIDE_PALETTES[colorScheme];
        const matchesConsensus = uChar === consChar;
        if (aminoAcidMode) {
          fillColor = activeAminoAcidPalette.colors[uChar] || '#4b5563';
          textColor = activeAminoAcidPalette.textColor;
          if (aminoAcidViewerSettings.dimConsensusMatches && matchesConsensus) {
            fillColor = '#14171d';
            textColor = '#6e7681';
          }
        } else if (nucleotidePalette) {
          fillColor = nucleotidePalette.colors[uChar] || '#1f242e';
          textColor = uChar === '-' ? '#4b5563' : nucleotidePalette.textColor;
        } else if (colorScheme === 'clustalx') {
          fillColor = CLUSTAL_COLORS[uChar] || '#1f242e';
          textColor = '#ffffff';
        } else if (colorScheme === 'difference') {
          if (uChar === consChar) {
            fillColor = '#14171d';
            textColor = '#4b5563';
          } else {
            fillColor = NUC_COLORS[uChar] || '#ff5722';
            textColor = '#ffffff';
          }
        }

        // Draw Cell Rect
        if (charWidth >= 4) {
          ctx.fillStyle = fillColor;
          ctx.fillRect(colX, rowY, charWidth, charHeight);

          // Grid line
          ctx.strokeStyle = '#1b1f27';
          ctx.lineWidth = 0.5;
          ctx.strokeRect(colX, rowY, charWidth, charHeight);
        } else {
          // Micro bird's eye mode
          ctx.fillStyle = fillColor;
          ctx.fillRect(colX, rowY, Math.max(1, charWidth), charHeight);
        }

        // Draw Character Glyph if zoomed in enough
        if (charWidth >= 9 && charHeight >= 12) {
          ctx.fillStyle = textColor;
          ctx.font = `${Math.min(charHeight - 4, 12)}px monospace`;
          ctx.textAlign = 'center';
          ctx.textBaseline = 'middle';
          const glyph =
            !aminoAcidMode && colorScheme === 'difference' && uChar === consChar ? '·' : uChar;
          ctx.fillText(glyph, colX + charWidth / 2, rowY + charHeight / 2 + 0.5);
        }

        // Stop Codon Indicator: Solid black box with asterisk ONLY on middle base (no letters on outer bases)
        const stopCodon = stopCodonsList.find(
          (sc) => sc.taxon === taxonName && c >= sc.start && c < sc.end
        );
        if (stopCodon) {
          ctx.fillStyle = '#000000';
          ctx.fillRect(colX, rowY, charWidth, charHeight);

          ctx.strokeStyle = '#ffffff';
          ctx.lineWidth = 0.75;
          ctx.strokeRect(colX + 0.5, rowY + 0.5, charWidth - 1, charHeight - 1);

          // Asterisk ONLY on the middle base of the 3-base stop codon, no letters on outer bases
          if (c === stopCodon.start + 1) {
            ctx.fillStyle = '#ffffff';
            ctx.font = `bold ${Math.min(charHeight, 15)}px monospace`;
            ctx.textAlign = 'center';
            ctx.textBaseline = 'middle';
            ctx.fillText('*', colX + charWidth / 2, rowY + charHeight / 2 + 1);
          }
        }

        // Live Diff Overlay: Ghost / Slate-Grey trimmed columns
        if (isTrimmedCol) {
          ctx.fillStyle = activeOverlay.fill;
          ctx.fillRect(colX, rowY, charWidth, charHeight);

          // Subtle diagonal hatch
          ctx.strokeStyle = 'rgba(255, 255, 255, 0.15)';
          ctx.lineWidth = 0.75;
          ctx.beginPath();
          ctx.moveTo(colX, rowY + charHeight);
          ctx.lineTo(colX + charWidth, rowY);
          ctx.stroke();
        }

        // HMM / Segment Masked Overlay: Purple/Indigo diagonal hatch
        const isMaskedSeg =
          effectiveViewMode === 'overlays' &&
          showDiffOverlay &&
          diff.masked_segments.some(
            (seg) => seg.taxon === taxonName && c >= seg.start && c < seg.end
          );

        if (isMaskedSeg && !isTrimmedCol) {
          ctx.fillStyle = activeOverlay.fill;
          ctx.fillRect(colX, rowY, charWidth, charHeight);

          ctx.strokeStyle = 'rgba(255, 255, 255, 0.15)';
          ctx.lineWidth = 0.75;
          ctx.beginPath();
          ctx.moveTo(colX, rowY);
          ctx.lineTo(colX + charWidth, rowY + charHeight);
          ctx.stroke();
        }
      }

      // Dropped Taxon Ghost Mask (in overlays view mode)
      if (effectiveViewMode === 'overlays' && showDiffOverlay && isDropped) {
        ctx.fillStyle = activeOverlay.fill;
        ctx.fillRect(taxaNameWidth, rowY, width - taxaNameWidth, charHeight);
      }
    }

    // 2. Draw Left Taxa Name Column (Pinned Sticky Left with Clean Divider)
    ctx.fillStyle = '#14171d';
    ctx.fillRect(0, 0, taxaNameWidth, height);

    // Left Border & Shadow Line
    ctx.strokeStyle = '#2d3545';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.moveTo(taxaNameWidth, 0);
    ctx.lineTo(taxaNameWidth, height);
    ctx.stroke();

    for (let r = startRow; r < endRow; r++) {
      const taxonName = activeTaxa[r];
      const isDropped = effectiveViewMode === 'overlays' && droppedTaxaSet.has(taxonName);
      const rowY = topHeaderHeight + (r * charHeight - scrollY);

      // Alternating row background
      ctx.fillStyle = r % 2 === 0 ? '#171b22' : '#14171d';
      ctx.fillRect(0, rowY, taxaNameWidth - 1, charHeight);

      // Row separator
      ctx.strokeStyle = '#1f242e';
      ctx.lineWidth = 0.5;
      ctx.strokeRect(0, rowY, taxaNameWidth - 1, charHeight);

      ctx.font = '11px monospace';
      ctx.textAlign = 'left';
      ctx.textBaseline = 'middle';

      if (effectiveViewMode === 'overlays' && showDiffOverlay && isDropped) {
        ctx.fillStyle = '#94a3b8';
        const label = `✖ ${taxonName}`;
        ctx.fillText(label, 10, rowY + charHeight / 2);
        const textWidth = ctx.measureText(label).width;
        ctx.strokeStyle = '#ef4444';
        ctx.lineWidth = 1.2;
        ctx.beginPath();
        ctx.moveTo(10, rowY + charHeight / 2);
        ctx.lineTo(10 + Math.min(textWidth, taxaNameWidth - 20), rowY + charHeight / 2);
        ctx.stroke();
      } else {
        ctx.fillStyle = '#dce6ff';
        ctx.fillText(taxonName, 10, rowY + charHeight / 2);
      }
    }

    // 3. Draw Sticky Top Ruler & PIS Track
    ctx.fillStyle = '#171b22';
    ctx.fillRect(0, 0, width, topHeaderHeight);

    // Top-Left Corner Box
    ctx.fillStyle = '#14171d';
    ctx.fillRect(0, 0, taxaNameWidth, topHeaderHeight);
    ctx.strokeStyle = '#2d3545';
    ctx.lineWidth = 1.5;
    ctx.strokeRect(0, 0, taxaNameWidth, topHeaderHeight);

    // Taxa / Site Corner Labels and alignment dimensions. Keeping these in the
    // pinned canvas corner leaves the toolbar entirely available for controls.
    ctx.fillStyle = '#8b949e';
    ctx.font = '10px font-mono, monospace';
    ctx.textAlign = 'left';
    ctx.fillText('TAXA / SITES', 10, 16);
    ctx.fillText('CONSENSUS', 10, 34);

    ctx.font = '9px font-mono, monospace';
    ctx.textAlign = 'right';
    ctx.fillText(
      `Original ${raw_alignment.num_taxa} × ${raw_alignment.length.toLocaleString()} bp`,
      taxaNameWidth - 10,
      16
    );
    ctx.fillText(
      `Trimmed ${trimmed_alignment.num_taxa} × ${(
        aminoAcidMode ? length : trimmed_alignment.length
      ).toLocaleString()} ${aminoAcidMode ? 'aa' : 'bp'}`,
      taxaNameWidth - 10,
      34
    );

    // Ruler Ticks and PIS Indicators
    for (let c = startCol; c < endCol; c++) {
      const colX = taxaNameWidth + (c * charWidth - scrollX);
      const isTrimmedCol =
        effectiveViewMode === 'overlays' && showDiffOverlay && trimmedColsSet.has(c);
      const isPis = !aminoAcidMode && effectiveViewMode === 'overlays' && pis_mask[c];

      // Protein positions are more compact, so label every 10 residues.
      const rulerInterval = aminoAcidMode ? 10 : 20;
      if ((c + 1) % rulerInterval === 0 || c === 0) {
        ctx.fillStyle = '#8b949e';
        ctx.font = '9px monospace';
        ctx.textAlign = 'center';
        ctx.fillText(`${c + 1}`, colX + charWidth / 2, 13);
      }

      // PIS Star / Indicator
      if (isPis) {
        ctx.fillStyle = '#f59e0b';
        ctx.font = '10px sans-serif';
        ctx.textAlign = 'center';
        ctx.fillText('★', colX + charWidth / 2, 23);
      }

      // Majority Consensus Track
      const consChar = displayConsensus[c] || '-';
      ctx.fillStyle = aminoAcidMode
        ? activeAminoAcidPalette.colors[consChar] || '#8b949e'
        : NUC_COLORS[consChar] || '#8b949e';
      ctx.font = '10px monospace';
      ctx.textAlign = 'center';
      ctx.fillText(consChar, colX + charWidth / 2, 35);

      // Top Ruler Overlay if column is trimmed
      if (isTrimmedCol) {
        ctx.fillStyle = activeOverlay.fill;
        ctx.fillRect(colX, 0, charWidth, topHeaderHeight);
      }
    }
  }, [
    raw_alignment,
    trimmed_alignment,
    diff,
    pis_mask,
    displayConsensus,
    scrollX,
    scrollY,
    charWidth,
    charHeight,
    colorScheme,
    showDiffOverlay,
    effectiveViewMode,
    overlayColor,
    numTaxa,
    length,
    activeTaxa,
    activeSeqs,
    taxaNameWidth,
    stopCodonsList,
    aminoAcidMode,
    aminoAcidViewerSettings.colorScheme,
    aminoAcidViewerSettings.dimConsensusMatches,
  ]);

  // Mouse Wheel Pan / Scroll
  const handleWheel = (e: React.WheelEvent) => {
    e.preventDefault();
    if (e.ctrlKey || e.metaKey) {
      if (e.deltaY < 0) {
        setCharWidth((w) => Math.min(32, w + 1));
      } else {
        setCharWidth((w) => Math.max(1, w - 1));
      }
    } else {
      setScrollX((prev) =>
        Math.max(0, Math.min(length * charWidth, prev + e.deltaX))
      );
      setScrollY((prev) =>
        Math.max(0, Math.min(numTaxa * charHeight, prev + e.deltaY))
      );
    }
  };

  const handleMouseDownCanvas = (e: React.MouseEvent<HTMLCanvasElement>) => {
    setIsDragging(true);
    lastMousePosRef.current = { x: e.clientX, y: e.clientY };
  };

  const handleMouseUpCanvas = () => {
    setIsDragging(false);
    lastMousePosRef.current = null;
  };

  // Mouse Move for Cell Hover Info and Panning
  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (isDragging && lastMousePosRef.current) {
      const dx = e.clientX - lastMousePosRef.current.x;
      const dy = e.clientY - lastMousePosRef.current.y;

      setScrollX((prev) =>
        Math.max(0, Math.min(length * charWidth, prev - dx))
      );
      setScrollY((prev) =>
        Math.max(0, Math.min(numTaxa * charHeight, prev - dy))
      );

      lastMousePosRef.current = { x: e.clientX, y: e.clientY };
      setHoverInfo(null);
      return;
    }

    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;

    if (mouseX <= taxaNameWidth || mouseY <= topHeaderHeight) {
      setHoverInfo(null);
      return;
    }

    const col = Math.floor((mouseX - taxaNameWidth + scrollX) / charWidth);
    const row = Math.floor((mouseY - topHeaderHeight + scrollY) / charHeight);

    if (col >= 0 && col < length && row >= 0 && row < numTaxa) {
      const taxon = activeTaxa[row];
      const base = activeSeqs[row]?.[col] || '-';
      const consensus = displayConsensus[col] || '-';
      const isPis =
        !aminoAcidMode && effectiveViewMode === 'overlays' && (pis_mask[col] || false);
      const isTrimmed =
        effectiveViewMode === 'overlays' && diff.trimmed_columns.includes(col);
      const isDroppedTaxon =
        effectiveViewMode === 'overlays' && diff.dropped_taxa.includes(taxon);

      // Trimming reason
      let trimReason: string | null = null;
      if (isTrimmed) {
        if (diff.column_reasons && diff.column_reasons[col]) {
          trimReason = diff.column_reasons[col];
        } else {
          trimReason = 'Trimmed by active quality filter';
        }
      }

      // Dropped reason
      let droppedReason: string | null = null;
      if (isDroppedTaxon) {
        if (diff.dropped_taxa_reasons && diff.dropped_taxa_reasons[taxon]) {
          droppedReason = diff.dropped_taxa_reasons[taxon];
        } else {
          droppedReason = 'Sample pruned due to quality criteria';
        }
      }

      // Stop codon check
      const stopCodon = stopCodonsList.find(
        (sc) => sc.taxon === taxon && col >= sc.start && col < sc.end
      );

      setHoverInfo({
        taxon,
        col: col + 1,
        base,
        consensus,
        isPis,
        isTrimmed,
        trimReason,
        isDroppedTaxon,
        droppedReason,
        isStopCodon: !!stopCodon,
        stopCodonSeq: stopCodon ? stopCodon.codon : null,
        screenX: e.clientX,
        screenY: e.clientY,
      });
    } else {
      setHoverInfo(null);
    }
  };

  const handleMouseLeave = () => {
    setHoverInfo(null);
    setIsDragging(false);
    lastMousePosRef.current = null;
  };

  return (
    <div className="flex-1 flex flex-col h-full bg-[#0e1014] overflow-hidden select-none relative">
      {/* Top Toolbar */}
      <div className="relative h-11 px-4 py-1 bg-[#14171d] border-b border-[#232833] flex items-center flex-none gap-4">
        {/* Locus Metadata & Status */}
        <div className="flex max-w-[340px] shrink-0 items-center gap-3 overflow-hidden text-xs">
          <span
            className="max-w-[180px] shrink-0 truncate font-mono text-sm font-semibold text-[#dce6ff]"
            title={raw_alignment.id}
          >
            {raw_alignment.id}
          </span>

          {/* This is always nucleotide alignment QC. Protein display never re-gates a locus. */}
          {diff.pass ? (
            <span className="inline-flex items-center gap-1 text-[11px] text-emerald-400 font-semibold bg-emerald-500/10 px-2 py-0.5 rounded border border-emerald-500/20 shrink-0">
              <CheckCircle className="w-3.5 h-3.5" />
              <span>Alignment QC: Pass</span>
            </span>
          ) : (
            <span
              className="inline-flex items-center gap-1 text-[11px] text-rose-400 font-semibold bg-rose-500/10 px-2 py-0.5 rounded border border-rose-500/20 cursor-help shrink-0"
              title={diff.fail_reasons?.join(', ') || 'Quality gate failed'}
            >
              <XCircle className="w-3.5 h-3.5 shrink-0" />
              <span>Alignment QC: Fail</span>
            </span>
          )}

        </div>

        {/* Display Mode & Knobs */}
        <div
          className="no-scrollbar min-w-0 flex-1 overflow-x-auto overflow-y-hidden"
          onScroll={() => setShowColorMenu(false)}
        >
        <div className="ml-auto flex w-max min-w-[792px] items-center justify-end gap-2.5 whitespace-nowrap">
          <button
            type="button"
            disabled={!aminoAcidAvailable && !aminoAcidViewerSettings.enabled}
            onClick={() =>
              onChangeAminoAcidViewerSettings({
                ...aminoAcidViewerSettings,
                enabled: !aminoAcidViewerSettings.enabled,
              })
            }
            className={`flex w-[126px] shrink-0 items-center justify-center gap-1.5 whitespace-nowrap rounded-md border px-2 py-1 text-xs font-medium transition-colors disabled:cursor-not-allowed disabled:opacity-35 ${
              aminoAcidMode
                ? 'border-violet-400/50 bg-violet-500/25 text-violet-100'
                : 'border-[#2d3545] bg-[#171b22] text-[#8b949e] hover:text-violet-300'
            }`}
            title={
              aminoAcidAvailable
                ? 'Switch between nucleotide and translated amino-acid display'
                : 'Amino-acid display requires an accepted ORF'
            }
          >
            <Languages className="h-3.5 w-3.5" />
            <span>{aminoAcidMode ? 'Amino Acids' : 'Nucleotides'}</span>
          </button>

          {/* Alignment View Mode Toggle: Overlays vs Final Trimmed */}
          <div className="grid w-[262px] shrink-0 grid-cols-[144px_112px] items-center p-0.5 bg-[#0e1014] rounded-lg border border-[#232833]">
            <button
              onClick={() => setViewMode('overlays')}
              disabled={aminoAcidMode}
              className={`flex shrink-0 items-center justify-center gap-1.5 whitespace-nowrap px-2 py-1 rounded-md text-xs font-medium transition-all ${
                effectiveViewMode === 'overlays'
                  ? 'bg-blue-600 text-white shadow-sm font-semibold'
                  : 'text-[#8b949e] hover:text-[#c9d1d9] disabled:cursor-not-allowed disabled:opacity-35'
              }`}
              title={
                aminoAcidMode
                  ? 'Raw overlays are nucleotide coordinates and are unavailable in amino-acid view'
                  : 'View raw alignment with live grey trimming overlays'
              }
            >
              <Eye className="w-3.5 h-3.5" />
              <span>Trimming Overlays</span>
            </button>

            <button
              onClick={() => setViewMode('final')}
              className={`flex shrink-0 items-center justify-center gap-1.5 whitespace-nowrap px-2 py-1 rounded-md text-xs font-medium transition-all ${
                effectiveViewMode === 'final'
                  ? 'bg-emerald-600 text-white shadow-sm font-semibold'
                  : 'text-[#8b949e] hover:text-[#c9d1d9]'
              }`}
              title="View final post-trimmed alignment without deleted columns/taxa"
            >
              <Sparkles className="w-3.5 h-3.5" />
              <span>Final Trimmed</span>
            </button>
          </div>

          {/* Keep this slot mounted so the controls never shift between views. */}
          <div
            className={`relative w-[110px] shrink-0 ${
              effectiveViewMode === 'overlays' ? '' : 'invisible pointer-events-none'
            }`}
          >
              <button
                onClick={(event) => {
                  if (!showColorMenu) {
                    const rect = event.currentTarget.getBoundingClientRect();
                    setOverlayMenuPosition({
                      left: Math.max(8, Math.min(window.innerWidth - 200, rect.right - 192)),
                      top: rect.bottom + 6,
                    });
                  }
                  setShowColorMenu((visible) => !visible);
                }}
                className="flex w-full items-center justify-center gap-1.5 whitespace-nowrap px-2 py-1 rounded-md bg-[#171b22] hover:bg-[#1f242e] text-[#c9d1d9] border border-[#232833] text-xs font-medium transition-colors"
                title="Change Trimming Overlay Color"
              >
                <div
                  className="w-3 h-3 rounded-full border border-white/30"
                  style={{ backgroundColor: OVERLAY_COLORS[overlayColor].stroke }}
                />
                <span>Overlay Color</span>
              </button>
          </div>

          {/* Color Scheme Selector */}
          <div className="flex w-[176px] shrink-0 items-center gap-1 overflow-hidden bg-[#171b22] px-2 py-1 rounded-md border border-[#232833]">
            <Palette className={`w-3.5 h-3.5 ${aminoAcidMode ? 'text-violet-400' : 'text-blue-400'}`} />
            {aminoAcidMode ? (
              <select
                value={aminoAcidViewerSettings.colorScheme ?? 'chemistry'}
                onChange={(event) =>
                  onChangeAminoAcidViewerSettings({
                    ...aminoAcidViewerSettings,
                    colorScheme: event.target.value as AminoAcidColorScheme,
                  })
                }
                className="min-w-0 flex-1 bg-transparent text-xs text-[#c9d1d9] outline-none cursor-pointer"
              >
                <option value="chemistry" className="bg-[#171b22]">Chemistry</option>
                <option value="colorblind" className="bg-[#171b22]">Color-blind Safe</option>
                <option value="zappo" className="bg-[#171b22]">Zappo</option>
                <option value="taylor" className="bg-[#171b22]">Taylor</option>
                <option value="hydrophobicity" className="bg-[#171b22]">Hydrophobicity</option>
                <option value="monochrome" className="bg-[#171b22]">Monochrome</option>
              </select>
            ) : (
              <select
                value={colorScheme}
                onChange={(e) => onChangeColorScheme(e.target.value as ColorScheme)}
                className="min-w-0 flex-1 bg-transparent text-xs text-[#c9d1d9] outline-none cursor-pointer"
              >
                <option value="nucleotide" className="bg-[#171b22]">Classic</option>
                <option value="nucleotide-colorblind" className="bg-[#171b22]">Color-blind Safe</option>
                <option value="nucleotide-high-contrast" className="bg-[#171b22]">High Contrast</option>
                <option value="nucleotide-pastel" className="bg-[#171b22]">Pastel</option>
                <option value="nucleotide-ocean" className="bg-[#171b22]">Ocean</option>
                <option value="nucleotide-sunset" className="bg-[#171b22]">Sunset</option>
                <option value="nucleotide-forest" className="bg-[#171b22]">Forest</option>
                <option value="clustalx" className="bg-[#171b22]">ClustalX</option>
                <option value="difference" className="bg-[#171b22]">Differences</option>
                <option value="monochrome" className="bg-[#171b22]">Monochrome</option>
              </select>
            )}
          </div>

          {/* Zoom controls */}
          <div className="flex w-[78px] shrink-0 items-center justify-center gap-1 bg-[#171b22] p-0.5 rounded-md border border-[#232833]">
            <button
              onClick={() => setCharWidth((w) => Math.max(1, w - 2))}
              className="p-1 text-[#8b949e] hover:text-[#c9d1d9] rounded hover:bg-[#232833]"
              title="Zoom Out"
            >
              <ZoomOut className="w-3.5 h-3.5" />
            </button>
            <span className="text-[11px] font-mono text-[#8b949e] px-1">{charWidth}px</span>
            <button
              onClick={() => setCharWidth((w) => Math.min(32, w + 2))}
              className="p-1 text-[#8b949e] hover:text-[#c9d1d9] rounded hover:bg-[#232833]"
              title="Zoom In"
            >
              <ZoomIn className="w-3.5 h-3.5" />
            </button>
          </div>
        </div>
        </div>

        {showColorMenu && effectiveViewMode === 'overlays' && (
          <div
            className="fixed z-50 w-48 space-y-1 rounded-lg border border-[#2d3545] bg-[#171b22] p-1.5 text-xs shadow-2xl"
            style={overlayMenuPosition}
          >
            {(Object.keys(OVERLAY_COLORS) as OverlayColorOption[]).map((key) => {
              const opt = OVERLAY_COLORS[key];
              const isSelected = overlayColor === key;
              return (
                <button
                  key={key}
                  onClick={() => {
                    setOverlayColor(key);
                    setShowColorMenu(false);
                  }}
                  className={`flex w-full items-center justify-between rounded px-2.5 py-1.5 text-left transition-colors ${
                    isSelected
                      ? 'bg-blue-600/20 font-medium text-[#dce6ff]'
                      : 'text-[#8b949e] hover:bg-[#1f242e] hover:text-[#c9d1d9]'
                  }`}
                >
                  <div className="flex items-center gap-2">
                    <span
                      className="h-3 w-3 rounded-full border border-white/20"
                      style={{ backgroundColor: opt.stroke }}
                    />
                    <span>{opt.label}</span>
                  </div>
                  {isSelected && <Check className="h-3.5 w-3.5 text-blue-400" />}
                </button>
              );
            })}
          </div>
        )}
      </div>

      {(diff.orf_evaluated ||
        diff.orf_reference_evaluated ||
        !diff.pass ||
        (aminoAcidViewerSettings.enabled && !aminoAcidAvailable)) && (
        <div className="no-scrollbar flex-none min-h-8 px-4 py-1.5 border-b border-[#232833] bg-[#11151b] flex items-center gap-4 whitespace-nowrap text-[10px] font-mono text-[#8b949e] overflow-x-auto">
          {(diff.orf_evaluated || diff.orf_reference_evaluated) && (
            <>
              <span className="text-emerald-400">
                Trimmed: {diff.new_taxa_count} taxa × {diff.new_length.toLocaleString()} bp ({diff.new_gap_percent.toFixed(1)}% gap)
              </span>
              {diff.dropped_taxa.length > 0 && (
                <span className="rounded border border-rose-500/20 bg-rose-500/10 px-1.5 py-0.5 text-rose-400">
                  -{diff.dropped_taxa.length} taxa
                </span>
              )}
              {diff.trimmed_columns.length > 0 && (
                <span className="rounded border border-amber-500/20 bg-amber-500/10 px-1.5 py-0.5 text-amber-400">
                  -{diff.trimmed_columns.length} cols
                </span>
              )}
              <span
                className={
                  diff.orf_candidate_found && diff.found_valid_orf
                    ? 'text-emerald-400 font-semibold'
                    : diff.orf_candidate_found
                      ? 'text-amber-400 font-semibold'
                      : 'text-rose-400 font-semibold'
                }
              >
                ORF candidate: {diff.orf_candidate_found ? (diff.found_valid_orf ? 'accepted' : 'rejected') : 'none'}
              </span>
              <span>Frame: {diff.orf_frame ? `${diff.orf_frame > 0 ? '+' : ''}${diff.orf_frame}` : '—'}</span>
              <span>Support: {diff.orf_support_count ?? 0} ({Number(diff.orf_support_percent ?? 0).toFixed(0)}%)</span>
              <span title={`AA conservation ${Number(diff.orf_amino_acid_conservation ?? 0).toFixed(1)}%; frame contrast ${Number(diff.orf_frame_contrast ?? 0).toFixed(1)}%`}>
                Coding evidence: {Number(diff.orf_coding_score ?? 0).toFixed(1)}/100
              </span>
              {diff.orf_reference_evaluated && (
                <span className={diff.orf_reference_matched ? 'text-blue-300' : 'text-rose-300'}>
                  Reference: {diff.orf_reference_matched
                    ? `matched (${Number(diff.orf_reference_identity ?? 0).toFixed(1)}% identity, ${Number(diff.orf_reference_coverage ?? 0).toFixed(1)}% coverage)`
                    : 'no match'}
                </span>
              )}
              {diff.orf_reference_matched && <span>Separated intron: {diff.orf_intron_length ?? 0} bp</span>}
            </>
          )}

          <span className="ml-auto flex shrink-0 items-center gap-3">
            {aminoAcidViewerSettings.enabled && !aminoAcidAvailable && (
              <span
                className="rounded border border-amber-500/20 bg-amber-500/10 px-2 py-0.5 font-sans font-medium text-amber-300"
                title="This locus must have an assessed, accepted ORF before it can be translated safely."
              >
                Amino-acid view unavailable for this locus
              </span>
            )}
            {!diff.pass && (
              <span
                className="rounded border border-rose-500/20 bg-rose-500/10 px-2 py-0.5 font-sans font-medium text-rose-300"
                title={diff.fail_reasons?.join(', ') || 'Quality gate failed'}
              >
                Fail reason: {diff.fail_reasons?.[0] || 'Quality gate failed'}
                {(diff.fail_reasons?.length ?? 0) > 1
                  ? ` (+${diff.fail_reasons.length - 1} more)`
                  : ''}
              </span>
            )}
          </span>
        </div>
      )}

      {/* Main Canvas Viewport */}
      <div
        ref={containerRef}
        onWheel={handleWheel}
        className={`flex-1 w-full h-full relative overflow-hidden bg-[#0e1014] ${
          isDragging ? 'cursor-grabbing' : 'cursor-crosshair'
        }`}
      >
        <canvas
          ref={canvasRef}
          onMouseMove={handleMouseMove}
          onMouseDown={handleMouseDownCanvas}
          onMouseUp={handleMouseUpCanvas}
          onMouseLeave={handleMouseLeave}
          className="w-full h-full block"
        />

        {/* Draggable Divider Handle between Taxon Names and Sequence Grid */}
        <div
          onMouseDown={handleDividerMouseDown}
          onDoubleClick={() => setCustomTaxaWidth(null)}
          className={`absolute top-0 bottom-0 w-3 -ml-1.5 cursor-col-resize z-30 group flex items-center justify-center hover:bg-blue-500/20 transition-colors ${
            isResizingDivider ? 'bg-blue-500/30' : ''
          }`}
          style={{ left: `${taxaNameWidth}px` }}
          title="Drag to resize taxon names column (Double-click to reset to auto-fit)"
        >
          <div className="w-0.5 h-full bg-[#2d3545] group-hover:bg-blue-400 transition-colors" />
        </div>

        {/* Hover Inspection Tooltip */}
        {hoverInfo && (
          <div
            className="fixed pointer-events-none z-50 bg-[#171b22] border border-[#2d3545] rounded-lg p-2.5 shadow-2xl text-xs text-[#c9d1d9] space-y-1 backdrop-blur-sm max-w-sm"
            style={{
              left: `${hoverInfo.screenX + 14}px`,
              top: `${hoverInfo.screenY + 14}px`,
            }}
          >
            <div className="font-semibold text-[#dce6ff] font-mono">{hoverInfo.taxon}</div>
            <div className="text-[11px] text-[#8b949e]">
              {aminoAcidMode ? 'AA' : 'Pos'}: <b className="text-[#dce6ff]">{hoverInfo.col}</b> | {aminoAcidMode ? 'Residue' : 'Base'}: <b className="text-cyan-400">{hoverInfo.base}</b> | Consensus: <b>{hoverInfo.consensus}</b>
            </div>
            {hoverInfo.isPis && (
              <div className="text-[10px] text-amber-400 font-medium">★ Parsimony-Informative Site</div>
            )}
            {hoverInfo.isStopCodon && (
              <div className="text-[10px] text-white bg-black px-1.5 py-0.5 rounded border border-white font-semibold flex items-center gap-1">
                <span className="font-mono font-bold">*</span> Stop Codon: {hoverInfo.stopCodonSeq} (Translates to STOP)
              </div>
            )}
            {hoverInfo.isTrimmed && effectiveViewMode === 'overlays' && (
              <div className="mt-1 pt-1 border-t border-slate-700/60 flex flex-col gap-0.5">
                <div className="text-amber-400 font-semibold flex items-center gap-1 text-[11px]">
                  <span>✂</span> Slated for Trimming:
                </div>
                <div className="text-slate-300 text-[10px] pl-3 font-mono">
                  {hoverInfo.trimReason || 'Trimmed by active quality filter'}
                </div>
              </div>
            )}
            {hoverInfo.isDroppedTaxon && effectiveViewMode === 'overlays' && (
              <div className="mt-1 pt-1 border-t border-rose-500/30 flex flex-col gap-0.5">
                <div className="text-rose-400 font-semibold flex items-center gap-1 text-[11px]">
                  <span>✖</span> Sample Pruned:
                </div>
                <div className="text-rose-300 text-[10px] pl-3 font-mono">
                  {hoverInfo.droppedReason || 'Low coverage / divergent outlier'}
                </div>
              </div>
            )}
          </div>
        )}
      </div>

      {/* Bottom Trim Status Inspector Bar */}
      <div className="h-7 px-4 bg-[#14171d] border-t border-[#232833] flex items-center justify-between text-xs text-[#8b949e] font-mono flex-none">
        <div className="flex items-center gap-4 truncate">
          {hoverInfo ? (
            <>
              <span><strong className="text-[#dce6ff]">{aminoAcidMode ? 'AA' : 'Pos'}:</strong> {hoverInfo.col}</span>
              <span className="truncate max-w-[200px]"><strong className="text-[#dce6ff]">Taxon:</strong> {hoverInfo.taxon}</span>
              <span><strong className="text-[#dce6ff]">State:</strong> <span className="font-bold text-white px-1 rounded bg-[#1f242e]">{hoverInfo.base}</span></span>
              <span><strong className="text-[#dce6ff]">Consensus:</strong> {hoverInfo.consensus}</span>
              {hoverInfo.isPis && <span className="text-amber-400 font-bold">★ PIS</span>}
              {hoverInfo.isStopCodon && (
                <span className="text-white bg-black px-1.5 py-0.5 rounded border border-white text-[10px] font-bold font-mono">
                  * STOP CODON ({hoverInfo.stopCodonSeq})
                </span>
              )}
              {hoverInfo.isTrimmed && effectiveViewMode === 'overlays' && (
                <span className="text-amber-300 bg-amber-500/10 px-2 py-0.5 rounded border border-amber-500/20 font-sans text-[11px] flex items-center gap-1">
                  <span>✂</span> <strong>Why Trimming:</strong> {hoverInfo.trimReason}
                </span>
              )}
            </>
          ) : (
            <span className="text-[#4b5563] italic flex items-center gap-1.5 font-sans text-[11px]">
              <Info className="w-3.5 h-3.5 text-blue-400" />
              {aminoAcidMode
                ? 'Hover over any residue to inspect the translated state and protein consensus'
                : 'Hover over any position or column to inspect nucleotide stats and trimming rationale'}
            </span>
          )}
        </div>
        <div className="flex items-center gap-3 text-[10px] text-[#6b7280]">
          <span>Scroll: Pan</span>
          <span>Ctrl+Scroll: Zoom</span>
        </div>
      </div>
    </div>
  );
};

export default MsaViewer;
