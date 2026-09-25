import React, { useRef, useEffect, useMemo, useState } from 'react';
import { TaxonOccupancy, AlignmentSummary } from '../types';
import { Grid } from 'lucide-react';

interface MatrixHeatmapProps {
  occupancy: TaxonOccupancy[];
  summaries: AlignmentSummary[];
  onSelectLocus: (id: string, filePath: string) => void;
}

/** Missing-data percent of one sample in a processed locus, or null when the sample is absent. */
function sampleMissingPercent(summary: AlignmentSummary, taxon: string): number | null {
  const basepairs = summary.retained_taxon_basepairs?.[taxon];
  if (basepairs === undefined) return null;
  if (summary.length <= 0) return 100;
  return 100 * (1 - basepairs / summary.length);
}

export const MatrixHeatmap: React.FC<MatrixHeatmapProps> = ({
  occupancy,
  summaries,
  onSelectLocus,
}) => {
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const containerRef = useRef<HTMLDivElement | null>(null);

  const [scrollX, setScrollX] = useState(0);
  const [scrollY, setScrollY] = useState(0);
  const [canvasSize, setCanvasSize] = useState({ width: 0, height: 0 });
  const [hoverInfo, setHoverInfo] = useState<{
    taxon: string;
    locusId: string;
    locusLength: number;
    missingPercent: number | null;
    screenX: number;
    screenY: number;
  } | null>(null);

  const numTaxa = occupancy.length;
  const numLoci = summaries.length;
  const hasData = numTaxa > 0 && numLoci > 0;
  // The browser build fetches the per-sample lists for this view, so they can
  // arrive after the summaries. A locus with samples but no list is not loaded.
  const sampleDataReady = useMemo(
    () =>
      summaries.every(
        (summary) => summary.num_taxa === 0 || (summary.retained_taxa?.length ?? 0) > 0
      ),
    [summaries]
  );

  const cellWidth = 4; // px per locus column
  const cellHeight = 16; // px per taxon row
  const leftTaxaWidth = 160;

  // Redraw when the window or the layout changes the canvas size.
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const observer = new ResizeObserver(() => {
      setCanvasSize({ width: canvas.clientWidth, height: canvas.clientHeight });
    });
    observer.observe(canvas);
    return () => observer.disconnect();
  }, [hasData]);

  // Render Canvas
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

    // Background
    ctx.fillStyle = '#0e1014';
    ctx.fillRect(0, 0, width, height);

    if (!hasData) return;

    const startCol = Math.max(0, Math.floor(scrollX / cellWidth));
    const endCol = Math.min(numLoci, startCol + Math.ceil((width - leftTaxaWidth) / cellWidth) + 1);

    const startRow = Math.max(0, Math.floor(scrollY / cellHeight));
    const endRow = Math.min(numTaxa, startRow + Math.ceil(height / cellHeight) + 1);

    // Draw Heatmap Cells: each sample's missing data in each processed locus
    if (sampleDataReady) {
      for (let r = startRow; r < endRow; r++) {
        const y = r * cellHeight - scrollY;
        const taxon = occupancy[r].taxon_name;

        for (let c = startCol; c < endCol; c++) {
          const x = leftTaxaWidth + c * cellWidth - scrollX;
          const missingPercent = sampleMissingPercent(summaries[c], taxon);

          if (missingPercent === null) {
            ctx.fillStyle = '#1f242e'; // absent
          } else if (missingPercent < 20) {
            ctx.fillStyle = '#10b981'; // green
          } else if (missingPercent < 50) {
            ctx.fillStyle = '#f59e0b'; // amber
          } else {
            ctx.fillStyle = '#ef4444'; // rose
          }

          ctx.fillRect(x, y, cellWidth - 0.5, cellHeight - 1);
        }
      }
    }

    // Draw Left Sticky Taxon Labels Background
    ctx.fillStyle = '#14171d';
    ctx.fillRect(0, 0, leftTaxaWidth, height);
    ctx.fillStyle = '#232833';
    ctx.fillRect(leftTaxaWidth - 1, 0, 1, height);

    // Draw Sticky Taxon Labels
    ctx.font = '10px monospace';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';

    for (let r = startRow; r < endRow; r++) {
      const y = r * cellHeight - scrollY;
      const taxon = occupancy[r];

      ctx.fillStyle = '#8b949e';
      ctx.fillText(
        taxon.taxon_name.length > 18
          ? taxon.taxon_name.slice(0, 17) + '…'
          : taxon.taxon_name,
        8,
        y + cellHeight / 2
      );
    }
  }, [occupancy, summaries, scrollX, scrollY, numTaxa, numLoci, hasData, sampleDataReady, canvasSize]);

  // React registers wheel listeners as passive, where preventDefault has no
  // effect, so the matrix scroll uses a native listener.
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;
    const handleWheel = (e: WheelEvent) => {
      e.preventDefault();
      const maxScrollX = Math.max(0, numLoci * cellWidth - container.clientWidth + leftTaxaWidth);
      const maxScrollY = Math.max(0, numTaxa * cellHeight - container.clientHeight);

      setScrollX((prev) => Math.max(0, Math.min(maxScrollX, prev + e.deltaX)));
      setScrollY((prev) => Math.max(0, Math.min(maxScrollY, prev + e.deltaY)));
    };
    container.addEventListener('wheel', handleWheel, { passive: false });
    return () => container.removeEventListener('wheel', handleWheel);
  }, [numLoci, numTaxa]);

  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const mouseX = e.clientX - rect.left;
    const mouseY = e.clientY - rect.top;

    if (mouseX <= leftTaxaWidth) {
      setHoverInfo(null);
      return;
    }

    const col = Math.floor((mouseX - leftTaxaWidth + scrollX) / cellWidth);
    const row = Math.floor((mouseY + scrollY) / cellHeight);

    if (sampleDataReady && col >= 0 && col < numLoci && row >= 0 && row < numTaxa) {
      const taxon = occupancy[row];
      const locus = summaries[col];
      setHoverInfo({
        taxon: taxon.taxon_name,
        locusId: locus.id,
        locusLength: locus.length,
        missingPercent: sampleMissingPercent(locus, taxon.taxon_name),
        screenX: e.clientX,
        screenY: e.clientY,
      });
    } else {
      setHoverInfo(null);
    }
  };

  const handleClick = (e: React.MouseEvent<HTMLCanvasElement>) => {
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;

    const mouseX = e.clientX - rect.left;
    if (mouseX <= leftTaxaWidth) return;

    const col = Math.floor((mouseX - leftTaxaWidth + scrollX) / cellWidth);
    if (col >= 0 && col < numLoci) {
      const locus = summaries[col];
      onSelectLocus(locus.id, locus.file_path);
    }
  };

  if (occupancy.length === 0 || summaries.length === 0) {
    return (
      <div className="flex-1 flex items-center justify-center bg-[#0e1014] text-[#8b949e]">
        No dataset loaded
      </div>
    );
  }

  return (
    <div className="flex-1 flex flex-col h-full bg-[#0e1014] overflow-hidden select-none">
      {/* Top Header */}
      <div className="h-11 px-4 bg-[#14171d] border-b border-[#232833] flex items-center justify-between flex-none">
        <div className="flex items-center gap-2">
          <Grid className="w-4 h-4 text-cyan-400" />
          <span className="font-semibold text-xs text-[#dce6ff]">
            Taxa Occupancy & Matrix Completeness ({occupancy.length} unique taxa × {numLoci.toLocaleString()} loci)
          </span>
        </div>

        <div className="flex items-center gap-3 text-[11px] text-[#8b949e]">
          {!sampleDataReady && <span className="text-cyan-400">Loading per-sample data…</span>}
          <div className="flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-sm bg-emerald-500 inline-block" />
            <span>Present (&lt;20% missing)</span>
          </div>
          <div className="flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-sm bg-amber-500 inline-block" />
            <span>Present (20-50% missing)</span>
          </div>
          <div className="flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-sm bg-rose-500 inline-block" />
            <span>Present (&gt;50% missing)</span>
          </div>
          <div className="flex items-center gap-1.5">
            <span className="w-2.5 h-2.5 rounded-sm bg-[#1f242e] inline-block" />
            <span>Absent</span>
          </div>
        </div>
      </div>

      {/* Main Occupancy Split View */}
      <div className="flex-1 flex overflow-hidden">
        {/* Left: Taxa Occupancy Ranking Table (renders every row in a scrolling list) */}
        <div className="w-80 border-r border-[#232833] bg-[#14171d] flex flex-col overflow-hidden flex-none">
          <div className="p-2.5 bg-[#171b22] border-b border-[#232833] text-xs font-semibold text-[#8b949e] font-mono">
            TAXON OCCUPANCY RANKING
          </div>
          <div className="flex-1 overflow-y-auto">
            <table className="w-full text-left text-xs border-collapse font-mono">
              <thead className="bg-[#1b2029] text-[#8b949e] text-[10px] sticky top-0">
                <tr>
                  <th className="py-2 px-3">TAXON</th>
                  <th className="py-2 px-3 text-right">LOCI</th>
                  <th className="py-2 px-3 text-right">OCCUPANCY</th>
                </tr>
              </thead>
              <tbody className="divide-y divide-[#1b1f27]">
                {occupancy.map((t) => (
                  <tr key={t.taxon_name} className="hover:bg-[#171b22] text-[#c9d1d9]">
                    <td className="py-1.5 px-3 truncate max-w-[120px]" title={t.taxon_name}>
                      {t.taxon_name}
                    </td>
                    <td className="py-1.5 px-3 text-right text-cyan-400 font-medium">
                      {t.present_loci_count}
                    </td>
                    <td className="py-1.5 px-3 text-right">
                      <div className="flex items-center justify-end gap-1.5">
                        <div className="w-10 bg-[#1f242e] h-1.5 rounded-full overflow-hidden">
                          <div
                            className="bg-cyan-500 h-full rounded-full"
                            style={{ width: `${t.present_loci_percent}%` }}
                          />
                        </div>
                        <span className="text-[11px]">{t.present_loci_percent.toFixed(0)}%</span>
                      </div>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>

        {/* Right: High-Performance Canvas Matrix Heatmap */}
        <div
          ref={containerRef}
          className="flex-1 relative overflow-hidden bg-[#0e1014] cursor-pointer"
        >
          <canvas
            ref={canvasRef}
            onMouseMove={handleMouseMove}
            onMouseLeave={() => setHoverInfo(null)}
            onClick={handleClick}
            className="w-full h-full block"
          />

          {/* Hover Tooltip */}
          {hoverInfo && (
            <div
              className="fixed pointer-events-none z-50 bg-[#171b22] border border-[#2d3545] rounded-lg p-2.5 shadow-2xl text-xs text-[#c9d1d9] space-y-1 backdrop-blur-sm"
              style={{
                left: `${hoverInfo.screenX + 12}px`,
                top: `${hoverInfo.screenY + 12}px`,
              }}
            >
              <div className="font-semibold text-[#dce6ff] font-mono">{hoverInfo.locusId}</div>
              <div className="text-[11px] text-[#8b949e]">Taxon: {hoverInfo.taxon}</div>
              <div className="text-[11px] text-[#8b949e]">
                Length: {hoverInfo.locusLength} bp |{' '}
                {hoverInfo.missingPercent === null
                  ? 'Sample absent'
                  : `Missing: ${hoverInfo.missingPercent.toFixed(1)}%`}
              </div>
              <div className="text-[10px] text-cyan-400 pt-0.5">Click to inspect in MSA viewer</div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default MatrixHeatmap;
