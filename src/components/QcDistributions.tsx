import React from 'react';
import { AlignmentSummary, DatasetOverview, TaxonOccupancy } from '../types';
import {
  Activity,
  BarChart2,
  CheckCircle2,
  Dna,
  Download,
  Layers,
  Loader2,
  Minus,
  Ruler,
  Scissors,
  UsersRound,
} from 'lucide-react';

interface QcDistributionsProps {
  overview: DatasetOverview;
  summaries: AlignmentSummary[];
  occupancy: TaxonOccupancy[];
  orfEnabled: boolean;
  onExportStats: () => Promise<string | null>;
}

type OccupancyMetric = 'loci' | 'basepairs';
type OccupancyOverlay = 'none' | 'trimming' | 'orf';

const compactNumberFormatter = new Intl.NumberFormat('en', {
  notation: 'compact',
  maximumFractionDigits: 1,
});

export const QcDistributions: React.FC<QcDistributionsProps> = React.memo(({
  overview,
  summaries,
  occupancy,
  orfEnabled,
  onExportStats,
}) => {
  const [isExporting, setIsExporting] = React.useState(false);
  const [exportNotice, setExportNotice] = React.useState<string | null>(null);
  const [occupancyMetric, setOccupancyMetric] = React.useState<OccupancyMetric>('loci');
  const [occupancyOverlay, setOccupancyOverlay] = React.useState<OccupancyOverlay>('trimming');

  React.useEffect(() => {
    if (!orfEnabled && occupancyOverlay === 'orf') setOccupancyOverlay('trimming');
  }, [orfEnabled, occupancyOverlay]);

  const handleExportStats = async () => {
    setIsExporting(true);
    setExportNotice(null);
    try {
      const filePath = await onExportStats();
      if (filePath) {
        setExportNotice(`Saved ${filePath.split(/[\\/]/).pop() || 'CSV'}`);
      }
    } catch (error) {
      setExportNotice(error instanceof Error ? error.message : String(error));
    } finally {
      setIsExporting(false);
    }
  };

  // Helper to generate 10-bin histogram data
  const lengthBins = React.useMemo(() => {
    const values = summaries.map((s) => s.length);
    if (values.length === 0) return [];
    const min = Math.min(...values);
    const max = Math.max(...values);
    const numBins = 12;
    const binSize = (max - min) / numBins || 1;
    const bins = Array.from({ length: numBins }, (_, i) => ({
      start: min + i * binSize,
      end: min + (i + 1) * binSize,
      count: 0,
    }));
    for (const v of values) {
      const idx = Math.min(numBins - 1, Math.floor((v - min) / binSize));
      if (idx >= 0 && idx < numBins) bins[idx].count++;
    }
    const maxCount = Math.max(...bins.map((b) => b.count), 1);
    return bins.map((b) => ({ ...b, heightPercent: (b.count / maxCount) * 100 }));
  }, [summaries]);

  const pisBins = React.useMemo(() => {
    const values = summaries.map((s) => s.pis_percent);
    if (values.length === 0) return [];
    const min = Math.min(...values);
    const max = Math.max(...values);
    const numBins = 12;
    const binSize = (max - min) / numBins || 1;
    const bins = Array.from({ length: numBins }, (_, i) => ({
      start: min + i * binSize,
      end: min + (i + 1) * binSize,
      count: 0,
    }));
    for (const v of values) {
      const idx = Math.min(numBins - 1, Math.floor((v - min) / binSize));
      if (idx >= 0 && idx < numBins) bins[idx].count++;
    }
    const maxCount = Math.max(...bins.map((b) => b.count), 1);
    return bins.map((b) => ({ ...b, heightPercent: (b.count / maxCount) * 100 }));
  }, [summaries]);

  const gapBins = React.useMemo(() => {
    const values = summaries.map((s) => s.gap_percent);
    if (values.length === 0) return [];
    const min = Math.min(...values);
    const max = Math.max(...values);
    const numBins = 12;
    const binSize = (max - min) / numBins || 1;
    const bins = Array.from({ length: numBins }, (_, i) => ({
      start: min + i * binSize,
      end: min + (i + 1) * binSize,
      count: 0,
    }));
    for (const v of values) {
      const idx = Math.min(numBins - 1, Math.floor((v - min) / binSize));
      if (idx >= 0 && idx < numBins) bins[idx].count++;
    }
    const maxCount = Math.max(...bins.map((b) => b.count), 1);
    return bins.map((b) => ({ ...b, heightPercent: (b.count / maxCount) * 100 }));
  }, [summaries]);

  const divBins = React.useMemo(() => {
    const values = summaries.map((s) => s.mean_divergence * 100);
    if (values.length === 0) return [];
    const min = Math.min(...values);
    const max = Math.max(...values);
    const numBins = 12;
    const binSize = (max - min) / numBins || 1;
    const bins = Array.from({ length: numBins }, (_, i) => ({
      start: min + i * binSize,
      end: min + (i + 1) * binSize,
      count: 0,
    }));
    for (const v of values) {
      const idx = Math.min(numBins - 1, Math.floor((v - min) / binSize));
      if (idx >= 0 && idx < numBins) bins[idx].count++;
    }
    const maxCount = Math.max(...bins.map((b) => b.count), 1);
    return bins.map((b) => ({ ...b, heightPercent: (b.count / maxCount) * 100 }));
  }, [summaries]);

  const additionalStats = React.useMemo(() => {
    const total = summaries.length;
    const rawColumns = summaries.reduce(
      (sum, summary) => sum + (summary.raw_length || summary.length),
      0
    );
    const retainedColumns = summaries.reduce((sum, summary) => sum + summary.length, 0);
    const matrixCells = summaries.reduce(
      (sum, summary) => sum + summary.num_taxa * summary.length,
      0
    );
    const totalVariableSites = summaries.reduce(
      (sum, summary) => sum + summary.variable_count,
      0
    );
    return {
      passRate: total > 0 ? (overview.passed_alignments / total) * 100 : 0,
      matrixCompleteness:
        matrixCells > 0
          ? (summaries.reduce((sum, summary) => sum + summary.total_basepairs, 0) /
              matrixCells) *
            100
          : 0,
      columnsRetained: rawColumns > 0 ? (retainedColumns / rawColumns) * 100 : 0,
      retainedColumns,
      rawColumns,
      meanVariablePercent:
        total > 0
          ? summaries.reduce((sum, summary) => sum + summary.variable_percent, 0) / total
          : 0,
      totalVariableSites,
      meanGcPercent:
        total > 0
          ? summaries.reduce((sum, summary) => sum + summary.gc_percent, 0) / total
          : 0,
      meanDivergencePercent:
        total > 0
          ? (summaries.reduce((sum, summary) => sum + summary.mean_divergence, 0) /
              total) *
            100
          : 0,
    };
  }, [summaries, overview.passed_alignments]);

  const occupancyRows = React.useMemo(() => {
    const trimmed = new Map<string, { loci: number; basepairs: number }>();
    const acceptedOrfs = new Map<string, { loci: number; basepairs: number }>();

    const addSample = (
      destination: Map<string, { loci: number; basepairs: number }>,
      taxon: string,
      basepairs: number
    ) => {
      const current = destination.get(taxon) ?? { loci: 0, basepairs: 0 };
      current.loci += 1;
      current.basepairs += basepairs;
      destination.set(taxon, current);
    };

    for (const summary of summaries) {
      // General and ORF output are independent: Catalog pass controls the
      // trimming overlay, while ORF acceptance controls the ORF overlay.
      const retainedTaxa = new Set(summary.pass ? summary.retained_taxa ?? [] : []);
      const acceptedOrf = Boolean(
        summary.orf_evaluated && summary.orf_candidate_found && summary.orf_valid
      );
      for (const taxon of retainedTaxa) {
        const basepairs = summary.retained_taxon_basepairs?.[taxon] ?? 0;
        addSample(trimmed, taxon, basepairs);
      }
      if (acceptedOrf) {
        for (const taxon of new Set(summary.orf_retained_taxa ?? [])) {
          const basepairs = summary.orf_retained_taxon_basepairs?.[taxon] ?? 0;
          addSample(acceptedOrfs, taxon, basepairs);
        }
      }
    }

    return occupancy
      .map((sample) => ({
        name: sample.taxon_name,
        rawLoci: sample.present_loci_count,
        rawBasepairs: sample.total_bp,
        trimmedLoci: trimmed.get(sample.taxon_name)?.loci ?? 0,
        trimmedBasepairs: trimmed.get(sample.taxon_name)?.basepairs ?? 0,
        orfLoci: acceptedOrfs.get(sample.taxon_name)?.loci ?? 0,
        orfBasepairs: acceptedOrfs.get(sample.taxon_name)?.basepairs ?? 0,
      }))
      .sort((left, right) =>
        occupancyMetric === 'loci'
          ? right.rawLoci - left.rawLoci || left.name.localeCompare(right.name)
          : right.rawBasepairs - left.rawBasepairs || left.name.localeCompare(right.name)
      );
  }, [occupancy, summaries, occupancyMetric]);

  const occupancyMaximum = React.useMemo(
    () =>
      Math.max(
        1,
        ...occupancyRows.map((row) =>
          occupancyMetric === 'loci' ? row.rawLoci : row.rawBasepairs
        )
      ),
    [occupancyRows, occupancyMetric]
  );

  if (summaries.length === 0) {
    return (
      <div className="flex-1 flex items-center justify-center bg-[#0e1014] text-[#8b949e]">
        No dataset loaded
      </div>
    );
  }

  return (
    <div className="flex-1 flex flex-col h-full bg-[#0e1014] overflow-y-auto p-6 select-none space-y-6">
      <div className="flex items-center justify-between gap-4">
        <div>
          <h2 className="text-sm font-semibold text-[#dce6ff]">QC Statistics</h2>
          <p className="mt-0.5 text-[11px] text-[#8b949e]">
            Current post-filter alignment measurements and pass/fail results
          </p>
        </div>
        <div className="flex items-center gap-2">
          {exportNotice && (
            <span className="flex items-center gap-1 text-[10px] text-emerald-300">
              <CheckCircle2 className="w-3 h-3" />
              {exportNotice}
            </span>
          )}
          <button
            onClick={handleExportStats}
            disabled={isExporting}
            className="flex items-center gap-1.5 rounded-md border border-blue-500/30 bg-blue-500/10 px-2.5 py-1.5 text-[11px] font-medium text-blue-300 hover:bg-blue-500/15 disabled:cursor-wait disabled:opacity-50"
          >
            {isExporting ? (
              <Loader2 className="w-3.5 h-3.5 animate-spin" />
            ) : (
              <Download className="w-3.5 h-3.5" />
            )}
            Export Alignment Stats CSV
          </button>
        </div>
      </div>

      {/* Top Metric Cards */}
      <div className="grid grid-cols-1 gap-4 sm:grid-cols-2 xl:grid-cols-5">
        <div className="bg-[#14171d] border border-[#232833] rounded-xl p-4 flex flex-col justify-between shadow-sm">
          <div className="flex items-center justify-between text-[#8b949e] text-xs">
            <span>Total Alignments</span>
            <Layers className="w-4 h-4 text-blue-400" />
          </div>
          <div className="mt-2 text-2xl font-bold text-[#dce6ff] font-mono">
            {overview.total_alignments.toLocaleString()}
          </div>
          <div className="mt-1 text-[11px] text-emerald-400 flex items-center gap-1 font-medium">
            <span>{overview.passed_alignments.toLocaleString()} passing gates</span>
          </div>
        </div>

        <div className="bg-[#14171d] border border-[#232833] rounded-xl p-4 flex flex-col justify-between shadow-sm">
          <div className="flex items-center justify-between text-[#8b949e] text-xs">
            <span>Mean Taxa / Locus</span>
            <UsersRound className="w-4 h-4 text-cyan-400" />
          </div>
          <div className="mt-2 text-2xl font-bold text-[#dce6ff] font-mono">
            {overview.mean_taxa.toFixed(1)}
          </div>
          <div className="mt-1 text-[11px] text-[#8b949e]">
            {overview.total_unique_taxa} total species
          </div>
        </div>

        <div className="bg-[#14171d] border border-[#232833] rounded-xl p-4 flex flex-col justify-between shadow-sm">
          <div className="flex items-center justify-between text-[#8b949e] text-xs">
            <span>Mean Locus Length</span>
            <Ruler className="w-4 h-4 text-purple-400" />
          </div>
          <div className="mt-2 text-2xl font-bold text-[#dce6ff] font-mono">
            {overview.mean_length.toFixed(0)} bp
          </div>
          <div className="mt-1 text-[11px] text-[#8b949e]">
            {(overview.total_matrix_basepairs / 1_000_000).toFixed(2)} MB retained sequence
          </div>
        </div>

        <div className="bg-[#14171d] border border-[#232833] rounded-xl p-4 flex flex-col justify-between shadow-sm">
          <div className="flex items-center justify-between text-[#8b949e] text-xs">
            <span>Mean PIS Sites</span>
            <BarChart2 className="w-4 h-4 text-amber-400" />
          </div>
          <div className="mt-2 text-2xl font-bold text-[#dce6ff] font-mono">
            {overview.mean_pis.toFixed(0)}
          </div>
          <div className="mt-1 text-[11px] text-amber-400 font-medium">
            {((overview.mean_pis / (overview.mean_length || 1)) * 100).toFixed(1)}% of matrix
          </div>
        </div>

        <div className="bg-[#14171d] border border-[#232833] rounded-xl p-4 flex flex-col justify-between shadow-sm">
          <div className="flex items-center justify-between text-[#8b949e] text-xs">
            <span>Mean Gap %</span>
            <Minus className="w-4 h-4 text-rose-400" />
          </div>
          <div className="mt-2 text-2xl font-bold text-[#dce6ff] font-mono">
            {overview.mean_gap_percent.toFixed(1)}%
          </div>
          <div className="mt-1 text-[11px] text-[#8b949e]">across current alignments</div>
        </div>
      </div>

      {/* Additional dataset-wide summaries derived from completed catalog metrics. */}
      <div className="grid grid-cols-1 gap-3 sm:grid-cols-2 xl:grid-cols-3">
        <div className="rounded-xl border border-[#232833] bg-[#14171d] p-3.5 shadow-sm">
          <div className="flex items-center justify-between text-[11px] text-[#8b949e]">
            <span>Pass Rate</span>
            <CheckCircle2 className="h-4 w-4 text-emerald-400" />
          </div>
          <div className="mt-1.5 font-mono text-xl font-bold text-[#dce6ff]">
            {additionalStats.passRate.toFixed(1)}%
          </div>
          <div className="mt-0.5 text-[10px] text-[#8b949e]">
            {overview.discarded_alignments.toLocaleString()} failing alignments
          </div>
        </div>

        <div className="rounded-xl border border-[#232833] bg-[#14171d] p-3.5 shadow-sm">
          <div className="flex items-center justify-between text-[11px] text-[#8b949e]">
            <span>Matrix Completeness</span>
            <Layers className="h-4 w-4 text-cyan-400" />
          </div>
          <div className="mt-1.5 font-mono text-xl font-bold text-[#dce6ff]">
            {additionalStats.matrixCompleteness.toFixed(1)}%
          </div>
          <div className="mt-0.5 text-[10px] text-[#8b949e]">non-missing cells after processing</div>
        </div>

        <div className="rounded-xl border border-[#232833] bg-[#14171d] p-3.5 shadow-sm">
          <div className="flex items-center justify-between text-[11px] text-[#8b949e]">
            <span>Columns Retained</span>
            <Scissors className="h-4 w-4 text-blue-400" />
          </div>
          <div className="mt-1.5 font-mono text-xl font-bold text-[#dce6ff]">
            {additionalStats.columnsRetained.toFixed(1)}%
          </div>
          <div className="mt-0.5 text-[10px] text-[#8b949e]">
            {additionalStats.retainedColumns.toLocaleString()} / {additionalStats.rawColumns.toLocaleString()} locus columns
          </div>
        </div>

        <div className="rounded-xl border border-[#232833] bg-[#14171d] p-3.5 shadow-sm">
          <div className="flex items-center justify-between text-[11px] text-[#8b949e]">
            <span>Mean Variable Sites</span>
            <BarChart2 className="h-4 w-4 text-amber-400" />
          </div>
          <div className="mt-1.5 font-mono text-xl font-bold text-[#dce6ff]">
            {additionalStats.meanVariablePercent.toFixed(1)}%
          </div>
          <div className="mt-0.5 text-[10px] text-[#8b949e]">
            {additionalStats.totalVariableSites.toLocaleString()} variable sites total
          </div>
        </div>

        <div className="rounded-xl border border-[#232833] bg-[#14171d] p-3.5 shadow-sm">
          <div className="flex items-center justify-between text-[11px] text-[#8b949e]">
            <span>Mean GC Content</span>
            <Dna className="h-4 w-4 text-green-400" />
          </div>
          <div className="mt-1.5 font-mono text-xl font-bold text-[#dce6ff]">
            {additionalStats.meanGcPercent.toFixed(1)}%
          </div>
          <div className="mt-0.5 text-[10px] text-[#8b949e]">across current alignments</div>
        </div>

        <div className="rounded-xl border border-[#232833] bg-[#14171d] p-3.5 shadow-sm">
          <div className="flex items-center justify-between text-[11px] text-[#8b949e]">
            <span>Mean Divergence</span>
            <Activity className="h-4 w-4 text-rose-400" />
          </div>
          <div className="mt-1.5 font-mono text-xl font-bold text-[#dce6ff]">
            {additionalStats.meanDivergencePercent.toFixed(2)}%
          </div>
          <div className="mt-0.5 text-[10px] text-[#8b949e]">mean distance from locus consensus</div>
        </div>
      </div>

      {/* Sample occupancy and retained-output overlays. */}
      <section className="rounded-xl border border-[#232833] bg-[#14171d] p-5 shadow-sm">
        <div className="flex flex-wrap items-start justify-between gap-4">
          <div>
            <h3 className="text-xs font-semibold uppercase tracking-wider text-[#dce6ff]">
              Sample Occupancy
            </h3>
            <p className="mt-1 max-w-2xl text-[10px] leading-relaxed text-[#8b949e]">
              Full bars show each sample's raw dataset coverage. Add an overlay to compare loci
              that remain after trimming and alignment gates, or the accepted-ORF subset.
            </p>
          </div>

          <div className="flex flex-wrap items-center justify-end gap-2">
            <div className="flex items-center rounded-lg border border-[#2d3545] bg-[#0e1014] p-0.5">
              <button
                type="button"
                onClick={() => setOccupancyMetric('loci')}
                className={`whitespace-nowrap rounded-md px-2.5 py-1 text-[10px] font-medium transition-colors ${
                  occupancyMetric === 'loci'
                    ? 'bg-blue-600 text-white'
                    : 'text-[#8b949e] hover:text-[#dce6ff]'
                }`}
              >
                Number of Loci
              </button>
              <button
                type="button"
                onClick={() => setOccupancyMetric('basepairs')}
                className={`whitespace-nowrap rounded-md px-2.5 py-1 text-[10px] font-medium transition-colors ${
                  occupancyMetric === 'basepairs'
                    ? 'bg-blue-600 text-white'
                    : 'text-[#8b949e] hover:text-[#dce6ff]'
                }`}
              >
                Base Pairs
              </button>
            </div>

            <div className="flex items-center rounded-lg border border-[#2d3545] bg-[#0e1014] p-0.5">
              <button
                type="button"
                onClick={() => setOccupancyOverlay('none')}
                className={`rounded-md px-2.5 py-1 text-[10px] font-medium transition-colors ${
                  occupancyOverlay === 'none'
                    ? 'bg-slate-600 text-white'
                    : 'text-[#8b949e] hover:text-[#dce6ff]'
                }`}
              >
                Raw Only
              </button>
              <button
                type="button"
                onClick={() => setOccupancyOverlay('trimming')}
                className={`rounded-md px-2.5 py-1 text-[10px] font-medium transition-colors ${
                  occupancyOverlay === 'trimming'
                    ? 'bg-emerald-600 text-white'
                    : 'text-[#8b949e] hover:text-[#dce6ff]'
                }`}
              >
                Trimming
              </button>
              <button
                type="button"
                disabled={!orfEnabled}
                onClick={() => setOccupancyOverlay('orf')}
                className={`rounded-md px-2.5 py-1 text-[10px] font-medium transition-colors disabled:cursor-not-allowed disabled:opacity-35 ${
                  occupancyOverlay === 'orf'
                    ? 'bg-violet-600 text-white'
                    : 'text-[#8b949e] hover:text-[#dce6ff]'
                }`}
                title={orfEnabled ? 'Show samples retained in accepted ORF loci' : 'Enable ORF analysis to use this overlay'}
              >
                ORF
              </button>
            </div>
          </div>
        </div>

        <div className="mt-4 flex items-center gap-4 text-[10px] text-[#8b949e]">
          <span className="flex items-center gap-1.5">
            <span className="h-2.5 w-2.5 rounded-sm bg-slate-500/55" /> Raw total
          </span>
          {occupancyOverlay !== 'none' && (
            <span className="flex items-center gap-1.5">
              <span
                className={`h-2.5 w-2.5 rounded-sm ${
                  occupancyOverlay === 'orf' ? 'bg-violet-500' : 'bg-emerald-500'
                }`}
              />
              {occupancyOverlay === 'orf' ? 'Accepted ORFs' : 'Remaining after trimming'}
            </span>
          )}
          <span className="ml-auto font-mono">
            Scale: 0–{occupancyMaximum.toLocaleString()} {occupancyMetric === 'loci' ? 'loci' : 'bp'}
          </span>
        </div>

        <div className="mt-3 max-h-[480px] overflow-y-auto rounded-lg border border-[#232833] bg-[#101319]">
          {occupancyRows.map((row) => {
            const rawValue = occupancyMetric === 'loci' ? row.rawLoci : row.rawBasepairs;
            const overlayValue =
              occupancyOverlay === 'orf'
                ? occupancyMetric === 'loci'
                  ? row.orfLoci
                  : row.orfBasepairs
                : occupancyMetric === 'loci'
                  ? row.trimmedLoci
                  : row.trimmedBasepairs;
            const rawWidth = (rawValue / occupancyMaximum) * 100;
            const overlayWidth = (overlayValue / occupancyMaximum) * 100;
            const unit = occupancyMetric === 'loci' ? 'loci' : 'bp';
            return (
              <div
                key={row.name}
                className="grid min-h-9 grid-cols-[minmax(130px,240px)_minmax(180px,1fr)_135px] items-center gap-3 border-b border-[#1c222c] px-3 last:border-b-0"
              >
                <span className="truncate font-mono text-[10px] text-[#c9d1d9]" title={row.name}>
                  {row.name}
                </span>
                <div
                  className="relative h-3 overflow-hidden rounded-sm bg-[#1b2029]"
                  title={`${row.name}: ${rawValue.toLocaleString()} raw ${unit}${
                    occupancyOverlay !== 'none'
                      ? `; ${overlayValue.toLocaleString()} ${
                          occupancyOverlay === 'orf' ? 'accepted ORF' : 'remaining'
                        } ${unit}`
                      : ''
                  }`}
                >
                  <div
                    className="absolute inset-y-0 left-0 bg-slate-500/55"
                    style={{ width: `${rawWidth}%` }}
                  />
                  {occupancyOverlay !== 'none' && (
                    <div
                      className={`absolute inset-y-0 left-0 ${
                        occupancyOverlay === 'orf' ? 'bg-violet-500' : 'bg-emerald-500'
                      }`}
                      style={{ width: `${overlayWidth}%` }}
                    />
                  )}
                </div>
                <span
                  className="text-right font-mono text-[10px] text-[#8b949e]"
                  title={`${overlayValue.toLocaleString()} / ${rawValue.toLocaleString()} ${unit}`}
                >
                  {occupancyOverlay === 'none'
                    ? compactNumberFormatter.format(rawValue)
                    : `${compactNumberFormatter.format(overlayValue)} / ${compactNumberFormatter.format(rawValue)}`}
                  {' '}{unit}
                </span>
              </div>
            );
          })}
        </div>
      </section>

      {/* 4 Multi-Metric Interactive Histograms */}
      <div className="grid grid-cols-1 gap-6 xl:grid-cols-2">
        {/* 1. Alignment Length Distribution */}
        <div className="bg-[#14171d] border border-[#232833] rounded-xl p-5 shadow-sm flex flex-col">
          <div className="flex items-center justify-between mb-4">
            <span className="font-semibold text-xs text-[#dce6ff] uppercase tracking-wider">
              Alignment Length Distribution (bp)
            </span>
            <span className="text-xs font-mono text-[#8b949e]">
              Mean: {overview.mean_length.toFixed(0)} bp
            </span>
          </div>
          <div className="h-44 flex items-end gap-1.5 pt-4 pb-2 border-b border-[#232833]">
            {lengthBins.map((bin, idx) => (
              <div key={idx} className="flex-1 flex flex-col items-center h-full justify-end group relative">
                <div
                  className="w-full bg-purple-500/70 hover:bg-purple-400 rounded-t transition-all duration-200"
                  style={{ height: `${Math.max(4, bin.heightPercent)}%` }}
                />
                {/* Tooltip */}
                <div className="absolute -top-7 opacity-0 group-hover:opacity-100 bg-[#1f242e] text-[10px] font-mono text-[#dce6ff] px-1.5 py-0.5 rounded border border-[#2d3545] pointer-events-none transition-opacity">
                  {bin.count} loci ({bin.start.toFixed(0)}-{bin.end.toFixed(0)} bp)
                </div>
              </div>
            ))}
          </div>
          <div className="flex justify-between text-[11px] font-mono text-[#8b949e] mt-2">
            <span>{lengthBins[0]?.start.toFixed(0)} bp</span>
            <span>{lengthBins[lengthBins.length - 1]?.end.toFixed(0)} bp</span>
          </div>
        </div>

        {/* 2. Parsimony-Informative Sites (%) */}
        <div className="bg-[#14171d] border border-[#232833] rounded-xl p-5 shadow-sm flex flex-col">
          <div className="flex items-center justify-between mb-4">
            <span className="font-semibold text-xs text-[#dce6ff] uppercase tracking-wider">
              Parsimony-Informative Sites (% of Length)
            </span>
            <span className="text-xs font-mono text-[#8b949e]">
              Mean: {overview.mean_pis.toFixed(0)} sites
            </span>
          </div>
          <div className="h-44 flex items-end gap-1.5 pt-4 pb-2 border-b border-[#232833]">
            {pisBins.map((bin, idx) => (
              <div key={idx} className="flex-1 flex flex-col items-center h-full justify-end group relative">
                <div
                  className="w-full bg-amber-500/70 hover:bg-amber-400 rounded-t transition-all duration-200"
                  style={{ height: `${Math.max(4, bin.heightPercent)}%` }}
                />
                <div className="absolute -top-7 opacity-0 group-hover:opacity-100 bg-[#1f242e] text-[10px] font-mono text-[#dce6ff] px-1.5 py-0.5 rounded border border-[#2d3545] pointer-events-none transition-opacity">
                  {bin.count} loci ({bin.start.toFixed(1)}-{bin.end.toFixed(1)}%)
                </div>
              </div>
            ))}
          </div>
          <div className="flex justify-between text-[11px] font-mono text-[#8b949e] mt-2">
            <span>{pisBins[0]?.start.toFixed(1)}%</span>
            <span>{pisBins[pisBins.length - 1]?.end.toFixed(1)}%</span>
          </div>
        </div>

        {/* 3. Gap Percentage Distribution */}
        <div className="bg-[#14171d] border border-[#232833] rounded-xl p-5 shadow-sm flex flex-col">
          <div className="flex items-center justify-between mb-4">
            <span className="font-semibold text-xs text-[#dce6ff] uppercase tracking-wider">
              Gap & Missing Data Distribution (%)
            </span>
            <span className="text-xs font-mono text-[#8b949e]">
              Mean: {overview.mean_gap_percent.toFixed(1)}%
            </span>
          </div>
          <div className="h-44 flex items-end gap-1.5 pt-4 pb-2 border-b border-[#232833]">
            {gapBins.map((bin, idx) => (
              <div key={idx} className="flex-1 flex flex-col items-center h-full justify-end group relative">
                <div
                  className="w-full bg-rose-500/70 hover:bg-rose-400 rounded-t transition-all duration-200"
                  style={{ height: `${Math.max(4, bin.heightPercent)}%` }}
                />
                <div className="absolute -top-7 opacity-0 group-hover:opacity-100 bg-[#1f242e] text-[10px] font-mono text-[#dce6ff] px-1.5 py-0.5 rounded border border-[#2d3545] pointer-events-none transition-opacity">
                  {bin.count} loci ({bin.start.toFixed(1)}-{bin.end.toFixed(1)}%)
                </div>
              </div>
            ))}
          </div>
          <div className="flex justify-between text-[11px] font-mono text-[#8b949e] mt-2">
            <span>0% Gaps</span>
            <span>100% Gaps</span>
          </div>
        </div>

        {/* 4. Genetic Divergence Distribution */}
        <div className="bg-[#14171d] border border-[#232833] rounded-xl p-5 shadow-sm flex flex-col">
          <div className="flex items-center justify-between mb-4">
            <span className="font-semibold text-xs text-[#dce6ff] uppercase tracking-wider">
              Pairwise Divergence from Consensus (%)
            </span>
            <span className="text-xs font-mono text-[#8b949e]">Paralog Cutoff Monitor</span>
          </div>
          <div className="h-44 flex items-end gap-1.5 pt-4 pb-2 border-b border-[#232833]">
            {divBins.map((bin, idx) => (
              <div key={idx} className="flex-1 flex flex-col items-center h-full justify-end group relative">
                <div
                  className="w-full bg-cyan-500/70 hover:bg-cyan-400 rounded-t transition-all duration-200"
                  style={{ height: `${Math.max(4, bin.heightPercent)}%` }}
                />
                <div className="absolute -top-7 opacity-0 group-hover:opacity-100 bg-[#1f242e] text-[10px] font-mono text-[#dce6ff] px-1.5 py-0.5 rounded border border-[#2d3545] pointer-events-none transition-opacity">
                  {bin.count} loci ({bin.start.toFixed(1)}-{bin.end.toFixed(1)}% div)
                </div>
              </div>
            ))}
          </div>
          <div className="flex justify-between text-[11px] font-mono text-[#8b949e] mt-2">
            <span>{divBins[0]?.start.toFixed(1)}%</span>
            <span>{divBins[divBins.length - 1]?.end.toFixed(1)}%</span>
          </div>
        </div>
      </div>
    </div>
  );
});

export default QcDistributions;
