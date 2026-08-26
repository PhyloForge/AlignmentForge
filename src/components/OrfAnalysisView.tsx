import React, { useEffect, useMemo, useRef, useState } from 'react';
import { ArrowUpDown, CheckCircle, Dna, Eye, Loader2, Search, XCircle, CheckSquare } from 'lucide-react';
import { AlignmentSummary, OrfSearchMode } from '../types';
import { shouldSkipOrfLocus } from '../parsers/clientParser';

interface OrfAnalysisViewProps {
  summaries: AlignmentSummary[];
  enabled: boolean;
  isProcessing: boolean;
  processingPercent: number;
  searchTerm: string;
  onSearchTermChange: (term: string) => void;
  statusFilter: 'all' | 'accepted' | 'discarded';
  onStatusFilterChange: (filter: 'all' | 'accepted' | 'discarded') => void;
  orfSearchMode?: OrfSearchMode;
  skipNonCodingOrf?: boolean;
  onSelectLocus: (id: string, filePath: string) => void;
  selectedPaths: Set<string>;
  onSelectPath: (path: string, selected: boolean) => void;
  onSelectAllPaths: () => void;
  onClearSelectedPaths: () => void;
}

const ROW_HEIGHT = 38;
const OVERSCAN = 15;
const COLUMN_COUNT = 13;

type OrfSortField =
  | 'id'
  | 'candidate'
  | 'reference'
  | 'frame'
  | 'range'
  | 'support'
  | 'retained'
  | 'length'
  | 'coding'
  | 'intron'
  | 'status';

function getCandidateStatus(
  summary: AlignmentSummary,
  skipNonCodingOrf: boolean,
  orfSearchMode: OrfSearchMode
): 'Accepted' | 'Discarded' | 'Not assessed' | 'Skipped' {
  if (skipNonCodingOrf && shouldSkipOrfLocus(summary.id, orfSearchMode)) {
    return 'Skipped';
  }
  if (summary.orf_evaluated) {
    return summary.orf_candidate_found && summary.orf_valid ? 'Accepted' : 'Discarded';
  }
  if (summary.orf_reference_evaluated && !summary.orf_reference_matched) {
    return 'Discarded';
  }
  return 'Not assessed';
}

const OrfAnalysisView: React.FC<OrfAnalysisViewProps> = ({
  summaries,
  enabled,
  isProcessing,
  processingPercent,
  searchTerm,
  onSearchTermChange,
  statusFilter,
  onStatusFilterChange,
  orfSearchMode = 'continuouscds',
  skipNonCodingOrf = true,
  onSelectLocus,
  selectedPaths,
  onSelectPath,
  onSelectAllPaths,
  onClearSelectedPaths,
}) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const [scrollTop, setScrollTop] = useState(0);
  const [containerHeight, setContainerHeight] = useState(600);

  const [sortField, setSortField] = useState<OrfSortField>('id');
  const [sortAsc, setSortAsc] = useState(true);

  const rows = useMemo(() => {
    const term = searchTerm.toLowerCase();
    const filtered = summaries.filter((summary) => {
      const matchesSearch =
        !term || summary.id.toLowerCase().includes(term) || summary.file_name.toLowerCase().includes(term);

      const status = getCandidateStatus(summary, skipNonCodingOrf, orfSearchMode);
      const matchesFilter = statusFilter === 'all' || status.toLowerCase() === statusFilter;
      return matchesSearch && matchesFilter;
    });

    const candidateOrder = { 'Not assessed': 0, Skipped: 1, Discarded: 2, Accepted: 3 };

    const getSortValue = (summary: AlignmentSummary) => {
      switch (sortField) {
        case 'id':
          return summary.id;
        case 'candidate':
          return candidateOrder[getCandidateStatus(summary, skipNonCodingOrf, orfSearchMode)];
        case 'reference':
          if (!summary.orf_reference_evaluated) return 0;
          return summary.orf_reference_matched ? 2 : 1;
        case 'frame':
          return summary.orf_frame ?? Number.NEGATIVE_INFINITY;
        case 'range':
          return summary.orf_start ?? Number.NEGATIVE_INFINITY;
        case 'support':
          return summary.orf_support_percent ?? Number.NEGATIVE_INFINITY;
        case 'retained':
          return summary.orf_retained_samples ?? Number.NEGATIVE_INFINITY;
        case 'length':
          return summary.orf_candidate_length_aa ?? Number.NEGATIVE_INFINITY;
        case 'coding':
          return summary.orf_coding_score ?? Number.NEGATIVE_INFINITY;
        case 'intron':
          return summary.orf_intron_length ?? Number.NEGATIVE_INFINITY;
        case 'status':
          return summary.pass ? 1 : 0;
      }
    };

    return filtered.sort((a, b) => {
      const valA = getSortValue(a);
      const valB = getSortValue(b);

      const comparison = typeof valA === 'string' ? valA.localeCompare(String(valB)) : (valA as number) - (valB as number);
      return sortAsc ? comparison : -comparison;
    });
  }, [summaries, searchTerm, statusFilter, sortField, sortAsc, skipNonCodingOrf, orfSearchMode]);

  const counts = useMemo(() => {
    let accepted = 0;
    let discarded = 0;
    for (const summary of summaries) {
      const status = getCandidateStatus(summary, skipNonCodingOrf, orfSearchMode);
      if (status === 'Accepted') accepted++;
      if (status === 'Discarded') discarded++;
    }
    return { all: summaries.length, accepted, discarded };
  }, [summaries, skipNonCodingOrf, orfSearchMode]);

  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;
    const ro = new ResizeObserver((entries) => {
      if (entries[0]) setContainerHeight(entries[0].contentRect.height);
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, [enabled]);

  useEffect(() => {
    setScrollTop(0);
    if (containerRef.current) {
      containerRef.current.scrollTop = 0;
    }
  }, [searchTerm, statusFilter, summaries.length, sortField, sortAsc]);

  const handleSort = (field: OrfSortField) => {
    if (sortField === field) {
      setSortAsc((prev) => !prev);
    } else {
      setSortField(field);
      setSortAsc(true);
    }
  };

  const sortableHeader = (label: string, field: OrfSortField, className: string = '') => (
    <th onClick={() => handleSort(field)} className={`cursor-pointer px-3 py-2.5 hover:text-[#c9d1d9] ${className}`}>
      <span className="inline-flex items-center gap-1">
        <span>{label}</span>
        <ArrowUpDown className="h-3 w-3" />
      </span>
    </th>
  );

  const allSelected = rows.length > 0 && rows.every(r => selectedPaths.has(r.file_path));

  const startIndex = Math.max(0, Math.floor(scrollTop / ROW_HEIGHT) - OVERSCAN);
  const endIndex = Math.min(rows.length, Math.ceil((scrollTop + containerHeight) / ROW_HEIGHT) + OVERSCAN);
  const visibleRows = rows.slice(startIndex, endIndex);
  const topPadding = startIndex * ROW_HEIGHT;
  const bottomPadding = Math.max(0, (rows.length - endIndex) * ROW_HEIGHT);

  if (!enabled) {
    return (
      <div className="flex-1 flex items-center justify-center bg-[#0e1014] p-8 text-center">
        <div className="max-w-md rounded-xl border border-pink-500/20 bg-[#14171d] p-8 shadow-xl">
          <Dna className="mx-auto h-10 w-10 text-pink-500" />
          <h2 className="mt-3 text-base font-semibold text-[#dce6ff]">ORF Analysis is off</h2>
          <p className="mt-2 text-xs leading-relaxed text-[#8b949e]">
            Enable it in the ORF Analysis sidebar when you want to extract candidate coding regions, inspect frames, or use exact-name exon references.
          </p>
        </div>
      </div>
    );
  }

  return (
    <div className="flex-1 min-w-0 flex flex-col bg-[#0e1014] overflow-hidden select-none">
      <div className="h-11 px-4 bg-[#14171d] border-b border-[#232833] flex items-center justify-between gap-4 flex-none">
        <div className="flex items-center gap-3">
          <div className="flex items-center gap-2 bg-[#171b22] px-2.5 py-1 rounded-md border border-[#232833] w-64">
            <Search className="w-3.5 h-3.5 text-[#8b949e]" />
            <input value={searchTerm} onChange={(event) => onSearchTermChange(event.target.value)} placeholder="Search ORF loci…" className="w-full bg-transparent text-xs text-[#c9d1d9] outline-none" />
          </div>
          <div className="flex items-center p-0.5 bg-[#0e1014] rounded-md border border-[#232833]">
            <button onClick={() => onStatusFilterChange('all')} className={`px-2.5 py-1 rounded text-xs ${statusFilter === 'all' ? 'bg-[#232a36] text-[#dce6ff] font-medium' : 'text-[#8b949e] hover:text-[#c9d1d9]'}`}>All ({counts.all})</button>
            <button onClick={() => onStatusFilterChange('accepted')} className={`px-2.5 py-1 rounded text-xs ${statusFilter === 'accepted' ? 'bg-emerald-500/20 text-emerald-300 font-medium' : 'text-[#8b949e] hover:text-emerald-400'}`}>Accepted ({counts.accepted})</button>
            <button onClick={() => onStatusFilterChange('discarded')} className={`px-2.5 py-1 rounded text-xs ${statusFilter === 'discarded' ? 'bg-rose-500/20 text-rose-300 font-medium' : 'text-[#8b949e] hover:text-rose-400'}`}>Discarded ({counts.discarded})</button>
          </div>
        </div>
        <div className="flex items-center gap-3 text-[11px] text-[#8b949e]">
          {isProcessing && (
            <span className="relative inline-flex min-w-[190px] items-center overflow-hidden rounded border border-cyan-500/30 bg-cyan-500/10 px-2 py-1 text-cyan-200">
              <span className="absolute inset-y-0 left-0 bg-cyan-500/25" style={{ width: `${Math.max(0, Math.min(100, processingPercent))}%` }} />
              <span className="relative z-10 inline-flex w-full items-center gap-1.5"><Loader2 className="h-3.5 w-3.5 animate-spin" />Processing ORFs…<span className="ml-auto font-mono">{Math.round(processingPercent)}%</span></span>
            </span>
          )}
          {selectedPaths.size > 0 && (
            <span className="flex items-center gap-2 mr-2">
              <span className="flex items-center gap-1.5 px-2 py-0.5 rounded-full bg-blue-500/10 border border-blue-500/20 text-blue-400 text-[11px] font-medium">
                <CheckSquare className="w-3 h-3" />
                <span>
                  {selectedPaths.size}{' '}
                  <span className="text-blue-400/70 font-normal">
                    selected
                  </span>
                </span>
              </span>
            </span>
          )}
          <span><b className="text-[#c9d1d9]">{rows.length}</b> of {summaries.length} loci</span>
        </div>
      </div>

      <div ref={containerRef} onScroll={(event) => setScrollTop(event.currentTarget.scrollTop)} className="flex-1 overflow-auto relative">
        <table className="min-w-max w-full border-collapse text-left text-xs">
          <thead className="sticky top-0 z-10 border-b border-[#232833] bg-[#171b22] font-mono text-[10px] text-[#8b949e]">
            <tr>
              <th className="py-2.5 px-3 w-10 text-center">
                <input
                  type="checkbox"
                  checked={allSelected}
                  onChange={allSelected ? onClearSelectedPaths : onSelectAllPaths}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-blue-500"
                />
              </th>
              {sortableHeader('LOCUS', 'id')}
              {sortableHeader('ORF CANDIDATE', 'candidate')}
              {sortableHeader('REFERENCE', 'reference')}
              {sortableHeader('FRAME', 'frame')}
              {sortableHeader('RAW RANGE', 'range')}
              {sortableHeader('SUPPORT', 'support')}
              {sortableHeader('RETAINED', 'retained')}
              {sortableHeader('ORF LENGTH', 'length')}
              {sortableHeader('CODING EVIDENCE', 'coding')}
              {sortableHeader('INTRON', 'intron')}
              {sortableHeader('CATALOG QC', 'status')}
              <th className="px-3 py-2.5 text-center">VIEW</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-[#1b1f27]">
            {topPadding > 0 && <tr><td style={{ height: `${topPadding}px` }} colSpan={COLUMN_COUNT} /></tr>}
            {visibleRows.map((summary) => {
              const candidateStatus = getCandidateStatus(summary, skipNonCodingOrf, orfSearchMode);
              const isChecked = selectedPaths.has(summary.file_path);
              return (
                <tr key={summary.file_path} onClick={() => onSelectLocus(summary.id, summary.file_path)} style={{ height: `${ROW_HEIGHT}px` }} className="cursor-pointer text-[#c9d1d9] hover:bg-[#14171d]">
                  <td className="w-10 text-center" onClick={(e) => {
                    e.stopPropagation();
                    onSelectPath(summary.file_path, !isChecked);
                  }}>
                    <input
                      type="checkbox"
                      checked={isChecked}
                      onChange={() => {}}
                      className="rounded bg-[#1f242e] border-[#2d3545] text-blue-500 cursor-pointer"
                    />
                  </td>
                  <td className="max-w-[220px] truncate px-3 py-2 font-mono font-semibold text-[#dce6ff]">{summary.id}</td>
                  <td className="px-3 py-2 font-mono text-[11px]">
                    {candidateStatus === 'Accepted' ? <span className="text-emerald-400">Accepted</span> : candidateStatus === 'Discarded' ? <span className="text-rose-400">Discarded</span> : <span className="text-[#6e7681]">{candidateStatus}</span>}
                  </td>
                  <td className="px-3 py-2 font-mono text-[11px]">
                    {!summary.orf_reference_evaluated ? <span className="text-[#6e7681]">Not used</span> : summary.orf_reference_matched ? <span className="text-blue-300" title={`Identity ${Number(summary.orf_reference_identity ?? 0).toFixed(1)}%; coverage ${Number(summary.orf_reference_coverage ?? 0).toFixed(1)}%`}>Matched</span> : <span className="text-rose-300">No match</span>}
                  </td>
                  <td className="px-3 py-2 font-mono text-[#8b949e]">{summary.orf_frame ? `${summary.orf_frame > 0 ? '+' : ''}${summary.orf_frame}` : '—'}</td>
                  <td className="px-3 py-2 font-mono text-[#8b949e]">{summary.orf_start != null && summary.orf_end != null ? `${summary.orf_start + 1}–${summary.orf_end}` : '—'}</td>
                  <td className="px-3 py-2 font-mono text-[#8b949e]">{summary.orf_evaluated ? `${summary.orf_support_count ?? 0} (${Number(summary.orf_support_percent ?? 0).toFixed(0)}%)` : '—'}</td>
                  <td className="px-3 py-2 font-mono text-[#8b949e]">{summary.orf_evaluated ? summary.orf_retained_samples ?? 0 : '—'}</td>
                  <td className="px-3 py-2 font-mono text-[#8b949e]">{summary.orf_evaluated ? `${summary.orf_candidate_length_aa ?? 0} aa` : '—'}</td>
                  <td className="px-3 py-2 font-mono text-[#8b949e]" title={`AA conservation ${Number(summary.orf_amino_acid_conservation ?? 0).toFixed(1)}%; frame contrast ${Number(summary.orf_frame_contrast ?? 0).toFixed(1)}%`}>{summary.orf_evaluated ? `${Number(summary.orf_coding_score ?? 0).toFixed(1)}/100` : '—'}</td>
                  <td className="px-3 py-2 font-mono text-[#8b949e]">{summary.orf_reference_matched ? `${summary.orf_intron_length ?? 0} bp` : '—'}</td>
                  <td className="px-3 py-2">{summary.pass ? <span className="inline-flex items-center gap-1 text-emerald-400"><CheckCircle className="h-3.5 w-3.5" />Pass</span> : <span title={summary.fail_reasons.join(', ')} className="inline-flex items-center gap-1 text-rose-400"><XCircle className="h-3.5 w-3.5" />Fail</span>}</td>
                  <td className="px-3 py-2 text-center"><button onClick={(event) => { event.stopPropagation(); onSelectLocus(summary.id, summary.file_path); }} className="rounded p-1 text-[#8b949e] hover:bg-[#232833] hover:text-[#dce6ff]" title="Open in MSA Viewer"><Eye className="h-3.5 w-3.5" /></button></td>
                </tr>
              );
            })}
            {bottomPadding > 0 && <tr><td style={{ height: `${bottomPadding}px` }} colSpan={COLUMN_COUNT} /></tr>}
          </tbody>
        </table>
      </div>
    </div>
  );
};

export { OrfAnalysisView };
