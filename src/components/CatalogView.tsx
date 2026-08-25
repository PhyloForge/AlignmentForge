import React, { useState, useMemo, useRef, useEffect, useCallback } from 'react';
import {
  Search,
  CheckCircle,
  XCircle,
  Eye,
  ArrowUpDown,
  CheckSquare,
  Square,
  ChevronLeft,
  ChevronRight,
  Loader2,
} from 'lucide-react';
import { AlignmentSummary, CatalogSortField, OrfSearchMode } from '../types';
import { shouldSkipOrfLocus } from '../parsers/clientParser';

export interface CatalogViewProps {
  summaries: AlignmentSummary[];
  isProcessing?: boolean;
  processingPercent?: number;
  selectedLocusId: string | null;
  onSelectLocus: (id: string, filePath: string) => void;
  selectedPaths: Set<string>;
  onToggleSelectPath: (path: string) => void;
  onSelectAllPaths: () => void;
  onClearSelectedPaths: () => void;
  searchTerm?: string;
  onSearchTermChange?: (term: string) => void;
  statusFilter?: 'all' | 'pass' | 'fail';
  onStatusFilterChange?: (filter: 'all' | 'pass' | 'fail') => void;
  sortField?: CatalogSortField;
  sortAsc?: boolean;
  onSortChange?: (field: CatalogSortField, asc: boolean) => void;
  orfEnabled?: boolean;
  orfSearchMode?: OrfSearchMode;
  skipNonCodingOrf?: boolean;
}

const ROW_HEIGHT = 38; // px per row
const OVERSCAN = 15; // extra rows above/below visible viewport

export const CatalogView: React.FC<CatalogViewProps> = React.memo(({
  summaries,
  isProcessing = false,
  processingPercent = 0,
  selectedLocusId,
  onSelectLocus,
  selectedPaths,
  onToggleSelectPath,
  onSelectAllPaths,
  onClearSelectedPaths,
  searchTerm: controlledSearchTerm,
  onSearchTermChange,
  statusFilter: controlledStatusFilter,
  onStatusFilterChange,
  sortField: controlledSortField,
  sortAsc: controlledSortAsc,
  onSortChange,
  orfEnabled = false,
  orfSearchMode = 'continuouscds',
  skipNonCodingOrf = true,
}) => {
  const [internalSearchTerm, setInternalSearchTerm] = useState('');
  const [internalStatusFilter, setInternalStatusFilter] = useState<'all' | 'pass' | 'fail'>('all');
  const [internalSortField, setInternalSortField] = useState<CatalogSortField>('id');
  const [internalSortAsc, setInternalSortAsc] = useState<boolean>(true);

  const searchTerm = controlledSearchTerm !== undefined ? controlledSearchTerm : internalSearchTerm;
  const statusFilter = controlledStatusFilter !== undefined ? controlledStatusFilter : internalStatusFilter;
  const sortField = controlledSortField !== undefined ? controlledSortField : internalSortField;
  const sortAsc = controlledSortAsc !== undefined ? controlledSortAsc : internalSortAsc;

  const setSearchTerm = (term: string) => {
    if (onSearchTermChange) onSearchTermChange(term);
    else setInternalSearchTerm(term);
  };

  const setStatusFilter = (filter: 'all' | 'pass' | 'fail') => {
    if (onStatusFilterChange) onStatusFilterChange(filter);
    else setInternalStatusFilter(filter);
  };

  const containerRef = useRef<HTMLDivElement | null>(null);
  const [scrollTop, setScrollTop] = useState(0);
  const [containerHeight, setContainerHeight] = useState(600);

  // Update container height
  useEffect(() => {
    const el = containerRef.current;
    if (!el) return;
    const observer = new ResizeObserver((entries) => {
      if (entries[0]) {
        setContainerHeight(entries[0].contentRect.height);
      }
    });
    observer.observe(el);
    return () => observer.disconnect();
  }, []);

  const handleScroll = (e: React.UIEvent<HTMLDivElement>) => {
    setScrollTop(e.currentTarget.scrollTop);
  };

  const handleSort = (field: CatalogSortField) => {
    const nextAsc = sortField === field ? !sortAsc : true;
    if (onSortChange) {
      onSortChange(field, nextAsc);
    } else {
      setInternalSortField(field);
      setInternalSortAsc(nextAsc);
    }
  };

  const filteredAndSorted = useMemo(() => {
    const searchLower = searchTerm.toLowerCase();
    const result = summaries.filter((item) => {
      const matchesSearch =
        searchLower === '' ||
        item.id.toLowerCase().includes(searchLower) ||
        item.file_name.toLowerCase().includes(searchLower);

      const matchesStatus =
        statusFilter === 'all' ||
        (statusFilter === 'pass' && item.pass) ||
        (statusFilter === 'fail' && !item.pass);

      return matchesSearch && matchesStatus;
    });

    result.sort((a, b) => {
      let valA = a[sortField];
      let valB = b[sortField];
      if (valA === undefined || valA === null) valA = -Infinity;
      if (valB === undefined || valB === null) valB = -Infinity;
      if (typeof valA === 'string') {
        return sortAsc
          ? (valA as string).localeCompare(valB as string)
          : (valB as string).localeCompare(valA as string);
      }
      if (typeof valA === 'boolean') {
        return sortAsc ? (valA === valB ? 0 : valA ? -1 : 1) : valA === valB ? 0 : valA ? 1 : -1;
      }
      return sortAsc ? (valA as number) - (valB as number) : (valB as number) - (valA as number);
    });

    return result;
  }, [summaries, searchTerm, statusFilter, sortField, sortAsc]);

  const totalRows = filteredAndSorted.length;
  const totalHeight = totalRows * ROW_HEIGHT;

  // Virtual Window Calculation
  const startIndex = Math.max(0, Math.floor(scrollTop / ROW_HEIGHT) - OVERSCAN);
  const endIndex = Math.min(totalRows, Math.ceil((scrollTop + containerHeight) / ROW_HEIGHT) + OVERSCAN);
  const visibleItems = filteredAndSorted.slice(startIndex, endIndex);
  const topPadding = startIndex * ROW_HEIGHT;
  const bottomPadding = Math.max(0, (totalRows - endIndex) * ROW_HEIGHT);

  const allSelected =
    filteredAndSorted.length > 0 &&
    filteredAndSorted.every((item) => selectedPaths.has(item.file_path));

  return (
    <div className="flex-1 flex flex-col h-full bg-[#0e1014] overflow-hidden select-none">
      {/* Top Filter and Search Bar */}
      <div className="h-11 px-4 bg-[#14171d] border-b border-[#232833] flex items-center justify-between gap-4 flex-none">
        <div className="flex items-center gap-3">
          {/* Search */}
          <div className="flex items-center gap-2 bg-[#171b22] px-2.5 py-1 rounded-md border border-[#232833] w-64">
            <Search className="w-3.5 h-3.5 text-[#8b949e]" />
            <input
              type="text"
              placeholder="Search loci or file name…"
              value={searchTerm}
              onChange={(e) => setSearchTerm(e.target.value)}
              className="bg-transparent text-xs text-[#c9d1d9] outline-none w-full"
            />
          </div>

          {/* Status Filter Tabs */}
          <div className="flex items-center p-0.5 bg-[#0e1014] rounded-md border border-[#232833]">
            <button
              onClick={() => setStatusFilter('all')}
              className={`px-2.5 py-1 rounded text-xs transition-colors ${
                statusFilter === 'all'
                  ? 'bg-[#232a36] text-[#dce6ff] font-medium'
                  : 'text-[#8b949e] hover:text-[#c9d1d9]'
              }`}
            >
              All ({summaries.length})
            </button>
            <button
              onClick={() => setStatusFilter('pass')}
              className={`px-2.5 py-1 rounded text-xs transition-colors ${
                statusFilter === 'pass'
                  ? 'bg-emerald-500/20 text-emerald-300 font-medium'
                  : 'text-[#8b949e] hover:text-emerald-400'
              }`}
            >
              Pass ({summaries.filter((s) => s.pass).length})
            </button>
            <button
              onClick={() => setStatusFilter('fail')}
              className={`px-2.5 py-1 rounded text-xs transition-colors ${
                statusFilter === 'fail'
                  ? 'bg-rose-500/20 text-rose-300 font-medium'
                  : 'text-[#8b949e] hover:text-rose-400'
              }`}
            >
              Fail ({summaries.filter((s) => !s.pass).length})
            </button>
          </div>
        </div>

        {/* Selection stats */}
        <div className="flex items-center gap-3 text-xs text-[#8b949e]">
          {isProcessing && (
            <span
              className="relative inline-flex min-w-[190px] items-center overflow-hidden rounded border border-cyan-500/30 bg-cyan-500/10 px-2 py-1 text-cyan-200"
              role="status"
              aria-label={`Processing all loci, ${Math.round(processingPercent)} percent complete`}
            >
              <span
                className="absolute inset-y-0 left-0 bg-cyan-500/25 transition-[width] duration-200 ease-out"
                style={{ width: `${Math.max(0, Math.min(100, processingPercent))}%` }}
              />
              <span className="relative z-10 inline-flex w-full items-center gap-1.5">
                <Loader2 className="w-3.5 h-3.5 animate-spin" />
                <span>Processing all loci…</span>
                <span className="ml-auto font-mono text-[10px] font-semibold">
                  {Math.round(processingPercent)}%
                </span>
              </span>
            </span>
          )}
          <span>
            <b className="text-[#c9d1d9]">{filteredAndSorted.length}</b> of{' '}
            {summaries.length} loci
          </span>
          <span className="text-[#232833]">|</span>
          <button
            onClick={allSelected ? onClearSelectedPaths : onSelectAllPaths}
            className="flex items-center gap-1.5 text-blue-400 hover:text-blue-300 transition-colors"
          >
            {allSelected ? <CheckSquare className="w-3.5 h-3.5" /> : <Square className="w-3.5 h-3.5" />}
            <span>{selectedPaths.size} selected</span>
          </button>
        </div>
      </div>

      {/* Virtual Table Container */}
      <div
        ref={containerRef}
        onScroll={handleScroll}
        className="flex-1 overflow-y-auto overflow-x-auto relative"
      >
        <table className="w-full min-w-max text-left text-xs border-collapse">
          <thead className="bg-[#171b22] text-[#8b949e] sticky top-0 z-20 border-b border-[#232833] font-mono text-[11px]">
            <tr>
              <th className="py-2.5 px-3 w-10 text-center">
                <input
                  type="checkbox"
                  checked={allSelected}
                  onChange={allSelected ? onClearSelectedPaths : onSelectAllPaths}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-blue-500"
                />
              </th>
              <th
                onClick={() => handleSort('id')}
                className="py-2.5 px-3 cursor-pointer hover:text-[#c9d1d9]"
              >
                <div className="flex items-center gap-1">
                  <span>LOCUS ID</span>
                  <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th
                onClick={() => handleSort('num_taxa')}
                className="py-2.5 px-3 cursor-pointer hover:text-[#c9d1d9]"
              >
                <div className="flex items-center gap-1">
                  <span>TAXA</span>
                  <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th
                onClick={() => handleSort('length')}
                className="py-2.5 px-3 cursor-pointer hover:text-[#c9d1d9]"
              >
                <div className="flex items-center gap-1">
                  <span>LENGTH</span>
                  <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th
                onClick={() => handleSort('variable_count')}
                className="py-2.5 px-3 cursor-pointer hover:text-[#c9d1d9]"
              >
                <div className="flex items-center gap-1">
                  <span>VARIABLE SITES (%)</span>
                  <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th
                onClick={() => handleSort('pis_count')}
                className="py-2.5 px-3 cursor-pointer hover:text-[#c9d1d9]"
              >
                <div className="flex items-center gap-1">
                  <span>PIS SITES (%)</span>
                  <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th
                onClick={() => handleSort('gap_percent')}
                className="py-2.5 px-3 cursor-pointer hover:text-[#c9d1d9]"
              >
                <div className="flex items-center gap-1">
                  <span>GAP %</span>
                  <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th
                onClick={() => handleSort('mean_divergence')}
                className="py-2.5 px-3 cursor-pointer hover:text-[#c9d1d9]"
              >
                <div className="flex items-center gap-1">
                  <span>DIVERGENCE</span>
                  <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th
                onClick={() => handleSort('gc_percent')}
                className="py-2.5 px-3 cursor-pointer hover:text-[#c9d1d9]"
              >
                <div className="flex items-center gap-1">
                  <span>GC %</span>
                  <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th
                onClick={() => handleSort('orf_candidate_found')}
                className="py-2.5 px-3 cursor-pointer hover:text-[#c9d1d9]"
              >
                <div className="flex items-center gap-1">
                  <span>ORF CANDIDATE</span>
                  <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th
                onClick={() => handleSort('pass')}
                className="py-2.5 px-3 cursor-pointer hover:text-[#c9d1d9]"
              >
                <div className="flex items-center gap-1">
                  <span>STATUS</span>
                  <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th className="py-2.5 px-3 w-16 text-center">VIEW</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-[#1b1f27]">
            {topPadding > 0 && (
              <tr>
                <td style={{ height: `${topPadding}px` }} colSpan={12} />
              </tr>
            )}
            {visibleItems.map((item) => {
              const isSelected = item.id === selectedLocusId;
              const isChecked = selectedPaths.has(item.file_path);

              return (
                <tr
                  key={item.id}
                  onClick={() => onSelectLocus(item.id, item.file_path)}
                  style={{ height: `${ROW_HEIGHT}px` }}
                  className={`cursor-pointer transition-colors ${
                    isSelected
                      ? 'bg-blue-600/20 text-[#ffffff] font-medium border-l-2 border-blue-500'
                      : 'hover:bg-[#14171d] text-[#c9d1d9]'
                  }`}
                >
                  <td
                    className="py-2 px-3 text-center"
                    onClick={(e) => {
                      e.stopPropagation();
                      onToggleSelectPath(item.file_path);
                    }}
                  >
                    <input
                      type="checkbox"
                      checked={isChecked}
                      onChange={() => {}}
                      className="rounded bg-[#1f242e] border-[#2d3545] text-blue-500 cursor-pointer"
                    />
                  </td>
                  <td className="py-2 px-3 font-mono font-semibold text-[#dce6ff] truncate max-w-[180px]">
                    {item.id}
                  </td>
                  <td className="py-2 px-3 font-mono text-[#8b949e]">{item.num_taxa}</td>
                  <td className="py-2 px-3 font-mono">{item.length.toLocaleString()} bp</td>
                  <td className="py-2 px-3 font-mono text-cyan-400">
                    {item.variable_count ?? 0} ({Number(item.variable_percent ?? 0).toFixed(1)}%)
                  </td>
                  <td className="py-2 px-3 font-mono font-medium text-orange-400">
                    {item.pis_count} ({Number(item.pis_percent ?? 0).toFixed(1)}%)
                  </td>
                  <td className="py-2 px-3 font-mono">
                    <span
                      className={
                        item.gap_percent > 50
                          ? 'text-rose-400 font-semibold'
                          : item.gap_percent > 20
                          ? 'text-amber-400'
                          : 'text-[#8b949e]'
                      }
                    >
                      {Number(item.gap_percent ?? 0).toFixed(1)}%
                    </span>
                  </td>
                  <td className="py-2 px-3 font-mono text-[#8b949e]">
                    {(item.mean_divergence * 100).toFixed(1)}%
                  </td>
                  <td className="py-2 px-3 font-mono text-[#8b949e]">
                    {Number(item.gc_percent ?? 0).toFixed(1)}%
                  </td>
                  <td className="py-2 px-3 font-mono text-[11px]">
                    {orfEnabled && skipNonCodingOrf && shouldSkipOrfLocus(item.id, orfSearchMode) ? (
                      <span className="text-[#6e7681]">Skipped</span>
                    ) : item.orf_evaluated ? (
                      item.orf_candidate_found && item.orf_valid ? (
                        <span className="text-emerald-400">Accepted</span>
                      ) : (
                        <span className="text-rose-400">Discarded</span>
                      )
                    ) : item.orf_reference_evaluated && !item.orf_reference_matched ? (
                      <span className="text-rose-400">Discarded</span>
                    ) : (
                      <span className="text-[#6e7681]">Not assessed</span>
                    )}
                  </td>
                  <td className="py-2 px-3">
                    {isProcessing ? (
                      <span className="inline-flex items-center gap-1 text-[11px] text-cyan-300 font-medium">
                        <Loader2 className="w-3.5 h-3.5 animate-spin" />
                        <span>Processing</span>
                      </span>
                    ) : item.pass ? (
                      <span className="inline-flex items-center gap-1 text-[11px] text-emerald-400 font-medium">
                        <CheckCircle className="w-3.5 h-3.5" />
                        <span>Pass</span>
                      </span>
                    ) : (
                      <span
                        className="inline-flex items-center gap-1 text-[11px] text-rose-400 font-medium cursor-help"
                        title={item.fail_reasons.join(', ')}
                      >
                        <XCircle className="w-3.5 h-3.5" />
                        <span>Fail</span>
                      </span>
                    )}
                  </td>
                  <td className="py-2 px-3 text-center">
                    <button
                      onClick={(e) => {
                        e.stopPropagation();
                        onSelectLocus(item.id, item.file_path);
                      }}
                      className="p-1 rounded text-[#8b949e] hover:text-[#dce6ff] hover:bg-[#232833] transition-colors"
                      title="Open in MSA Viewer"
                    >
                      <Eye className="w-3.5 h-3.5" />
                    </button>
                  </td>
                </tr>
              );
            })}
            {bottomPadding > 0 && (
              <tr>
                <td style={{ height: `${bottomPadding}px` }} colSpan={12} />
              </tr>
            )}
          </tbody>
        </table>
      </div>
    </div>
  );
});

export default CatalogView;
