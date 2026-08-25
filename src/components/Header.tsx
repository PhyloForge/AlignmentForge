import React from 'react';
import {
  FolderOpen,
  Table,
  Grid,
  BarChart2,
  Eye,
  Download,
  Layers,
  Sparkles,
  Sun,
  Moon,
  Dna,
} from 'lucide-react';
import { ViewMode } from '../types';
import iconDark from '../assets/icon-dark.png';
import iconLight from '../assets/icon-light.png';

interface HeaderProps {
  currentPath: string | null;
  totalAlignments: number;
  passedAlignments: number;
  activeView: ViewMode;
  onSelectView: (view: ViewMode) => void;
  onOpenDirectory: () => void;
  onOpenExportModal: (mode: 'batch' | 'concatenate' | 'group') => void;
  isDarkTheme: boolean;
  onToggleTheme: () => void;
}

export const Header: React.FC<HeaderProps> = ({
  currentPath,
  totalAlignments,
  passedAlignments,
  activeView,
  onSelectView,
  onOpenDirectory,
  onOpenExportModal,
  isDarkTheme,
  onToggleTheme,
}) => {
  return (
    <header className="flex flex-wrap items-center justify-between gap-4 px-4 py-2.5 bg-[#14171d] border-b border-[#232833] select-none">
      {/* Brand & Folder Selector */}
      <div className="flex items-center gap-3">
        <div className="flex items-center gap-2.5 font-bold text-[15px] tracking-wide text-[#dce6ff]">
          <img
            src={isDarkTheme ? iconLight : iconDark}
            alt="AlignmentForge"
            className="w-8 h-8 rounded-[7px] shadow-md border border-white/10 object-cover"
          />
          <span>AlignmentForge</span>
          <span className="text-[11px] font-normal text-[#8b949e] px-1.5 py-0.5 rounded bg-[#1f242e] border border-[#2d3545]">
            v0.1.0
          </span>
        </div>

        <div className="h-4 w-[1px] bg-[#232833] mx-1" />

        <button
          onClick={onOpenDirectory}
          className="flex items-center gap-2 px-3 py-1.5 rounded-md bg-[#1b2029] hover:bg-[#232a36] text-[#c9d1d9] border border-[#2d3545] text-xs font-medium transition-colors shadow-sm"
          title="Open Alignments Folder"
        >
          <FolderOpen className="w-3.5 h-3.5 text-blue-400" />
          <span>{currentPath ? currentPath.split('/').pop() || currentPath : 'Open Directory…'}</span>
        </button>
      </div>

      {/* View Switcher Tabs */}
      <div className="flex items-center p-0.5 bg-[#0e1014] rounded-lg border border-[#232833]">
        <button
          onClick={() => onSelectView('catalog')}
          className={`flex items-center gap-1.5 px-3 py-1.5 rounded-md text-xs font-medium transition-all ${
            activeView === 'catalog'
              ? 'bg-[#232a36] text-[#dce6ff] shadow-sm font-semibold'
              : 'text-[#8b949e] hover:text-[#c9d1d9]'
          }`}
        >
          <Table className="w-3.5 h-3.5 text-indigo-400" />
          <span>Catalog</span>
        </button>

        <button
          onClick={() => onSelectView('orf')}
          className={`flex items-center gap-1.5 px-3 py-1.5 rounded-md text-xs font-medium transition-all ${
            activeView === 'orf'
              ? 'bg-[#232a36] text-[#dce6ff] shadow-sm font-semibold'
              : 'text-[#8b949e] hover:text-[#c9d1d9]'
          }`}
        >
          <Dna className="w-3.5 h-3.5 text-pink-400" />
          <span>ORF Analysis</span>
        </button>

        <button
          onClick={() => onSelectView('msa')}
          className={`flex items-center gap-1.5 px-3 py-1.5 rounded-md text-xs font-medium transition-all ${
            activeView === 'msa'
              ? 'bg-[#232a36] text-[#dce6ff] shadow-sm font-semibold'
              : 'text-[#8b949e] hover:text-[#c9d1d9]'
          }`}
        >
          <Eye className="w-3.5 h-3.5 text-emerald-400" />
          <span>MSA Viewer</span>
        </button>

        <button
          onClick={() => onSelectView('matrix')}
          className={`flex items-center gap-1.5 px-3 py-1.5 rounded-md text-xs font-medium transition-all ${
            activeView === 'matrix'
              ? 'bg-[#232a36] text-[#dce6ff] shadow-sm font-semibold'
              : 'text-[#8b949e] hover:text-[#c9d1d9]'
          }`}
        >
          <Grid className="w-3.5 h-3.5 text-cyan-400" />
          <span>Matrix</span>
        </button>

        <button
          onClick={() => onSelectView('qc')}
          className={`flex items-center gap-1.5 px-3 py-1.5 rounded-md text-xs font-medium transition-all ${
            activeView === 'qc'
              ? 'bg-[#232a36] text-[#dce6ff] shadow-sm font-semibold'
              : 'text-[#8b949e] hover:text-[#c9d1d9]'
          }`}
        >
          <BarChart2 className="w-3.5 h-3.5 text-amber-400" />
          <span>QC Stats</span>
        </button>
      </div>

      {/* Batch Export & Theme */}
      <div className="flex items-center gap-2">
        {/* Batch Export */}
        <button
          onClick={() => onOpenExportModal('batch')}
          disabled={totalAlignments === 0}
          className="flex items-center gap-1.5 px-2.5 py-1.5 rounded-md bg-[#1d273b] hover:bg-[#25334d] text-blue-300 border border-blue-500/30 text-xs font-medium transition-colors shadow-sm disabled:opacity-40 disabled:cursor-not-allowed"
          title="Batch Export Trimmed Alignments"
        >
          <Download className="w-3.5 h-3.5" />
          <span>Export</span>
        </button>

        {/* Concatenate Supermatrix */}
        <button
          onClick={() => onOpenExportModal('concatenate')}
          disabled={totalAlignments === 0}
          className="flex items-center gap-1.5 px-2.5 py-1.5 rounded-md bg-[#192b23] hover:bg-[#203a2e] text-emerald-300 border border-emerald-500/30 text-xs font-medium transition-colors shadow-sm disabled:opacity-40 disabled:cursor-not-allowed"
          title="Concatenate into Partitioned Supermatrix"
        >
          <Layers className="w-3.5 h-3.5" />
          <span>Supermatrix</span>
        </button>

        {/* Gene Batch */}
        <button
          onClick={() => onOpenExportModal('group')}
          disabled={totalAlignments === 0}
          className="flex items-center gap-1.5 px-2.5 py-1.5 rounded-md bg-[#2d1b38] hover:bg-[#3f254f] text-fuchsia-300 border border-fuchsia-500/30 text-xs font-medium transition-colors shadow-sm disabled:opacity-40 disabled:cursor-not-allowed"
          title="Concatenate Exons by Gene"
        >
          <Sparkles className="w-3.5 h-3.5" />
          <span>Gene Batch</span>
        </button>

        {/* Theme toggle */}
        <button
          onClick={onToggleTheme}
          className="p-1.5 rounded-md bg-[#171b22] hover:bg-[#232833] text-[#8b949e] hover:text-[#c9d1d9] border border-[#232833] transition-colors"
          title="Toggle Theme"
        >
          {isDarkTheme ? <Sun className="w-3.5 h-3.5" /> : <Moon className="w-3.5 h-3.5" />}
        </button>
      </div>
    </header>
  );
};
