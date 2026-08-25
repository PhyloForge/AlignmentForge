import React, { useState, useEffect, useCallback, useMemo, useRef } from 'react';
import {
  AlignmentSummary,
  AlignmentViewResponse,
  CatalogSortField,
  ColorScheme,
  AminoAcidViewerSettings,
  DatasetOverview,
  TaxonOccupancy,
  TrimmingRecipe,
  ViewMode,
} from './types';
import {
  getAlignment,
  getPresets,
  exportAlignmentStatsCsv,
  exportFilterConfig,
  isTauri,
  loadFilterConfig,
  loadDirectoryFromFiles,
  openDirectoryDialog,
  recalculateCatalog,
  scanDirectory,
} from './tauriClient';
import { reevaluateAlignmentView, reevaluateSummaries } from './parsers/clientParser';
import { Header } from './components/Header';
import { FilterSidebar } from './components/FilterSidebar';
import { OrfAnalysisSidebar } from './components/OrfAnalysisSidebar';
import { OrfAnalysisView } from './components/OrfAnalysisView';
import { CatalogView } from './components/CatalogView';
import { MatrixHeatmap } from './components/MatrixHeatmap';
import { QcDistributions } from './components/QcDistributions';
import { MsaViewer } from './components/MsaViewer';
import { ExportModal } from './components/ExportModal';
import {
  ArrowLeft,
  FolderOpen,
  UploadCloud,
  FileCode,
  CheckCircle2,
  Sparkles,
  Loader2,
  XCircle,
} from 'lucide-react';
import iconDark from './assets/icon-dark.png';
import iconLight from './assets/icon-light.png';

function sequenceProcessingRecipeKey(recipe: TrimmingRecipe): string {
  const processingFields: Partial<TrimmingRecipe> = { ...recipe };
  delete processingFields.name;
  delete processingFields.description;
  delete processingFields.assess_alignment;
  delete processingFields.min_taxa;
  delete processingFields.min_taxa_occupancy_percent;
  delete processingFields.min_length;
  delete processingFields.max_gap_percent;
  delete processingFields.min_pis_count;
  delete processingFields.min_pis_percent;
  delete processingFields.min_variable_count;
  delete processingFields.min_variable_percent;
  delete processingFields.fail_if_no_orf;
  return JSON.stringify(processingFields);
}

function assessmentRecipeKey(recipe: TrimmingRecipe): string {
  return JSON.stringify({
    assess_alignment: recipe.assess_alignment,
    min_taxa: recipe.min_taxa,
    min_taxa_occupancy_percent: recipe.min_taxa_occupancy_percent,
    min_length: recipe.min_length,
    max_gap_percent: recipe.max_gap_percent,
    min_pis_count: recipe.min_pis_count,
    min_pis_percent: recipe.min_pis_percent,
    min_variable_count: recipe.min_variable_count,
    min_variable_percent: recipe.min_variable_percent,
    fail_if_no_orf: recipe.fail_if_no_orf,
  });
}

export const App: React.FC = () => {
  const [currentPath, setCurrentPath] = useState<string | null>(null);
  const [summaries, setSummaries] = useState<AlignmentSummary[]>([]);
  const [overview, setOverview] = useState<DatasetOverview>({
    total_alignments: 0,
    passed_alignments: 0,
    discarded_alignments: 0,
    total_unique_taxa: 0,
    mean_taxa: 0,
    mean_length: 0,
    mean_gap_percent: 0,
    mean_pis: 0,
    total_matrix_basepairs: 0,
  });
  const [occupancy, setOccupancy] = useState<TaxonOccupancy[]>([]);

  const [activeView, setActiveView] = useState<ViewMode>('catalog');
  const [msaSidebarContext, setMsaSidebarContext] = useState<'catalog' | 'orf'>('catalog');
  const [selectedLocusId, setSelectedLocusId] = useState<string | null>(null);
  const [selectedFilePath, setSelectedFilePath] = useState<string | null>(null);
  const [viewData, setViewData] = useState<AlignmentViewResponse | null>(null);
  const [alignmentLoadError, setAlignmentLoadError] = useState<string | null>(null);

  const [recipes, setRecipes] = useState<TrimmingRecipe[]>([]);
  const [activeRecipe, setActiveRecipe] = useState<TrimmingRecipe>({
    name: 'AlignmentForge Default',
    description: 'Standard balanced phylogenomic filtering pipeline',
    replace_n_with_gap: true,
    ambiguity_strategy: 'keep',
    remove_gap_only_columns: true,
    trim_similarity: true,
    similarity_threshold: 0.4,
    trim_hmm: false,
    hmm_min_posterior: 0.45,
    hmm_min_segment_length: 8,
      hmm_min_island_length: 20,
    trim_segments: false,
    segment_window_size: 100,
    segment_threshold: 0.45,
    enable_orf: false,
    auto_shift_frame: true,
    auto_flip_reverse: true,
    stop_codon_action: 'removesample',
  macse_trim_terminal: true,
  macse_max_internal_sample: 3,
  macse_max_internal_locus: 10,
  max_stop_codons_sample: 2,
  max_stop_codons_locus: 5,
    genetic_code: 'standard',
    orf_search_mode: 'continuouscds',
    orf_min_shared_support_percent: 90.0,
    orf_min_segment_aa: 35,
    orf_min_coding_score: 40.0,
    exclude_uce: true,
    fail_if_no_orf: true,
    orf_use_references: false,
    orf_reference_sequences: {},
    trim_external: true,
    min_external_percent: 50.0,
    codon_preserving: false,
    trim_columns: false,
    min_column_gap_percent: 60.0,
    count_n_as_gap: true,
    enable_statistical_columns: false,
    stat_col_method: 'trimalsimilarity',
    stat_col_similarity_threshold: 0.35,
    stat_col_window_size: 3,
    stat_col_heuristic: 'custom',
    stat_col_min_block_length: 5,
    stat_col_max_nonconserved: 4,
    stat_col_gap_treatment: 'half',
    stat_col_entropy_threshold: 1.5,
    trim_coverage: true,
    min_coverage_bp: 60,
    min_coverage_percent: 50.0,
    relative_width: 'sample',
    min_sample_locus_occupancy_percent: 0.0,
    assess_alignment: true,
    min_taxa: 4,
    min_taxa_occupancy_percent: 50.0,
    min_length: 100,
    max_gap_percent: 50.0,
    min_pis_count: 0,
    min_pis_percent: 0.0,
    min_variable_count: 0,
    min_variable_percent: 0.0,
  });

  const [catalogSearchTerm, setCatalogSearchTerm] = useState<string>('');
  const [catalogStatusFilter, setCatalogStatusFilter] = useState<'all' | 'pass' | 'fail'>('all');
  const [orfStatusFilter, setOrfStatusFilter] = useState<'all' | 'accepted' | 'discarded'>('all');
  const [catalogSortField, setCatalogSortField] = useState<CatalogSortField>('id');
  const [catalogSortAsc, setCatalogSortAsc] = useState<boolean>(true);

  const [selectedPaths, setSelectedPaths] = useState<Set<string>>(new Set());
  const [colorScheme, setColorScheme] = useState<ColorScheme>('nucleotide');
  const [aminoAcidViewerSettings, setAminoAcidViewerSettings] =
    useState<AminoAcidViewerSettings>({
      enabled: false,
      hideTerminalStops: false,
      colorScheme: 'chemistry',
      dimConsensusMatches: false,
    });
  const [showDiffOverlay, setShowDiffOverlay] = useState<boolean>(true);
  const [exportModalMode, setExportModalMode] = useState<'batch' | 'concatenate' | null>(null);
  const [isDarkTheme, setIsDarkTheme] = useState<boolean>(true);
    useEffect(() => {
    if (isDarkTheme) {
      document.documentElement.classList.add('dark');
    } else {
      document.documentElement.classList.remove('dark');
    }
  }, [isDarkTheme]);
  const [isDragging, setIsDragging] = useState<boolean>(false);
  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [isCatalogProcessing, setIsCatalogProcessing] = useState<boolean>(false);
  const [catalogProcessingPercent, setCatalogProcessingPercent] = useState<number>(0);
  const [referenceRevision, setReferenceRevision] = useState(0);
  const [loadingProgress, setLoadingProgress] = useState<{
    current: number;
    total: number;
    percent: number;
    fileName: string;
  } | null>(null);

  const fileInputRef = useRef<HTMLInputElement>(null);
  const rawSummariesRef = useRef<AlignmentSummary[]>([]);
  const catalogRequestIdRef = useRef(0);
  const activeRecipeRef = useRef(activeRecipe);
  activeRecipeRef.current = activeRecipe;

  const processingRecipeKey = sequenceProcessingRecipeKey(activeRecipe);
  const processingRecipe = useMemo(() => activeRecipe, [processingRecipeKey]);
  const gatingRecipeKey = assessmentRecipeKey(activeRecipe);

  // Initialize presets & Tauri progress event listener
  useEffect(() => {
    getPresets().then((presetList) => {
      setRecipes(presetList);
      if (presetList.length > 0) {
        setActiveRecipe(presetList[0]);
      }
    });

    // Listen to native Tauri scan_progress events
    let unlisten: (() => void) | undefined;
    if (isTauri) {
      import('@tauri-apps/api/event').then(({ listen }) => {
        listen<any>('scan_progress', (event) => {
          const payload = event.payload;
          setLoadingProgress({
            current: payload.current,
            total: payload.total,
            percent: payload.percent,
            fileName: payload.file_name,
          });
        }).then((unsub) => {
          unlisten = unsub;
        });
      });
    }

    return () => {
      if (unlisten) unlisten();
    };
  }, []);

  const handleOpenDirectory = async () => {
    if (isTauri) {
      const targetDir = await openDirectoryDialog();
      if (targetDir) {
        handleLoadPath(targetDir);
      }
    } else {
      fileInputRef.current?.click();
    }
  };

  const handleLoadPath = async (dirPath: string) => {
    setIsLoading(true);
    setLoadingProgress({ current: 0, total: 0, percent: 0, fileName: 'Discovering alignment files…' });
    setCurrentPath(dirPath);
    try {
      const res = await scanDirectory(dirPath);
      rawSummariesRef.current = res.summaries;
      setSummaries(res.summaries);
      setOverview(res.overview);
      setOccupancy(res.occupancy);

      if (res.summaries.length > 0) {
        const first = res.summaries[0];
        setSelectedLocusId(first.id);
        setSelectedFilePath(first.file_path);
      }
    } catch (e: any) {
      console.error('Failed to scan directory:', e);
      alert(`Could not load directory: ${e?.message || e}`);
    } finally {
      setIsLoading(false);
      setLoadingProgress(null);
    }
  };

  const handleFilesSelected = async (files: FileList | File[]) => {
    if (!files || files.length === 0) return;
    setIsLoading(true);
    setLoadingProgress({ current: 0, total: files.length, percent: 0, fileName: 'Reading files…' });
    try {
      const { dirName, scanResponse } = await loadDirectoryFromFiles(
        files,
        activeRecipe,
        (current, total, fileName) => {
          const percent = total > 0 ? (current / total) * 100 : 0;
          setLoadingProgress({ current, total, percent, fileName });
        }
      );
      setCurrentPath(dirName);
      rawSummariesRef.current = scanResponse.summaries;
      setSummaries(scanResponse.summaries);
      setOverview(scanResponse.overview);
      setOccupancy(scanResponse.occupancy);

      if (scanResponse.summaries.length > 0) {
        const first = scanResponse.summaries[0];
        setSelectedLocusId(first.id);
        setSelectedFilePath(first.file_path);
      }
    } catch (e: any) {
      console.error('Failed to load files:', e);
      alert(`Could not load alignments: ${e?.message || e}`);
    } finally {
      setIsLoading(false);
      setLoadingProgress(null);
    }
  };

  const handleSelectLocus = useCallback(
    (id: string, filePath: string) => {
      if (selectedFilePath === filePath) return;
      setSelectedLocusId(id);
      setSelectedFilePath(filePath);
    },
    [selectedFilePath]
  );

  // Alignment previews are viewer-only. Catalog rows are updated exclusively by the
  // dataset-wide recalculation below, so opening an alignment is never required to
  // compute (or able to overwrite) its table status.
  useEffect(() => {
    if (activeView !== 'msa' || !selectedFilePath) {
      setViewData(null);
      setAlignmentLoadError(null);
      return;
    }

    setViewData(null);
    setAlignmentLoadError(null);
    let isCancelled = false;
    const totalUniqueTaxa = overview.total_unique_taxa || occupancy.length || 0;
    getAlignment(selectedFilePath, processingRecipe, totalUniqueTaxa)
      .then((data) => {
        if (!isCancelled) {
          setViewData(
            reevaluateAlignmentView(data, activeRecipeRef.current, totalUniqueTaxa)
          );
        }
      })
      .catch((err) => {
        if (!isCancelled) {
          console.error('Failed to load alignment:', err);
          setAlignmentLoadError(err instanceof Error ? err.message : String(err));
        }
      });

    return () => {
      isCancelled = true;
    };
  }, [activeView, selectedFilePath, processingRecipe, overview.total_unique_taxa, occupancy.length]);

  // Instant in-memory re-gating whenever recipe assessment thresholds change.
  // Uses the latest completed pipeline metrics for immediate slider feedback.
  useEffect(() => {
    if (rawSummariesRef.current.length === 0) return;
    const totalUniqueTaxa = overview.total_unique_taxa || occupancy.length || 0;
    const { summaries: updated, overview: updatedOverview } = reevaluateSummaries(
      rawSummariesRef.current,
      activeRecipe,
      totalUniqueTaxa
    );
    setSummaries(updated);
    setOverview(updatedOverview);
  }, [activeRecipe]);

  // Assessment controls only change pass/fail. Keep an open preview synchronized
  // without re-running any sequence transformations.
  useEffect(() => {
    const totalUniqueTaxa = overview.total_unique_taxa || occupancy.length || 0;
    setViewData((current) =>
      current ? reevaluateAlignmentView(current, activeRecipe, totalUniqueTaxa) : current
    );
  }, [gatingRecipeKey, overview.total_unique_taxa, occupancy.length]);

  // Run the full trimming pipeline across every alignment whenever the recipe changes.
  // A request generation prevents slower, stale calculations from overwriting the
  // newest recipe. This job runs regardless of which view is open.
  useEffect(() => {
    const requestId = ++catalogRequestIdRef.current;
    if (isLoading || !currentPath || rawSummariesRef.current.length === 0) {
      setIsCatalogProcessing(false);
      setCatalogProcessingPercent(0);
      return;
    }

    const paths = rawSummariesRef.current.map((s) => s.file_path);
    const totalUniqueTaxa = overview.total_unique_taxa || occupancy.length || 0;
    setIsCatalogProcessing(true);
    setCatalogProcessingPercent(0);

    const timer = window.setTimeout(() => {
      recalculateCatalog(paths, processingRecipe, totalUniqueTaxa, (progress) => {
        if (catalogRequestIdRef.current === requestId) {
          const nextPercent = Math.max(0, Math.min(100, progress.percent));
          setCatalogProcessingPercent((currentPercent) => Math.max(currentPercent, nextPercent));
        }
      })
        .then((res) => {
          if (catalogRequestIdRef.current !== requestId || res.summaries.length === 0) return;
          rawSummariesRef.current = res.summaries;
          const latestGating = reevaluateSummaries(
            res.summaries,
            activeRecipeRef.current,
            totalUniqueTaxa
          );
          setSummaries(latestGating.summaries);
          setOverview(latestGating.overview);
        })
        .catch((err) => {
          if (catalogRequestIdRef.current === requestId) {
            console.error('Failed to recalculate catalog:', err);
          }
        })
        .finally(() => {
          if (catalogRequestIdRef.current === requestId) {
            setIsCatalogProcessing(false);
            setCatalogProcessingPercent(0);
          }
        });
    }, 150);

    return () => {
      window.clearTimeout(timer);
    };
  }, [processingRecipe, referenceRevision, currentPath, isLoading, overview.total_unique_taxa, occupancy.length]);

  /*
   * Keep this state reset close to loading so a dataset switch cannot leave a stale
   * processing badge behind while the new directory is being indexed.
   */
  useEffect(() => {
    if (isLoading) {
      setViewData(null);
    }
  }, [isLoading]);

  const handleToggleSelectPath = useCallback((path: string) => {
    setSelectedPaths((prev) => {
      const next = new Set(prev);
      if (next.has(path)) next.delete(path);
      else next.add(path);
      return next;
    });
  }, []);

  const handleSelectAllPaths = useCallback(() => {
    setSelectedPaths(new Set(summaries.map((s) => s.file_path)));
  }, [summaries]);

  const handleClearSelectedPaths = useCallback(() => {
    setSelectedPaths(new Set());
  }, []);

  const handleSelectLocusAndSwitchView = useCallback(
    (id: string, filePath: string, sidebarContext: 'catalog' | 'orf') => {
      setViewData(null);
      setAlignmentLoadError(null);
      setMsaSidebarContext(sidebarContext);
      handleSelectLocus(id, filePath);
      setActiveView('msa');
    },
    [handleSelectLocus]
  );

  const handleSelectHeaderView = useCallback(
    (view: ViewMode) => {
      if (view === 'msa' && activeView !== 'msa') {
        setMsaSidebarContext(activeView === 'orf' ? 'orf' : 'catalog');
      }
      setActiveView(view);
    },
    [activeView]
  );

  // Drag & drop handlers
  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(true);
  };

  const handleDragLeave = () => {
    setIsDragging(false);
  };

  const handleDrop = async (e: React.DragEvent) => {
    e.preventDefault();
    setIsDragging(false);
    if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
      handleFilesSelected(e.dataTransfer.files);
    }
  };

  return (
    <div
      className={`h-screen w-screen flex flex-col overflow-hidden ${isDarkTheme ? 'dark' : ''}`}
      onDragOver={handleDragOver}
      onDragLeave={handleDragLeave}
      onDrop={handleDrop}
    >
      {/* Hidden File / Folder Input */}
      <input
        ref={fileInputRef}
        type="file"
        // @ts-ignore
        webkitdirectory="true"
        directory="true"
        multiple
        className="hidden"
        onChange={(e) => {
          if (e.target.files) {
            handleFilesSelected(e.target.files);
          }
        }}
      />

      {/* Drag Overlay */}
      {isDragging && (
        <div className="fixed inset-0 z-50 bg-blue-600/30 backdrop-blur-sm border-4 border-dashed border-blue-400 flex items-center justify-center pointer-events-none">
          <div className="bg-[#14171d] border border-blue-500 rounded-2xl p-8 text-center shadow-2xl space-y-3">
            <UploadCloud className="w-16 h-16 text-blue-400 mx-auto animate-bounce" />
            <h2 className="text-xl font-bold text-[#dce6ff]">Drop Alignment Folder to Open</h2>
            <p className="text-xs text-[#8b949e]">Accepts PHYLIP (.phy), FASTA (.fa), NEXUS (.nex)</p>
          </div>
        </div>
      )}

      {/* Top Main Navigation Header */}
      <Header
        currentPath={currentPath}
        totalAlignments={overview.total_alignments}
        passedAlignments={overview.passed_alignments}
        activeView={activeView}
        onSelectView={handleSelectHeaderView}
        onOpenDirectory={handleOpenDirectory}
        onOpenExportModal={setExportModalMode}
        isDarkTheme={isDarkTheme}
        onToggleTheme={() => setIsDarkTheme((t) => !t)}
      />

      {/* Main Content Area */}
      {summaries.length === 0 && !isLoading ? (
        /* Empty State Screen */
        <div className="flex-1 flex flex-col items-center justify-center bg-[#0e1014] text-[#c9d1d9] p-8 select-none">
          <div className="max-w-xl w-full bg-[#14171d] border border-[#232833] rounded-2xl p-10 text-center shadow-2xl space-y-6">
            <img
              src={isDarkTheme ? iconLight : iconDark}
              alt="AlignmentForge"
              className="w-32 h-32 rounded-3xl shadow-2xl border border-white/10 mx-auto object-cover"
            />

            <div className="space-y-2">
              <h1 className="text-2xl font-bold text-[#dce6ff] tracking-tight">
                Welcome to AlignmentForge
              </h1>
              <p className="text-xs text-[#8b949e] leading-relaxed">
                High-performance visualization, quality control, and recipe-driven trimming for thousands of phylogenomic alignments.
              </p>
            </div>

            <div className="pt-2 flex flex-col items-center gap-3">
              <button
                onClick={handleOpenDirectory}
                className="px-6 py-3 bg-blue-600 hover:bg-blue-500 text-white rounded-xl font-semibold text-sm shadow-lg shadow-blue-500/20 hover:shadow-blue-500/40 transition-all flex items-center gap-2"
              >
                <FolderOpen className="w-4 h-4" />
                <span>Open Alignments Directory</span>
              </button>
              <span className="text-[11px] text-[#8b949e]">
                or drag and drop your alignment folder here
              </span>
            </div>

            <div className="pt-6 border-t border-[#232833] grid grid-cols-3 gap-3 text-left">
              <div className="p-3 bg-[#171b22] border border-[#232833] rounded-lg space-y-1">
                <div className="text-[11px] font-semibold text-[#dce6ff] flex items-center gap-1.5">
                  <FileCode className="w-3.5 h-3.5 text-cyan-400" />
                  <span>Universal Formats</span>
                </div>
                <div className="text-[10px] text-[#8b949e]">
                  Relaxed PHYLIP, FASTA, NEXUS sequential & interleaved
                </div>
              </div>

              <div className="p-3 bg-[#171b22] border border-[#232833] rounded-lg space-y-1">
                <div className="text-[11px] font-semibold text-[#dce6ff] flex items-center gap-1.5">
                  <Sparkles className="w-3.5 h-3.5 text-purple-400" />
                  <span>Real-Time Trimming</span>
                </div>
                <div className="text-[10px] text-[#8b949e]">
                  Instant live visual previews across 8 quality and trimming filters
                </div>
              </div>

              <div className="p-3 bg-[#171b22] border border-[#232833] rounded-lg space-y-1">
                <div className="text-[11px] font-semibold text-[#dce6ff] flex items-center gap-1.5">
                  <CheckCircle2 className="w-3.5 h-3.5 text-emerald-400" />
                  <span>Non-Destructive</span>
                </div>
                <div className="text-[10px] text-[#8b949e]">
                  Original files are never modified; batch export to folder
                </div>
              </div>
            </div>
          </div>
        </div>
      ) : (
        /* Main Workspace Body */
        <div className="flex-1 flex overflow-hidden relative">
          {/* Left Filter Sidebar */}
          {activeView === 'orf' || (activeView === 'msa' && msaSidebarContext === 'orf') ? (
            <OrfAnalysisSidebar
              recipe={activeRecipe}
              onChangeRecipe={setActiveRecipe}
              onResetRecipe={() => {
                const found = recipes.find((r) => r.name === activeRecipe.name);
                if (found) setActiveRecipe(found);
              }}
              onReferencesChanged={() => setReferenceRevision((revision) => revision + 1)}
              aminoAcidViewerSettings={aminoAcidViewerSettings}
              onChangeAminoAcidViewerSettings={setAminoAcidViewerSettings}
            />
          ) : (
            <FilterSidebar
              recipes={recipes}
              recipe={activeRecipe}
              onChangeRecipe={setActiveRecipe}
              onResetRecipe={() => {
                const found = recipes.find((r) => r.name === activeRecipe.name);
                if (found) setActiveRecipe(found);
              }}
              onExportFilters={() => exportFilterConfig(activeRecipe)}
              onLoadFilters={async () => {
                const loaded = await loadFilterConfig();
                if (!loaded) return null;

                setRecipes((currentRecipes) => {
                  const alreadyListed = currentRecipes.some(
                    (candidate) => candidate.name === loaded.recipe.name
                  );
                  if (!alreadyListed) return [...currentRecipes, loaded.recipe];
                  return currentRecipes.map((candidate) =>
                    candidate.name === loaded.recipe.name ? loaded.recipe : candidate
                  );
                });
                setActiveRecipe(loaded.recipe);
                return loaded.filePath;
              }}
              passedCount={overview.passed_alignments}
              totalCount={overview.total_alignments}
            />
          )}

          {/* Center Viewport Switcher */}
          <main className="flex-1 flex flex-col overflow-hidden">
            {activeView === 'catalog' && (
              <CatalogView
                summaries={summaries}
                isProcessing={isCatalogProcessing}
                processingPercent={catalogProcessingPercent}
                selectedLocusId={selectedLocusId}
                onSelectLocus={(id: string, filePath: string) =>
                  handleSelectLocusAndSwitchView(id, filePath, 'catalog')
                }
                selectedPaths={selectedPaths}
                onToggleSelectPath={handleToggleSelectPath}
                onSelectAllPaths={handleSelectAllPaths}
                onClearSelectedPaths={handleClearSelectedPaths}
                searchTerm={catalogSearchTerm}
                onSearchTermChange={setCatalogSearchTerm}
                statusFilter={catalogStatusFilter}
                onStatusFilterChange={setCatalogStatusFilter}
                sortField={catalogSortField}
                sortAsc={catalogSortAsc}
                onSortChange={(field, asc) => {
                  setCatalogSortField(field);
                  setCatalogSortAsc(asc);
                }}
                orfEnabled={activeRecipe.enable_orf}
                orfSearchMode={activeRecipe.orf_search_mode}
                skipNonCodingOrf={activeRecipe.exclude_uce}
              />
            )}

            {activeView === 'matrix' && (
              <MatrixHeatmap
                occupancy={occupancy}
                summaries={summaries}
                onSelectLocus={(id: string, filePath: string) =>
                  handleSelectLocusAndSwitchView(id, filePath, 'catalog')
                }
              />
            )}

            {activeView === 'orf' && (
              <OrfAnalysisView
                selectedPaths={selectedPaths}
                onSelectPath={(path: string, selected: boolean) => {
                  const next = new Set(selectedPaths);
                  if (selected) next.add(path);
                  else next.delete(path);
                  setSelectedPaths(next);
                }}
                onSelectAllPaths={() => {
                  const next = new Set(summaries.map((s) => s.file_path));
                  setSelectedPaths(next);
                }}
                onClearSelectedPaths={() => setSelectedPaths(new Set())}
                summaries={summaries}
                enabled={activeRecipe.enable_orf}
                isProcessing={isCatalogProcessing}
                processingPercent={catalogProcessingPercent}
                searchTerm={catalogSearchTerm}
                onSearchTermChange={setCatalogSearchTerm}
                statusFilter={orfStatusFilter}
                onStatusFilterChange={setOrfStatusFilter}
                orfSearchMode={activeRecipe.orf_search_mode}
                skipNonCodingOrf={activeRecipe.exclude_uce}
                onSelectLocus={(id: string, filePath: string) =>
                  handleSelectLocusAndSwitchView(id, filePath, 'orf')
                }
              />
            )}

            {activeView === 'qc' && (
              <QcDistributions
                overview={overview}
                summaries={summaries}
                occupancy={occupancy}
                orfEnabled={activeRecipe.enable_orf}
                onExportStats={() => exportAlignmentStatsCsv(summaries)}
              />
            )}

            {activeView === 'msa' && (
              viewData ? (
                <MsaViewer
                  viewData={viewData}
                  colorScheme={colorScheme}
                  onChangeColorScheme={setColorScheme}
                  showDiffOverlay={showDiffOverlay}
                  onToggleDiffOverlay={() => setShowDiffOverlay((v) => !v)}
                  geneticCode={activeRecipe.genetic_code}
                  aminoAcidViewerSettings={aminoAcidViewerSettings}
                  onChangeAminoAcidViewerSettings={setAminoAcidViewerSettings}
                />
              ) : (
                <div className="flex-1 flex items-center justify-center bg-[#0e1014] text-[#8b949e]">
                  <div className="flex flex-col items-center gap-3 text-center px-6">
                    {alignmentLoadError ? (
                      <>
                        <XCircle className="w-9 h-9 text-rose-400" />
                        <div>
                          <p className="text-sm font-medium text-rose-300">Could not process this alignment</p>
                          <p className="text-xs text-[#8b949e] mt-1 max-w-lg">{alignmentLoadError}</p>
                        </div>
                      </>
                    ) : (
                      <>
                        <Loader2 className="w-9 h-9 text-emerald-400 animate-spin" />
                        <div>
                          <p className="text-sm font-medium text-[#dce6ff]">Processing alignment…</p>
                          <p className="text-xs text-[#8b949e] mt-1">
                            Applying the active recipe without blocking catalog processing.
                          </p>
                        </div>
                      </>
                    )}
                    <button
                      onClick={() => setActiveView(msaSidebarContext)}
                      className="mt-2 inline-flex items-center gap-1.5 px-3 py-1.5 rounded-md bg-[#1b2029] hover:bg-[#232a36] text-[#c9d1d9] border border-[#2d3545] text-xs font-medium transition-colors"
                    >
                      <ArrowLeft className="w-3.5 h-3.5" />
                      Back to {msaSidebarContext === 'orf' ? 'ORF Analysis' : 'Catalog'}
                    </button>
                  </div>
                </div>
              )
            )}
          </main>
        </div>
      )}

      {/* Modern Live Loading Progress Bar Modal */}
      {isLoading && (
        <div className="fixed inset-0 z-50 bg-black/75 backdrop-blur-md flex items-center justify-center p-6 select-none">
          <div className="bg-[#14171d] border border-[#2d3545] rounded-2xl p-8 max-w-md w-full shadow-2xl space-y-5 text-center">
            <div className="w-12 h-12 rounded-xl bg-blue-500/10 border border-blue-500/20 text-blue-400 flex items-center justify-center mx-auto">
              <Loader2 className="w-6 h-6 animate-spin" />
            </div>

            <div className="space-y-1">
              <h3 className="text-base font-semibold text-[#dce6ff]">
                Indexing Alignment Dataset
              </h3>
              <p className="text-xs text-[#8b949e]">
                Parsing alignments, calculating occupancy, and computing quality metrics in parallel…
              </p>
            </div>

            {/* Progress Track & Fill */}
            <div className="space-y-2">
              <div className="w-full bg-[#1b2029] h-3 rounded-full overflow-hidden border border-[#232833] p-0.5">
                <div
                  className="bg-gradient-to-r from-blue-600 via-cyan-500 to-emerald-400 h-full rounded-full transition-all duration-150 ease-out shadow-sm"
                  style={{
                    width: `${loadingProgress && loadingProgress.total > 0 ? Math.max(3, loadingProgress.percent) : 100}%`,
                  }}
                />
              </div>

              {/* Progress Counters & Current File */}
              <div className="flex items-center justify-between text-[11px] font-mono text-[#8b949e]">
                <span>
                  {loadingProgress && loadingProgress.total > 0
                    ? `${loadingProgress.current.toLocaleString()} / ${loadingProgress.total.toLocaleString()} files`
                    : 'Discovering files in directory…'}
                </span>
                <span className="text-cyan-400 font-semibold">
                  {loadingProgress && loadingProgress.total > 0 ? `${loadingProgress.percent.toFixed(1)}%` : ''}
                </span>
              </div>
            </div>

            {loadingProgress?.fileName && (
              <div className="text-[10px] font-mono text-[#8b949e] truncate bg-[#0e1014] px-3 py-1.5 rounded border border-[#232833]">
                {loadingProgress.fileName}
              </div>
            )}
          </div>
        </div>
      )}

      {/* Export / Supermatrix Modal */}
      {exportModalMode && (
        <ExportModal
          mode={exportModalMode}
          onClose={() => setExportModalMode(null)}
          selectedPaths={Array.from(selectedPaths)}
          allPaths={summaries.map((s) => s.file_path)}
          recipe={activeRecipe}
        />
      )}
    </div>
  );
};

export default App;
