import React, { useState, useEffect, useCallback, useMemo, useRef } from 'react';
import {
  AlignmentSummary,
  AlignmentViewResponse,
  CatalogSortField,
  ColorScheme,
  AminoAcidViewerSettings,
  DatasetOverview,
  ParseFailure,
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
import { recipeWithoutOrfAnalysis } from './summaries';
import { DEFAULT_RECIPE } from './defaultRecipe';
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
  AlertTriangle,
  X,
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
  });
}

function applyCatalogAssessmentToView(
  view: AlignmentViewResponse,
  summary: AlignmentSummary | undefined
): AlignmentViewResponse {
  if (!summary) return view;
  return {
    ...view,
    diff: {
      ...view.diff,
      pass: summary.pass,
      fail_reasons: summary.fail_reasons,
    },
  };
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
  const [parseFailures, setParseFailures] = useState<ParseFailure[]>([]);
  /** Per-sample retention, loaded only while the QC view needs it. */
  const [retentionSummaries, setRetentionSummaries] = useState<AlignmentSummary[] | null>(null);
  const [showParseFailures, setShowParseFailures] = useState<boolean>(false);
  const [loadError, setLoadError] = useState<string | null>(null);
  const [catalogError, setCatalogError] = useState<string | null>(null);

  const [recipes, setRecipes] = useState<TrimmingRecipe[]>([]);
  const [activeRecipe, setActiveRecipe] = useState<TrimmingRecipe>(DEFAULT_RECIPE);

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
  // Held here so changing the recipe, which reloads the alignment view, does
  // not clear the highlighted sample.
  const [selectedTaxon, setSelectedTaxon] = useState<string | null>(null);
  /** File paths in the order the open list view is showing them. */
  const [visibleOrder, setVisibleOrder] = useState<string[]>([]);
  const [exportModalMode, setExportModalMode] = useState<'batch' | 'concatenate' | 'group' | null>(null);
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
  const summariesRef = useRef<AlignmentSummary[]>([]);
  const catalogRequestIdRef = useRef(0);
  const activeRecipeRef = useRef(activeRecipe);
  activeRecipeRef.current = activeRecipe;
  summariesRef.current = summaries;

  const processingRecipeKey = sequenceProcessingRecipeKey(activeRecipe);
  // Hold the recipe identity steady while only assessment thresholds change, so
  // the dataset-wide recalculation does not restart for a gating-only edit.
  // A `useMemo` keyed on the string would have to return `activeRecipe` without
  // listing it as a dependency, which goes stale as soon as a new field is added
  // to the recipe; comparing the key against the stored recipe cannot.
  const processingRecipeRef = useRef(activeRecipe);
  if (sequenceProcessingRecipeKey(processingRecipeRef.current) !== processingRecipeKey) {
    processingRecipeRef.current = activeRecipe;
  }
  const processingRecipe = processingRecipeRef.current;
  const viewerProcessingRecipe = useMemo(
    () =>
      msaSidebarContext === 'orf'
        ? processingRecipe
        : recipeWithoutOrfAnalysis(processingRecipe),
    [msaSidebarContext, processingRecipe]
  );
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
    } else {
      // Browser environment: Check if ?run= is in URL for auto-loading
      const urlParams = new URLSearchParams(window.location.search);
      const runParam = urlParams.get('run');
      if (runParam) {
        import('./tauriClient').then(({ loadDirectoryFromUrl }) => {
          setIsLoading(true);
          setLoadingProgress(null);
          loadDirectoryFromUrl(runParam, processingRecipe, (cur, tot, file) => {
            setLoadingProgress({ current: cur, total: tot, percent: (cur / tot) * 100, fileName: file });
          }).then(({ dirName, scanResponse }) => {
            setCurrentPath(dirName);
            setSummaries(scanResponse.summaries);
            rawSummariesRef.current = scanResponse.summaries;
            setOverview(scanResponse.overview);
            setOccupancy(scanResponse.occupancy);
            setParseFailures(scanResponse.parse_failures ?? []);
            if (scanResponse.summaries.length > 0) {
              const first = scanResponse.summaries[0];
              setSelectedLocusId(first.id);
              setSelectedFilePath(first.file_path);
            }
            setReferenceRevision((r) => r + 1);
            setIsLoading(false);
            setLoadingProgress(null);
          }).catch((err: unknown) => {
            console.error('Failed to load example from URL:', err);
            setLoadError(
              `Could not load the example dataset: ${err instanceof Error ? err.message : String(err)}`
            );
            setIsLoading(false);
            setLoadingProgress(null);
          });
        });
      }
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
    setLoadError(null);
    setLoadingProgress({ current: 0, total: 0, percent: 0, fileName: 'Discovering alignment files…' });
    setCurrentPath(dirPath);
    try {
      const res = await scanDirectory(dirPath);
      rawSummariesRef.current = res.summaries;
      setSummaries(res.summaries);
      setOverview(res.overview);
      setOccupancy(res.occupancy);
      setParseFailures(res.parse_failures ?? []);

      if (res.summaries.length > 0) {
        const first = res.summaries[0];
        setSelectedLocusId(first.id);
        setSelectedFilePath(first.file_path);
      }
    } catch (e: unknown) {
      console.error('Failed to scan directory:', e);
      setLoadError(
        `Could not load directory: ${e instanceof Error ? e.message : String(e)}`
      );
      setCurrentPath(null);
    } finally {
      setIsLoading(false);
      setLoadingProgress(null);
    }
  };

  const handleFilesSelected = async (files: FileList | File[]) => {
    if (!files || files.length === 0) return;
    setIsLoading(true);
    setLoadError(null);
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
      setParseFailures(scanResponse.parse_failures ?? []);

      if (scanResponse.summaries.length > 0) {
        const first = scanResponse.summaries[0];
        setSelectedLocusId(first.id);
        setSelectedFilePath(first.file_path);
      }
    } catch (e: unknown) {
      console.error('Failed to load files:', e);
      setLoadError(
        `Could not load alignments: ${e instanceof Error ? e.message : String(e)}`
      );
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
  // Opening the viewer without a locus selected would otherwise leave the
  // spinner running with nothing loading.
  useEffect(() => {
    if (activeView !== 'msa' || selectedFilePath || summaries.length === 0) return;
    const first = summaries[0];
    setSelectedLocusId(first.id);
    setSelectedFilePath(first.file_path);
  }, [activeView, selectedFilePath, summaries]);

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
    getAlignment(selectedFilePath, viewerProcessingRecipe, totalUniqueTaxa)
      .then((data) => {
        if (!isCancelled) {
          setViewData(
            applyCatalogAssessmentToView(
              data,
              summariesRef.current.find((summary) => summary.file_path === selectedFilePath)
            )
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
  }, [activeView, selectedFilePath, viewerProcessingRecipe, overview.total_unique_taxa, occupancy.length]);

  // Keep an open preview synchronized to Catalog QC without applying ORF sample
  // loss or ORF candidate failure to the alignment assessment badge.
  useEffect(() => {
    setViewData((current) => {
      if (!current) return current;
      const summary = summaries.find(
        (candidate) => candidate.file_path === current.raw_alignment.file_path
      );
      return applyCatalogAssessmentToView(current, summary);
    });
  }, [summaries, gatingRecipeKey]);

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
    setCatalogError(null);

    const timer = window.setTimeout(() => {
      recalculateCatalog(paths, activeRecipeRef.current, totalUniqueTaxa, (progress) => {
        if (catalogRequestIdRef.current === requestId) {
          const nextPercent = Math.max(0, Math.min(100, progress.percent));
          setCatalogProcessingPercent((currentPercent) => Math.max(currentPercent, nextPercent));
        }
      })
        .then((res) => {
          if (catalogRequestIdRef.current !== requestId || res.summaries.length === 0) return;
          rawSummariesRef.current = res.summaries;
          setSummaries(res.summaries);
          setOverview(res.overview);
          setCatalogError(null);
        })
        .catch((err) => {
          if (catalogRequestIdRef.current === requestId) {
            console.error('Failed to recalculate catalog:', err);
            setCatalogError(err instanceof Error ? err.message : String(err));
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
  }, [
    processingRecipe,
    gatingRecipeKey,
    referenceRevision,
    currentPath,
    isLoading,
    overview.total_unique_taxa,
    occupancy.length,
  ]);

  /*
   * Keep this state reset close to loading so a dataset switch cannot leave a stale
   * processing badge behind while the new directory is being indexed.
   */
  useEffect(() => {
    if (isLoading) {
      setViewData(null);
    }
  }, [isLoading]);

  useEffect(() => {
    setShowParseFailures(parseFailures.length > 0);
  }, [parseFailures]);

  // The occupancy chart is the only view that needs per-sample retention, so it
  // is fetched when that view opens rather than carried on every recalculation.
  useEffect(() => {
    if (activeView !== 'qc' || isTauri || summaries.length === 0) {
      setRetentionSummaries(null);
      return;
    }
    let cancelled = false;
    import('./engine/enginePool')
      .then(({ enginePool }) => enginePool.retentionDetails())
      .then((details) => {
        if (cancelled || details.size === 0) return;
        setRetentionSummaries(
          summaries.map((summary) => {
            const detail = details.get(summary.file_path);
            return detail ? { ...summary, ...detail } : summary;
          })
        );
      })
      .catch((error) => console.warn('Could not load retention details:', error));
    return () => {
      cancelled = true;
    };
  }, [activeView, summaries]);

  /** Adds or removes a sample from the dataset-wide discard list. */
  const handleToggleDiscardTaxon = useCallback((taxon: string) => {
    setActiveRecipe((current) => {
      const discarded = current.discarded_taxa ?? [];
      const next = discarded.includes(taxon)
        ? discarded.filter((name) => name !== taxon)
        : [...discarded, taxon].sort();
      return { ...current, discarded_taxa: next };
    });
  }, []);

  const handleVisibleOrderChange = useCallback((filePaths: string[]) => {
    // Replace only on a real change, so a list re-render cannot loop.
    setVisibleOrder((current) =>
      current.length === filePaths.length && current.every((path, i) => path === filePaths[i])
        ? current
        : filePaths
    );
  }, []);

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
      setSelectedTaxon(null);
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

  // Stepping through alignments follows the order of the list the viewer was
  // opened from, so Next matches what the user just saw.
  const currentOrderIndex = selectedFilePath ? visibleOrder.indexOf(selectedFilePath) : -1;
  const stepAlignment = useCallback(
    (delta: number) => {
      if (currentOrderIndex < 0) return;
      const nextPath = visibleOrder[currentOrderIndex + delta];
      if (!nextPath) return;
      const summary = summariesRef.current.find((item) => item.file_path === nextPath);
      if (!summary) return;
      setViewData(null);
      setAlignmentLoadError(null);
      setSelectedTaxon(null);
      setSelectedLocusId(summary.id);
      setSelectedFilePath(summary.file_path);
    },
    [currentOrderIndex, visibleOrder]
  );

  // Drag & drop handlers.
  // `dragleave` fires when the pointer crosses into a child element, so count
  // enter and leave events instead of clearing the overlay on the first leave.
  const dragDepthRef = useRef(0);

  const handleDragEnter = (e: React.DragEvent) => {
    e.preventDefault();
    dragDepthRef.current += 1;
    setIsDragging(true);
  };

  const handleDragOver = (e: React.DragEvent) => {
    e.preventDefault();
  };

  const handleDragLeave = (e: React.DragEvent) => {
    e.preventDefault();
    dragDepthRef.current = Math.max(0, dragDepthRef.current - 1);
    if (dragDepthRef.current === 0) setIsDragging(false);
  };

  /**
   * Walks a dropped directory entry and collects every file inside it.
   * `DataTransfer.files` holds one entry with no content for a dropped folder,
   * so the folder has to be expanded through the entry API instead.
   */
  const collectFilesFromEntry = async (
    entry: FileSystemEntry,
    pathPrefix: string
  ): Promise<File[]> => {
    if (entry.isFile) {
      const fileEntry = entry as FileSystemFileEntry;
      const file = await new Promise<File | null>((resolve) =>
        fileEntry.file(resolve, () => resolve(null))
      );
      if (!file) return [];
      // Preserve the folder name so the dataset is labelled correctly.
      Object.defineProperty(file, 'webkitRelativePath', {
        value: `${pathPrefix}${file.name}`,
        configurable: true,
      });
      return [file];
    }

    if (entry.isDirectory) {
      const reader = (entry as FileSystemDirectoryEntry).createReader();
      const children: FileSystemEntry[] = [];
      // readEntries returns at most 100 entries per call.
      for (;;) {
        const batch = await new Promise<FileSystemEntry[]>((resolve) =>
          reader.readEntries(resolve, () => resolve([]))
        );
        if (batch.length === 0) break;
        children.push(...batch);
      }
      const nested = await Promise.all(
        children.map((child) => collectFilesFromEntry(child, `${pathPrefix}${entry.name}/`))
      );
      return nested.flat();
    }

    return [];
  };

  const handleDrop = async (e: React.DragEvent) => {
    e.preventDefault();
    dragDepthRef.current = 0;
    setIsDragging(false);

    const items = Array.from(e.dataTransfer.items ?? []);
    const entries = items
      .map((item) => (item.kind === 'file' ? item.webkitGetAsEntry?.() ?? null : null))
      .filter((entry): entry is FileSystemEntry => entry !== null);

    if (entries.length > 0) {
      try {
        const collected = await Promise.all(
          entries.map((entry) => collectFilesFromEntry(entry, ''))
        );
        const files = collected.flat();
        if (files.length > 0) {
          handleFilesSelected(files);
          return;
        }
      } catch (error) {
        console.error('Failed to read dropped folder:', error);
      }
    }

    // Plain file drops, and browsers without the entry API.
    if (e.dataTransfer.files && e.dataTransfer.files.length > 0) {
      handleFilesSelected(e.dataTransfer.files);
    }
  };

  return (
    <div
      className={`h-screen w-screen flex flex-col overflow-hidden ${isDarkTheme ? 'dark' : ''}`}
      onDragEnter={handleDragEnter}
      onDragOver={handleDragOver}
      onDragLeave={handleDragLeave}
      onDrop={handleDrop}
    >
      {/* Hidden File / Folder Input */}
      <input
        ref={fileInputRef}
        type="file"
        // Directory picking is not in the standard input typings.
        // @ts-expect-error - non-standard attributes for folder selection
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
            <p className="text-xs text-[#8b949e]">
              Accepts PHYLIP, FASTA, and NEXUS files
            </p>
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

      {/* Dataset could not be loaded */}
      {loadError && (
        <div className="shrink-0 bg-rose-500/10 border-b border-rose-500/30 px-4 py-2 text-xs text-rose-200 flex items-start gap-2">
          <XCircle className="w-4 h-4 shrink-0 mt-0.5 text-rose-400" />
          <div className="flex-1 min-w-0">
            <div className="font-semibold text-rose-300">Could not open that folder</div>
            <div className="text-rose-200/80 break-words">{loadError}</div>
          </div>
          <button
            onClick={() => setLoadError(null)}
            className="shrink-0 text-rose-300 hover:text-rose-100 transition-colors"
            aria-label="Dismiss"
          >
            <X className="w-4 h-4" />
          </button>
        </div>
      )}

      {catalogError && (
        <div className="shrink-0 bg-rose-500/10 border-b border-rose-500/30 px-4 py-2 text-xs text-rose-200 flex items-start gap-2">
          <XCircle className="w-4 h-4 shrink-0 mt-0.5 text-rose-400" />
          <div className="flex-1 min-w-0">
            <div className="font-semibold text-rose-300">Catalog calculation failed</div>
            <div className="text-rose-200/80 break-words">
              The displayed results use the last successful settings. {catalogError}
            </div>
          </div>
          <button
            onClick={() => setCatalogError(null)}
            className="shrink-0 text-rose-300 hover:text-rose-100 transition-colors"
            aria-label="Dismiss"
          >
            <X className="w-4 h-4" />
          </button>
        </div>
      )}

      {/* Files that could not be read */}
      {showParseFailures && parseFailures.length > 0 && (
        <div className="shrink-0 bg-amber-500/10 border-b border-amber-500/30 px-4 py-2 text-xs text-amber-200 flex items-start gap-2">
          <AlertTriangle className="w-4 h-4 shrink-0 mt-0.5 text-amber-400" />
          <div className="flex-1 min-w-0 space-y-1">
            <div className="font-semibold text-amber-300">
              {parseFailures.length === 1
                ? '1 file could not be read and was left out'
                : `${parseFailures.length} files could not be read and were left out`}
            </div>
            <ul className="space-y-0.5 max-h-24 overflow-y-auto font-mono text-[11px] text-amber-200/80">
              {parseFailures.slice(0, 20).map((failure) => (
                <li key={failure.file_path} className="truncate">
                  {failure.file_name}: {failure.error}
                </li>
              ))}
            </ul>
            {parseFailures.length > 20 && (
              <div className="text-[11px] text-amber-200/60">
                and {parseFailures.length - 20} more
              </div>
            )}
          </div>
          <button
            onClick={() => setShowParseFailures(false)}
            className="shrink-0 text-amber-300 hover:text-amber-100 transition-colors"
            aria-label="Dismiss"
          >
            <X className="w-4 h-4" />
          </button>
        </div>
      )}

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
                onVisibleOrderChange={handleVisibleOrderChange}
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
                onVisibleOrderChange={handleVisibleOrderChange}
              />
            )}

            {activeView === 'qc' && (
              <QcDistributions
                overview={overview}
                summaries={retentionSummaries ?? summaries}
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
                  discardedTaxa={activeRecipe.discarded_taxa}
                  onToggleDiscardTaxon={handleToggleDiscardTaxon}
                  selectedTaxon={selectedTaxon}
                  onSelectTaxon={setSelectedTaxon}
                  onPreviousAlignment={
                    currentOrderIndex > 0 ? () => stepAlignment(-1) : undefined
                  }
                  onNextAlignment={
                    currentOrderIndex >= 0 && currentOrderIndex < visibleOrder.length - 1
                      ? () => stepAlignment(1)
                      : undefined
                  }
                  alignmentPosition={
                    currentOrderIndex >= 0
                      ? { index: currentOrderIndex + 1, total: visibleOrder.length }
                      : undefined
                  }
                />
              ) : (
                <div className="flex-1 flex items-center justify-center bg-[#0e1014] text-[#8b949e]">
                  <div className="flex flex-col items-center gap-3 text-center px-6">
                    {summaries.length === 0 ? (
                      <>
                        <XCircle className="w-9 h-9 text-[#8b949e]" />
                        <div>
                          <p className="text-sm font-medium text-[#dce6ff]">No alignment selected</p>
                          <p className="text-xs text-[#8b949e] mt-1">
                            Open a folder, then choose a locus in the Catalog.
                          </p>
                        </div>
                      </>
                    ) : alignmentLoadError ? (
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
                    ? `${Math.round(loadingProgress.current).toLocaleString()} / ${loadingProgress.total.toLocaleString()} files`
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
