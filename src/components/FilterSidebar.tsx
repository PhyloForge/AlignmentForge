import React from 'react';
import {
  Star, Eraser, Crop, Filter, Activity, Dna, Columns, BarChart2,
  Scissors,
  CheckCircle2,
  AlertTriangle,
  RotateCcw,
  SlidersHorizontal,
  ChevronDown,
  ChevronRight,
  ShieldAlert,
  Download,
  Upload,
  Loader2,
} from 'lucide-react';
import {
  TrimmingRecipe,
  AmbiguityStrategy,
  RelativeWidth,
  OrfSearchMode,
  StatisticalColumnMethod,
  StatColGapTreatment,
  TrimalHeuristic,
} from '../types';

interface FilterSidebarProps {
  recipes: TrimmingRecipe[];
  recipe: TrimmingRecipe;
  onChangeRecipe: (updated: TrimmingRecipe) => void;
  onResetRecipe: () => void;
  onExportFilters: () => Promise<string | null>;
  onLoadFilters: () => Promise<string | null>;
  passedCount: number;
  totalCount: number;
}

export const FilterSidebar: React.FC<FilterSidebarProps> = ({
  recipes,
  recipe,
  onChangeRecipe,
  onResetRecipe,
  onExportFilters,
  onLoadFilters,
  passedCount,
  totalCount,
}) => {
  const [collapsedSections, setCollapsedSections] = React.useState<Record<string, boolean>>({});
  const [configAction, setConfigAction] = React.useState<'export' | 'load' | null>(null);
  const [configNotice, setConfigNotice] = React.useState<{
    kind: 'success' | 'error';
    message: string;
  } | null>(null);
  const [sampleOccupancyDraft, setSampleOccupancyDraft] = React.useState(
    recipe.min_sample_locus_occupancy_percent
  );
  const [orfSupportDraft, setOrfSupportDraft] = React.useState(
    recipe.orf_min_shared_support_percent ?? 90
  );
  const [orfSegmentDraft, setOrfSegmentDraft] = React.useState(
    recipe.orf_min_segment_aa ?? 35
  );
  const [orfCodingDraft, setOrfCodingDraft] = React.useState(
    recipe.orf_min_coding_score ?? 40
  );
  const lastSampleOccupancyThreshold = React.useRef(
    recipe.min_sample_locus_occupancy_percent > 0
      ? recipe.min_sample_locus_occupancy_percent
      : 10
  );

  React.useEffect(() => {
    setSampleOccupancyDraft(recipe.min_sample_locus_occupancy_percent);
    if (recipe.min_sample_locus_occupancy_percent > 0) {
      lastSampleOccupancyThreshold.current = recipe.min_sample_locus_occupancy_percent;
    }
  }, [recipe.min_sample_locus_occupancy_percent]);

  React.useEffect(() => {
    setOrfSupportDraft(recipe.orf_min_shared_support_percent ?? 90);
  }, [recipe.orf_min_shared_support_percent]);

  React.useEffect(() => {
    setOrfSegmentDraft(recipe.orf_min_segment_aa ?? 35);
  }, [recipe.orf_min_segment_aa]);

  React.useEffect(() => {
    setOrfCodingDraft(recipe.orf_min_coding_score ?? 40);
  }, [recipe.orf_min_coding_score]);

  const toggleSection = (key: string) => {
    setCollapsedSections((prev) => ({ ...prev, [key]: !prev[key] }));
  };

  const update = <K extends keyof TrimmingRecipe>(key: K, value: TrimmingRecipe[K]) => {
    onChangeRecipe({ ...recipe, [key]: value });
  };

  const commitSampleOccupancy = (value: number) => {
    if (value > 0) {
      lastSampleOccupancyThreshold.current = value;
    }
    if (value !== recipe.min_sample_locus_occupancy_percent) {
      update('min_sample_locus_occupancy_percent', value);
    }
  };

  const commitOrfSupport = (value: number) => {
    if (value !== recipe.orf_min_shared_support_percent) {
      update('orf_min_shared_support_percent', value);
    }
  };

  const commitOrfSegment = (value: number) => {
    if (value !== recipe.orf_min_segment_aa) {
      update('orf_min_segment_aa', value);
    }
  };

  const commitOrfCoding = (value: number) => {
    if (value !== recipe.orf_min_coding_score) {
      update('orf_min_coding_score', value);
    }
  };

  const toggleSampleOccupancy = (enabled: boolean) => {
    const nextValue = enabled ? lastSampleOccupancyThreshold.current : 0;
    if (!enabled && sampleOccupancyDraft > 0) {
      lastSampleOccupancyThreshold.current = sampleOccupancyDraft;
    }
    setSampleOccupancyDraft(nextValue);
    commitSampleOccupancy(nextValue);
  };

  const handleToggleOrf = (checked: boolean) => {
    if (checked) {
      // Automatically disable column-based filters to protect coding triplet frames
      onChangeRecipe({
        ...recipe,
        enable_orf: true,
        trim_columns: false,
        enable_statistical_columns: false,
        fail_if_no_orf: false,
      });
    } else {
      onChangeRecipe({
        ...recipe,
        enable_orf: false,
      });
    }
  };

  const passPercent = totalCount > 0 ? ((passedCount / totalCount) * 100).toFixed(1) : '100';

  const runConfigAction = async (action: 'export' | 'load') => {
    setConfigAction(action);
    setConfigNotice(null);
    try {
      const filePath = action === 'export' ? await onExportFilters() : await onLoadFilters();
      if (filePath) {
        const fileName = filePath.split(/[\\/]/).pop() || filePath;
        setConfigNotice({
          kind: 'success',
          message: action === 'export' ? `Exported ${fileName}` : `Loaded ${fileName}`,
        });
      }
    } catch (error) {
      setConfigNotice({
        kind: 'error',
        message: error instanceof Error ? error.message : String(error),
      });
    } finally {
      setConfigAction(null);
    }
  };

  return (
    <aside className="w-80 min-h-0 flex-none bg-[#14171d] border-r border-[#232833] flex flex-col h-full overflow-hidden select-none">
      {/* Sidebar Header with pass indicator */}
      <div className="flex-none p-3.5 border-b border-[#232833] bg-[#171b22]/70 flex items-center justify-between">
        <div className="flex items-center gap-2">
          <SlidersHorizontal className="w-4 h-4 text-blue-400" />
          <span className="font-semibold text-xs text-[#dce6ff] uppercase tracking-wider">
            Filtering Pipeline
          </span>
        </div>
        <button
          onClick={onResetRecipe}
          className="text-[#8b949e] hover:text-[#c9d1d9] text-[11px] flex items-center gap-1 transition-colors"
          title="Reset to recipe defaults"
        >
          <RotateCcw className="w-3 h-3" />
          <span>Reset</span>
        </button>
      </div>

      {/* Dataset Impact Badge */}
      <div className="flex-none px-3.5 py-2.5 bg-[#0e1014] border-b border-[#232833] flex items-center justify-between">
        <span className="text-[11px] text-[#8b949e]">Dataset Retention</span>
        <div className="flex items-center gap-2">
          <div className="w-20 bg-[#1f242e] h-2 rounded-full overflow-hidden">
            <div
              className="bg-emerald-500 h-full rounded-full transition-all duration-300"
              style={{ width: `${passPercent}%` }}
            />
          </div>
          <span className="font-mono text-xs font-semibold text-emerald-400">
            {passedCount}/{totalCount} ({passPercent}%)
          </span>
        </div>
      </div>

      {/* Filter configuration import / export */}
      <div className="flex-none px-3.5 py-2 bg-[#101319] border-b border-[#232833] space-y-1.5">
        <div className="grid grid-cols-2 gap-2">
          <button
            onClick={() => runConfigAction('load')}
            disabled={configAction !== null}
            className="flex items-center justify-center gap-1.5 px-2 py-1.5 rounded bg-[#192b23] hover:bg-[#203a2e] text-emerald-300 border border-emerald-500/30 text-[11px] font-medium transition-colors disabled:opacity-50 disabled:cursor-wait"
            title="Load filters from an editable AlignmentForge TOML file"
          >
            {configAction === 'load' ? (
              <Loader2 className="w-3 h-3 animate-spin" />
            ) : (
              <Upload className="w-3 h-3" />
            )}
            <span>Load Config</span>
          </button>
          <button
            onClick={() => runConfigAction('export')}
            disabled={configAction !== null}
            className="flex items-center justify-center gap-1.5 px-2 py-1.5 rounded bg-[#1d273b] hover:bg-[#25334d] text-blue-300 border border-blue-500/30 text-[11px] font-medium transition-colors disabled:opacity-50 disabled:cursor-wait"
            title="Save the current filters as an editable AlignmentForge TOML file"
          >
            {configAction === 'export' ? (
              <Loader2 className="w-3 h-3 animate-spin" />
            ) : (
              <Download className="w-3 h-3" />
            )}
            <span>Export Config</span>
          </button>
        </div>
        <div className="flex items-center gap-2 rounded border border-[#2d3545] bg-[#171b22] px-2 py-1.5">
          <span className="shrink-0 text-[10px] font-medium uppercase tracking-wide text-[#8b949e]">
            Preset
          </span>
          <select
            value={recipe.name}
            onChange={(event) => {
              const selected = recipes.find((candidate) => candidate.name === event.target.value);
              if (selected) onChangeRecipe(selected);
            }}
            className="min-w-0 flex-1 cursor-pointer bg-transparent text-[11px] font-medium text-[#c9d1d9] outline-none"
            aria-label="Filter preset"
          >
            {recipes.map((preset) => (
              <option
                key={preset.name}
                value={preset.name}
                className="bg-[#171b22] text-[#c9d1d9]"
              >
                {preset.name}
              </option>
            ))}
          </select>
        </div>
        
        {configNotice && (
          <div
            className={`flex items-start gap-1.5 rounded px-2 py-1.5 text-[10px] ${
              configNotice.kind === 'success'
                ? 'bg-emerald-500/10 text-emerald-300'
                : 'bg-red-500/10 text-red-300'
            }`}
          >
            {configNotice.kind === 'success' ? (
              <CheckCircle2 className="w-3 h-3 shrink-0 mt-0.5" />
            ) : (
              <ShieldAlert className="w-3 h-3 shrink-0 mt-0.5" />
            )}
            <span className="break-words min-w-0">{configNotice.message}</span>
          </div>
        )}
      </div>

      {/* Scrollable Filters List */}
      <div className="min-h-0 flex-1 overflow-y-auto overscroll-contain px-3.5 py-3 flex flex-col gap-4 text-xs">
        {/* Locus Pass / Fail Criteria */}
        <div
          className="shrink-0 border border-emerald-500/30 rounded-lg bg-[#171b22]/60 overflow-hidden shadow-sm"
          style={{ order: 0 }}
        >
          <button
            onClick={() => toggleSection('gating')}
            className="w-full px-3 py-2 bg-[#1b2029]/80 flex items-center justify-between text-[#dce6ff] font-medium"
          >
            <div className="flex items-center gap-2">
              <span className="w-5 h-5 rounded flex items-center justify-center text-emerald-400"><Star className="w-3.5 h-3.5" /></span>
              <span className="font-semibold text-emerald-300">Locus Pass / Fail Criteria</span>
            </div>
            {collapsedSections['gating'] ? (
              <ChevronRight className="w-3.5 h-3.5 text-[#8b949e]" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5 text-[#8b949e]" />
            )}
          </button>

          {!collapsedSections['gating'] && (
            <div className="p-3 space-y-3">
              <label className="flex items-center justify-between cursor-pointer">
                <span className="text-[#c9d1d9] font-medium">Enable Quality Gating</span>
                <input
                  type="checkbox"
                  checked={recipe.assess_alignment}
                  onChange={(e) => update('assess_alignment', e.target.checked)}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-emerald-500 focus:ring-0 cursor-pointer"
                />
              </label>

              {recipe.assess_alignment && (
                <>
                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Min Taxon Occupancy %</span>
                      <span className="font-mono text-emerald-400 font-semibold">
                        {recipe.min_taxa_occupancy_percent ?? 50}%
                      </span>
                    </div>
                    <input
                      type="range"
                      min="0"
                      max="100"
                      step="5"
                      value={recipe.min_taxa_occupancy_percent ?? 50}
                      onChange={(e) => update('min_taxa_occupancy_percent', parseFloat(e.target.value))}
                      className="w-full accent-emerald-500"
                    />
                    <div className="text-[10px] text-[#8b949e] mt-0.5">
                      {(recipe.min_taxa_occupancy_percent ?? 50) === 0
                        ? 'Occupancy filter disabled'
                        : `Fails loci missing >${100 - (recipe.min_taxa_occupancy_percent ?? 50)}% of dataset taxa`}
                    </div>
                  </div>

                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Min Surviving Taxa (Count)</span>
                      <span className="font-mono text-emerald-400 font-semibold">
                        {recipe.min_taxa} taxa
                      </span>
                    </div>
                    <input
                      type="range"
                      min="2"
                      max="50"
                      step="1"
                      value={recipe.min_taxa}
                      onChange={(e) => update('min_taxa', parseInt(e.target.value))}
                      className="w-full accent-emerald-500"
                    />
                  </div>

                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Min Post-Trim Length</span>
                      <span className="font-mono text-emerald-400 font-semibold">
                        {recipe.min_length} bp
                      </span>
                    </div>
                    <input
                      type="range"
                      min="50"
                      max="1000"
                      step="25"
                      value={recipe.min_length}
                      onChange={(e) => update('min_length', parseInt(e.target.value))}
                      className="w-full accent-emerald-500"
                    />
                  </div>

                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Max Overall Gap %</span>
                      <span className="font-mono text-emerald-400 font-semibold">
                        {recipe.max_gap_percent}%
                      </span>
                    </div>
                    <input
                      type="range"
                      min="10"
                      max="90"
                      step="5"
                      value={recipe.max_gap_percent}
                      onChange={(e) => update('max_gap_percent', parseFloat(e.target.value))}
                      className="w-full accent-emerald-500"
                    />
                  </div>

                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Min Variable Sites</span>
                      <span className="font-mono text-cyan-400 font-semibold">
                        {recipe.min_variable_count ?? 0}
                      </span>
                    </div>
                    <input
                      type="range"
                      min="0"
                      max="1000"
                      step="1"
                      value={recipe.min_variable_count ?? 0}
                      onChange={(e) => update('min_variable_count', parseInt(e.target.value))}
                      className="w-full accent-cyan-500"
                    />
                    <div className="text-[10px] text-[#8b949e] mt-0.5">
                      {(recipe.min_variable_count ?? 0) === 0
                        ? 'Variable-site count filter disabled'
                        : 'Minimum sites with at least two resolved nucleotide states'}
                    </div>
                  </div>

                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Min Variable-Site Proportion</span>
                      <span className="font-mono text-cyan-400 font-semibold">
                        {(recipe.min_variable_percent ?? 0).toFixed(1)}%
                      </span>
                    </div>
                    <input
                      type="range"
                      min="0"
                      max="100"
                      step="0.5"
                      value={recipe.min_variable_percent ?? 0}
                      onChange={(e) => update('min_variable_percent', parseFloat(e.target.value))}
                      className="w-full accent-cyan-500"
                    />
                    <div className="text-[10px] text-[#8b949e] mt-0.5">
                      {(recipe.min_variable_percent ?? 0) === 0
                        ? 'Variable-site proportion filter disabled'
                        : 'Gaps, N/n, ?, and ambiguity codes are treated as missing'}
                    </div>
                  </div>

                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Min Parsimony-Informative Sites</span>
                      <span className="font-mono text-emerald-400 font-semibold">
                        {recipe.min_pis_count ?? 0}
                      </span>
                    </div>
                    <input
                      type="range"
                      min="0"
                      max="1000"
                      step="1"
                      value={recipe.min_pis_count ?? 0}
                      onChange={(e) => update('min_pis_count', parseInt(e.target.value))}
                      className="w-full accent-emerald-500"
                    />
                    <div className="text-[10px] text-[#8b949e] mt-0.5">
                      {(recipe.min_pis_count ?? 0) === 0 ? 'PIS count filter disabled' : 'Minimum post-trim PIS count'}
                    </div>
                  </div>

                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Min Parsimony-Informative Proportion</span>
                      <span className="font-mono text-emerald-400 font-semibold">
                        {(recipe.min_pis_percent ?? 0).toFixed(1)}%
                      </span>
                    </div>
                    <input
                      type="range"
                      min="0"
                      max="100"
                      step="0.5"
                      value={recipe.min_pis_percent ?? 0}
                      onChange={(e) => update('min_pis_percent', parseFloat(e.target.value))}
                      className="w-full accent-emerald-500"
                    />
                    <div className="text-[10px] text-[#8b949e] mt-0.5">
                      {(recipe.min_pis_percent ?? 0) === 0 ? 'PIS proportion filter disabled' : 'Minimum post-trim PIS percentage'}
                    </div>
                  </div>
                </>
              )}
            </div>
          )}
        </div>

        {/* Display 2: Sanitation */}
        <div
          className="shrink-0 border border-[#232833] rounded-lg bg-[#171b22]/40 overflow-hidden"
          style={{ order: 2 }}
        >
          <button
            onClick={() => toggleSection('sanitation')}
            className="w-full px-3 py-2 bg-[#1b2029]/60 flex items-center justify-between text-[#c9d1d9] font-medium"
          >
            <div className="flex items-center gap-2">
              <span className="w-5 h-5 rounded flex items-center justify-center text-blue-400"><Eraser className="w-3.5 h-3.5" /></span>
              <span>Sanitation & Ambiguities</span>
            </div>
            {collapsedSections['sanitation'] ? (
              <ChevronRight className="w-3.5 h-3.5 text-[#8b949e]" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5 text-[#8b949e]" />
            )}
          </button>

          {!collapsedSections['sanitation'] && (
            <div className="p-3 space-y-3">
              <label className="flex items-center gap-2 cursor-pointer">
                <input
                  type="checkbox"
                  checked={recipe.replace_n_with_gap}
                  onChange={(e) => update('replace_n_with_gap', e.target.checked)}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-blue-500 focus:ring-0"
                />
                <span className="text-[#c9d1d9]">Replace 'N' / '?' with '-'</span>
              </label>

              <label className="flex items-center gap-2 cursor-pointer">
                <input
                  type="checkbox"
                  checked={recipe.remove_gap_only_columns ?? true}
                  onChange={(e) => update('remove_gap_only_columns', e.target.checked)}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-blue-500 focus:ring-0"
                />
                <span className="text-[#c9d1d9]">Remove Gap-Only Columns</span>
              </label>
              <p className="-mt-2 pl-6 text-[10px] leading-snug text-[#8b949e]">
                Removes columns containing only gaps, N, or ? across every sample.
              </p>

              <div>
                <label className="text-[11px] text-[#8b949e] block mb-1">IUPAC Ambiguity Mode</label>
                <select
                  value={recipe.ambiguity_strategy}
                  onChange={(e) =>
                    update('ambiguity_strategy', e.target.value as AmbiguityStrategy)
                  }
                  className="w-full bg-[#1b2029] border border-[#2d3545] rounded px-2 py-1 text-xs text-[#c9d1d9] outline-none"
                >
                  <option value="keep">Keep as IUPAC codes</option>
                  <option value="fixedstandard">Standard A/T Priority</option>
                  <option value="majoritybase">Majority Column Base</option>
                  <option value="converttogap">Convert to Gaps</option>
                </select>
              </div>
            </div>
          )}
        </div>

        {/* Display 5: Profile HMM Segment Cleaner (TAPIR-Style) */}
        <div
          className="shrink-0 border border-[#232833] rounded-lg bg-[#171b22]/40 overflow-hidden"
          style={{ order: 5 }}
        >
          <button
            onClick={() => toggleSection('hmm')}
            className="w-full px-3 py-2 bg-[#1b2029]/60 flex items-center justify-between text-[#c9d1d9] font-medium"
          >
            <div className="flex items-center gap-2">
              <span className="w-5 h-5 rounded flex items-center justify-center text-indigo-400"><Activity className="w-3.5 h-3.5" /></span>
              <span>Profile HMM Cleaner</span>
            </div>
            {collapsedSections['hmm'] ? (
              <ChevronRight className="w-3.5 h-3.5 text-[#8b949e]" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5 text-[#8b949e]" />
            )}
          </button>

          {!collapsedSections['hmm'] && (
            <div className="p-3 space-y-3">
              <label className="flex items-center justify-between cursor-pointer">
                <span className="text-[#c9d1d9]">Profile HMM Posterior Masking</span>
                <input
                  type="checkbox"
                  checked={recipe.trim_hmm ?? false}
                  onChange={(e) => update('trim_hmm', e.target.checked)}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-indigo-500 focus:ring-0"
                />
              </label>

              {(recipe.trim_hmm ?? false) && (
                <>
                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Min Posterior Match Confidence</span>
                      <span className="font-mono text-indigo-400 font-semibold">
                        {((recipe.hmm_min_posterior ?? 0.45) * 100).toFixed(0)}%
                      </span>
                    </div>
                    <input
                      type="range"
                      min="0.10"
                      max="0.90"
                      step="0.05"
                      value={recipe.hmm_min_posterior ?? 0.45}
                      onChange={(e) => update('hmm_min_posterior', parseFloat(e.target.value))}
                      className="w-full accent-indigo-500"
                    />
                    <div className="text-[10px] text-[#8b949e] mt-0.5">
                      Residues with posterior match confidence below threshold are flagged as misaligned.
                    </div>
                  </div>

                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Min Aberrant Segment Length</span>
                      <span className="font-mono text-indigo-400 font-semibold">
                        {recipe.hmm_min_segment_length ?? 8} bp
                      </span>
                    </div>
                    <input
                      type="range"
                      min="4"
                      max="30"
                      step="1"
                      value={recipe.hmm_min_segment_length ?? 8}
                      onChange={(e) => update('hmm_min_segment_length', parseInt(e.target.value))}
                      className="w-full accent-indigo-500"
                    />
                    <div className="text-[10px] text-[#8b949e] mt-0.5">
                      Requires at least this many contiguous low-confidence residues to trigger masking.
                    </div>
                  </div>

                  <div className="pt-2">
                    <div className="flex items-center justify-between mb-1 text-[11px]">
                      <span className="text-[#8b949e]">Min Distance Between Masked</span>
                      <span className="font-mono text-indigo-400 font-semibold">
                        {recipe.hmm_min_island_length ?? 20} bp
                      </span>
                    </div>
                    <input
                      type="range"
                      min="0"
                      max="100"
                      step="1"
                      value={recipe.hmm_min_island_length ?? 20}
                      onChange={(e) => update('hmm_min_island_length', parseInt(e.target.value))}
                      className="w-full accent-indigo-500"
                    />
                    <div className="text-[10px] text-[#8b949e] mt-0.5">
                      Merges masked segments separated by fewer than this many good bases.
                    </div>
                  </div>
                </>
              )}
            </div>
          )}
        </div>

        {/* Display 6: Sliding Window Segment Masking */}
        <div
          className="shrink-0 border border-[#232833] rounded-lg bg-[#171b22]/40 overflow-hidden"
          style={{ order: 6 }}
        >
          <button
            onClick={() => toggleSection('segments')}
            className="w-full px-3 py-2 bg-[#1b2029]/60 flex items-center justify-between text-[#c9d1d9] font-medium"
          >
            <div className="flex items-center gap-2">
              <span className="w-5 h-5 rounded flex items-center justify-center text-amber-400"><Scissors className="w-3.5 h-3.5" /></span>
              <span>Sliding Window Segment Mask</span>
            </div>
            {collapsedSections['segments'] ? (
              <ChevronRight className="w-3.5 h-3.5 text-[#8b949e]" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5 text-[#8b949e]" />
            )}
          </button>

          {!collapsedSections['segments'] && (
            <div className="p-3 space-y-3">
              <label className="flex items-center justify-between cursor-pointer">
                <span className="text-[#c9d1d9]">Mask Divergent Windows</span>
                <input
                  type="checkbox"
                  checked={recipe.trim_segments}
                  onChange={(e) => update('trim_segments', e.target.checked)}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-amber-500 focus:ring-0"
                />
              </label>

              {recipe.trim_segments && (
                <>
                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Window Size (bp)</span>
                      <span className="font-mono text-amber-400 font-semibold">
                        {recipe.segment_window_size} bp
                      </span>
                    </div>
                    <input
                      type="range"
                      min="30"
                      max="300"
                      step="10"
                      value={recipe.segment_window_size}
                      onChange={(e) => update('segment_window_size', parseInt(e.target.value))}
                      className="w-full accent-amber-500"
                    />
                  </div>

                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Window Divergence Cutoff</span>
                      <span className="font-mono text-amber-400 font-semibold">
                        {(recipe.segment_threshold * 100).toFixed(0)}%
                      </span>
                    </div>
                    <input
                      type="range"
                      min="0.20"
                      max="0.80"
                      step="0.01"
                      value={recipe.segment_threshold}
                      onChange={(e) => update('segment_threshold', parseFloat(e.target.value))}
                      className="w-full accent-amber-500"
                    />
                  </div>
                </>
              )}
            </div>
          )}
        </div>

        {/* Display 7: Candidate Open Reading Frame Extraction & Codon QC (Exons) */}
        <div
          className="hidden shrink-0 border border-[#232833] rounded-lg bg-[#171b22]/40 overflow-hidden"
          style={{ order: 7 }}
        >
          <button
            onClick={() => toggleSection('orf')}
            className="w-full px-3 py-2 bg-[#1b2029]/60 flex items-center justify-between text-[#c9d1d9] font-medium"
          >
            <div className="flex items-center gap-2">
              <span className="w-5 h-5 rounded flex items-center justify-center text-emerald-400"><Dna className="w-3.5 h-3.5" /></span>
              <span>Candidate ORF Extraction & Codons</span>
            </div>
            {collapsedSections['orf'] ? (
              <ChevronRight className="w-3.5 h-3.5 text-[#8b949e]" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5 text-[#8b949e]" />
            )}
          </button>

          {!collapsedSections['orf'] && (
            <div className="p-3 space-y-3">
              <label className="flex items-center justify-between cursor-pointer">
                <span className="text-[#c9d1d9] font-medium">Extract Candidate Coding Region</span>
                <input
                  type="checkbox"
                  checked={recipe.enable_orf ?? false}
                  onChange={(e) => handleToggleOrf(e.target.checked)}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-emerald-500 focus:ring-0 cursor-pointer"
                />
              </label>

              {(recipe.enable_orf ?? false) && (
                <>
                  <div className="bg-amber-500/10 border border-amber-500/20 rounded-md p-2.5 flex items-start gap-2 text-[11px] text-amber-300">
                    <AlertTriangle className="w-4 h-4 text-amber-400 shrink-0 mt-0.5" />
                    <div>
                      <span className="font-semibold block text-amber-200 mb-0.5">Candidate Extraction Active</span>
                      Continuous CDS mode is for exon-only alignments. Candidate ORF Extraction finds one shared region in an exon with attached intronic or flanking sequence. It combines stop-codon support with internal protein-profile evidence, but is not independent proof of gene annotation and does not splice multiple exons.
                    </div>
                  </div>

                  <div>
                    <span className="text-[11px] text-[#8b949e] block mb-1">ORF Search Mode</span>
                    <select
                      value={recipe.orf_search_mode ?? 'continuouscds'}
                      onChange={(e) =>
                        update('orf_search_mode', e.target.value as OrfSearchMode)
                      }
                      className="w-full bg-[#13161c] border border-[#2d3545] rounded px-2 py-1 text-xs text-[#c9d1d9] focus:outline-none focus:border-emerald-500"
                    >
                      <option value="continuouscds">Continuous CDS (exon-only)</option>
                      <option value="bestsharedsegment">Candidate ORF Extraction (exon + flanks)</option>
                    </select>
                  </div>

                  {(recipe.orf_search_mode ?? 'continuouscds') === 'continuouscds' ? (
                    <label className="flex items-center justify-between cursor-pointer text-[11px]">
                      <span className="text-[#c9d1d9]">Trim to Open Reading Frame (Lock Triplets)</span>
                      <input
                        type="checkbox"
                        checked={recipe.auto_shift_frame ?? true}
                        onChange={(e) => update('auto_shift_frame', e.target.checked)}
                        className="rounded bg-[#1f242e] border-[#2d3545] text-emerald-500 focus:ring-0"
                      />
                    </label>
                  ) : (
                    <>
                      <div>
                        <div className="flex justify-between text-[11px] mb-1">
                          <span className="text-[#8b949e]">Minimum Sample Support</span>
                          <span className="font-mono text-emerald-400 font-semibold">
                            {orfSupportDraft}%
                          </span>
                        </div>
                        <input
                          type="range"
                          min="10"
                          max="100"
                          step="5"
                          value={orfSupportDraft}
                          onChange={(e) => setOrfSupportDraft(parseFloat(e.target.value))}
                          onMouseUp={(e) => commitOrfSupport(parseFloat(e.currentTarget.value))}
                          onTouchEnd={(e) => commitOrfSupport(parseFloat(e.currentTarget.value))}
                          onKeyUp={(e) => commitOrfSupport(parseFloat(e.currentTarget.value))}
                          className="w-full accent-emerald-500"
                        />
                      </div>

                      <div>
                        <div className="flex justify-between text-[11px] mb-1">
                          <span className="text-[#8b949e]">Minimum Segment Length</span>
                          <span className="font-mono text-emerald-400 font-semibold">
                            {orfSegmentDraft} aa
                          </span>
                        </div>
                        <input
                          type="range"
                          min="5"
                          max="300"
                          step="5"
                          value={orfSegmentDraft}
                          onChange={(e) => setOrfSegmentDraft(parseInt(e.target.value))}
                          onMouseUp={(e) => commitOrfSegment(parseInt(e.currentTarget.value))}
                          onTouchEnd={(e) => commitOrfSegment(parseInt(e.currentTarget.value))}
                          onKeyUp={(e) => commitOrfSegment(parseInt(e.currentTarget.value))}
                          className="w-full accent-emerald-500"
                        />
                      </div>

                      <div>
                        <div className="flex justify-between text-[11px] mb-1">
                          <span className="text-[#8b949e]">Minimum Coding Evidence</span>
                          <span className="font-mono text-emerald-400 font-semibold">
                            {orfCodingDraft}/100
                          </span>
                        </div>
                        <input
                          type="range"
                          min="0"
                          max="100"
                          step="5"
                          value={orfCodingDraft}
                          onChange={(e) => setOrfCodingDraft(parseFloat(e.target.value))}
                          onMouseUp={(e) => commitOrfCoding(parseFloat(e.currentTarget.value))}
                          onTouchEnd={(e) => commitOrfCoding(parseFloat(e.currentTarget.value))}
                          onKeyUp={(e) => commitOrfCoding(parseFloat(e.currentTarget.value))}
                          className="w-full accent-emerald-500"
                        />
                        <p className="mt-1 text-[9px] leading-snug text-[#6e7681]">
                          Uses amino-acid conservation, separation from competing frames, synonymous change, and codon-position periodicity.
                        </p>
                      </div>
                    </>
                  )}


                  <label className="flex items-center justify-between cursor-pointer text-[11px]">
                    <span className="text-[#c9d1d9]">Auto-Flip Reverse Strand Exons</span>
                    <input
                      type="checkbox"
                      checked={recipe.auto_flip_reverse ?? true}
                      onChange={(e) => update('auto_flip_reverse', e.target.checked)}
                      className="rounded bg-[#1f242e] border-[#2d3545] text-emerald-500 focus:ring-0"
                    />
                  </label>

                  <label className="flex items-center justify-between cursor-pointer text-[11px]">
                    <span className="text-[#c9d1d9]">Skip UCE / Explicitly Non-Coding Loci</span>
                    <input
                      type="checkbox"
                      checked={recipe.exclude_uce ?? true}
                      onChange={(e) => update('exclude_uce', e.target.checked)}
                      className="rounded bg-[#1f242e] border-[#2d3545] text-emerald-500 focus:ring-0"
                    />
                  </label>
                  <p className="text-[10px] leading-snug text-[#8b949e]">
                    Continuous CDS also skips names marked intron, supercontig, or flanking.
                    Candidate ORF Extraction processes those names but still skips UCE, non-coding,
                    and intergenic alignments.
                  </p>

                  <div>
                    <span className="text-[11px] text-[#8b949e] block mb-1">Sample Stop Codon Action</span>
                    <select
                      value={recipe.stop_codon_action ?? 'removesample'}
                      onChange={(e) => update('stop_codon_action', e.target.value as any)}
                      className="w-full bg-[#13161c] border border-[#2d3545] rounded px-2 py-1 text-xs text-[#c9d1d9] focus:outline-none focus:border-emerald-500"
                    >
                      <option value="removesample">Prune Sample (Remove if internal stop codon)</option>
                      <option value="maskcodon">Mask Codon (Replace stops with '---')</option>
                      <option value="keep">Keep (Annotate only)</option>
                    </select>
                  </div>

                  <div>
                    <span className="text-[11px] text-[#8b949e] block mb-1">Genetic Code Table</span>
                    <select
                      value={recipe.genetic_code ?? 'standard'}
                      onChange={(e) => update('genetic_code', e.target.value as any)}
                      className="w-full bg-[#13161c] border border-[#2d3545] rounded px-2 py-1 text-xs text-[#c9d1d9] focus:outline-none focus:border-emerald-500"
                    >
                      <option value="standard">Standard Nuclear (Table 1)</option>
                      <option value="vertebratemitochondrial">Vertebrate Mitochondrial (Table 2)</option>
                      <option value="invertebratemitochondrial">Invertebrate Mitochondrial (Table 5)</option>
                    </select>
                  </div>
                </>
              )}
            </div>
          )}
        </div>

        {/* Display 3: External Edge Trimming */}
        <div
          className="shrink-0 border border-[#232833] rounded-lg bg-[#171b22]/40 overflow-hidden"
          style={{ order: 3 }}
        >
          <button
            onClick={() => toggleSection('external')}
            className="w-full px-3 py-2 bg-[#1b2029]/60 flex items-center justify-between text-[#c9d1d9] font-medium"
          >
            <div className="flex items-center gap-2">
              <span className="w-5 h-5 rounded flex items-center justify-center text-cyan-400"><Crop className="w-3.5 h-3.5" /></span>
              <span>Ragged Edge Trimming</span>
            </div>
            {collapsedSections['external'] ? (
              <ChevronRight className="w-3.5 h-3.5 text-[#8b949e]" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5 text-[#8b949e]" />
            )}
          </button>

          {!collapsedSections['external'] && (
            <div className="p-3 space-y-3">
              <label className="flex items-center justify-between cursor-pointer">
                <span className="text-[#c9d1d9]">Trim 5' and 3' Ragged Ends</span>
                <input
                  type="checkbox"
                  checked={recipe.trim_external}
                  onChange={(e) => update('trim_external', e.target.checked)}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-cyan-500 focus:ring-0"
                />
              </label>

              {recipe.trim_external && (
                <>
                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Minimum Taxa Occupancy</span>
                      <span className="font-mono text-cyan-400 font-semibold">
                        {recipe.min_external_percent}%
                      </span>
                    </div>
                    <input
                      type="range"
                      min="10"
                      max="100"
                      step="5"
                      value={recipe.min_external_percent}
                      onChange={(e) => update('min_external_percent', parseFloat(e.target.value))}
                      className="w-full accent-cyan-500"
                    />
                  </div>

                  <label className="flex items-center gap-2 cursor-pointer pt-1">
                    <input
                      type="checkbox"
                      checked={recipe.enable_orf ? true : recipe.codon_preserving}
                      disabled={recipe.enable_orf}
                      onChange={(e) => update('codon_preserving', e.target.checked)}
                      className="rounded bg-[#1f242e] border-[#2d3545] text-cyan-500 focus:ring-0 disabled:opacity-50 disabled:cursor-not-allowed"
                    />
                    <span className="text-[#c9d1d9]">Codon-Preserving Frame Snapping</span>
                  </label>
                </>
              )}
            </div>
          )}
        </div>

        {/* Display 8: Column Gap Filter */}
        <div
          className="shrink-0 border border-[#232833] rounded-lg bg-[#171b22]/40 overflow-hidden"
          style={{ order: 8 }}
        >
          <button
            onClick={() => toggleSection('columns')}
            className="w-full px-3 py-2 bg-[#1b2029]/60 flex items-center justify-between text-[#c9d1d9] font-medium"
          >
            <div className="flex items-center gap-2">
              <span className="w-5 h-5 rounded flex items-center justify-center text-rose-400"><Columns className="w-3.5 h-3.5" /></span>
              <span>Column Gap Filter</span>
            </div>
            {collapsedSections['columns'] ? (
              <ChevronRight className="w-3.5 h-3.5 text-[#8b949e]" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5 text-[#8b949e]" />
            )}
          </button>

          {!collapsedSections['columns'] && (
            <div className="p-3 space-y-3">
              {recipe.enable_orf && (
                <div className="bg-amber-500/10 border border-amber-500/20 rounded-md p-2 flex items-start gap-2 text-[10px] text-amber-300">
                  <AlertTriangle className="w-3.5 h-3.5 text-amber-400 shrink-0 mt-0.5" />
                  <span>
                    Candidate ORF extraction is active: Column gap trimming is automatically bypassed on loci processed as coding. {recipe.exclude_uce ? ((recipe.orf_search_mode ?? 'continuouscds') === 'bestsharedsegment' ? 'It still runs on skipped UCE, non-coding, and intergenic loci.' : 'It still runs on skipped UCE, intron, supercontig, flanking, non-coding, and intergenic loci.') : ''}
                  </span>
                </div>
              )}

              <label className="flex items-center justify-between cursor-pointer">
                <span className="text-[#c9d1d9]">Remove High-Gap Columns</span>
                <input
                  type="checkbox"
                  checked={recipe.trim_columns}
                  onChange={(e) => update('trim_columns', e.target.checked)}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-rose-500 focus:ring-0"
                />
              </label>

              {recipe.trim_columns && (
                <>
                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Max Allowed Column Gap %</span>
                      <span className="font-mono text-rose-400 font-semibold">
                        {recipe.min_column_gap_percent}%
                      </span>
                    </div>
                    <input
                      type="range"
                      min="10"
                      max="100"
                      step="5"
                      value={recipe.min_column_gap_percent}
                      onChange={(e) =>
                        update('min_column_gap_percent', parseFloat(e.target.value))
                      }
                      className="w-full accent-rose-500"
                    />
                  </div>

                  <label className="flex items-center gap-2 cursor-pointer pt-1">
                    <input
                      type="checkbox"
                      checked={recipe.count_n_as_gap}
                      onChange={(e) => update('count_n_as_gap', e.target.checked)}
                      className="rounded bg-[#1f242e] border-[#2d3545] text-rose-500 focus:ring-0"
                    />
                    <span className="text-[#c9d1d9]">Count 'N' and '?' as gaps</span>
                  </label>
                </>
              )}
            </div>
          )}
        </div>

        {/* Display 9: Statistical Column Trimming (trimAl & Gblocks) */}
        <div
          className="shrink-0 border border-[#232833] rounded-lg bg-[#171b22]/40 overflow-hidden"
          style={{ order: 9 }}
        >
          <button
            onClick={() => toggleSection('stat_columns')}
            className="w-full px-3 py-2 bg-[#1b2029]/60 flex items-center justify-between text-[#c9d1d9] font-medium"
          >
            <div className="flex items-center gap-2">
              <span className="w-5 h-5 rounded flex items-center justify-center text-indigo-400"><BarChart2 className="w-3.5 h-3.5" /></span>
              <span>Statistical Column Trimming</span>
            </div>
            {collapsedSections['stat_columns'] ? (
              <ChevronRight className="w-3.5 h-3.5 text-[#8b949e]" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5 text-[#8b949e]" />
            )}
          </button>

          {!collapsedSections['stat_columns'] && (
            <div className="p-3 space-y-3">
              {recipe.enable_orf && (
                <div className="bg-amber-500/10 border border-amber-500/20 rounded-md p-2 flex items-start gap-2 text-[10px] text-amber-300">
                  <AlertTriangle className="w-3.5 h-3.5 text-amber-400 shrink-0 mt-0.5" />
                  <span>
                    Candidate ORF extraction is active: Statistical column trimming is automatically bypassed on loci processed as coding. {recipe.exclude_uce ? ((recipe.orf_search_mode ?? 'continuouscds') === 'bestsharedsegment' ? 'It still runs on skipped UCE, non-coding, and intergenic loci.' : 'It still runs on skipped UCE, intron, supercontig, flanking, non-coding, and intergenic loci.') : ''}
                  </span>
                </div>
              )}

              <label className="flex items-center justify-between cursor-pointer">
                <span className="text-[#c9d1d9]">Enable Statistical Trimming</span>
                <input
                  type="checkbox"
                  checked={recipe.enable_statistical_columns}
                  onChange={(e) => update('enable_statistical_columns', e.target.checked)}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-indigo-500 focus:ring-0 cursor-pointer"
                />
              </label>

              {recipe.enable_statistical_columns && (
                <>
                  <div>
                    <label className="text-[11px] text-[#8b949e] block mb-1">Trimming Algorithm</label>
                    <select
                      value={recipe.stat_col_method || 'trimalsimilarity'}
                      onChange={(e) =>
                        update('stat_col_method', e.target.value as StatisticalColumnMethod)
                      }
                      className="w-full bg-[#1b2029] border border-[#2d3545] rounded px-2 py-1 text-xs text-[#c9d1d9] outline-none"
                    >
                      <option value="trimalsimilarity">trimAl Similarity & Consistency</option>
                      <option value="gblocksblocks">Gblocks Conserved Blocks</option>
                      <option value="entropy">Shannon Information Entropy</option>
                    </select>
                  </div>

                  {recipe.stat_col_method === 'trimalsimilarity' && (
                    <>
                      <div>
                        <label className="text-[11px] text-[#8b949e] block mb-1">trimAl Mode / Heuristic</label>
                        <select
                          value={recipe.stat_col_heuristic || 'custom'}
                          onChange={(e) =>
                            update('stat_col_heuristic', e.target.value as TrimalHeuristic)
                          }
                          className="w-full bg-[#1b2029] border border-[#2d3545] rounded px-2 py-1 text-xs text-[#c9d1d9] outline-none"
                        >
                          <option value="custom">Custom Similarity Cutoff</option>
                          <option value="gappyout">trimAl Gappyout (Gap Elbow)</option>
                          <option value="strict">trimAl Strict</option>
                          <option value="strictplus">trimAl StrictPlus (High Conservation)</option>
                        </select>
                      </div>

                      {recipe.stat_col_heuristic === 'custom' && (
                        <div>
                          <div className="flex justify-between text-[11px] mb-1">
                            <span className="text-[#8b949e]">Min Column Similarity</span>
                            <span className="font-mono text-indigo-400 font-semibold">
                              {(recipe.stat_col_similarity_threshold ?? 0.35).toFixed(2)}
                            </span>
                          </div>
                          <input
                            type="range"
                            min="0.10"
                            max="0.90"
                            step="0.05"
                            value={recipe.stat_col_similarity_threshold ?? 0.35}
                            onChange={(e) =>
                              update('stat_col_similarity_threshold', parseFloat(e.target.value))
                            }
                            className="w-full accent-indigo-500"
                          />
                        </div>
                      )}

                      <div>
                        <div className="flex justify-between text-[11px] mb-1">
                          <span className="text-[#8b949e]">Smoothing Window Size</span>
                          <span className="font-mono text-indigo-400 font-semibold">
                            {recipe.stat_col_window_size ?? 3} bp
                          </span>
                        </div>
                        <input
                          type="range"
                          min="1"
                          max="9"
                          step="2"
                          value={recipe.stat_col_window_size ?? 3}
                          onChange={(e) =>
                            update('stat_col_window_size', parseInt(e.target.value))
                          }
                          className="w-full accent-indigo-500"
                        />
                      </div>
                    </>
                  )}

                  {recipe.stat_col_method === 'gblocksblocks' && (
                    <>
                      <div>
                        <div className="flex justify-between text-[11px] mb-1">
                          <span className="text-[#8b949e]">Min Conserved Block Length</span>
                          <span className="font-mono text-indigo-400 font-semibold">
                            {recipe.stat_col_min_block_length ?? 5} bp
                          </span>
                        </div>
                        <input
                          type="range"
                          min="3"
                          max="20"
                          step="1"
                          value={recipe.stat_col_min_block_length ?? 5}
                          onChange={(e) =>
                            update('stat_col_min_block_length', parseInt(e.target.value))
                          }
                          className="w-full accent-indigo-500"
                        />
                        <div className="text-[10px] text-[#8b949e] mt-0.5">
                          Removes isolated conserved spikes &lt; {recipe.stat_col_min_block_length ?? 5} bp
                        </div>
                      </div>

                      <div>
                        <label className="text-[11px] text-[#8b949e] block mb-1">Gap Allowance in Blocks</label>
                        <select
                          value={recipe.stat_col_gap_treatment || 'half'}
                          onChange={(e) =>
                            update('stat_col_gap_treatment', e.target.value as StatColGapTreatment)
                          }
                          className="w-full bg-[#1b2029] border border-[#2d3545] rounded px-2 py-1 text-xs text-[#c9d1d9] outline-none"
                        >
                          <option value="none">None (0% Gaps Allowed)</option>
                          <option value="half">With Half (&le; 50% Gaps Allowed)</option>
                          <option value="all">All (Gaps Allowed in Blocks)</option>
                        </select>
                      </div>
                    </>
                  )}

                  {recipe.stat_col_method === 'entropy' && (
                    <>
                      <div>
                        <div className="flex justify-between text-[11px] mb-1">
                          <span className="text-[#8b949e]">Max Column Shannon Entropy</span>
                          <span className="font-mono text-indigo-400 font-semibold">
                            {(recipe.stat_col_entropy_threshold ?? 1.5).toFixed(2)} bits
                          </span>
                        </div>
                        <input
                          type="range"
                          min="0.4"
                          max="2.0"
                          step="0.1"
                          value={recipe.stat_col_entropy_threshold ?? 1.5}
                          onChange={(e) =>
                            update('stat_col_entropy_threshold', parseFloat(e.target.value))
                          }
                          className="w-full accent-indigo-500"
                        />
                        <div className="text-[10px] text-[#8b949e] mt-0.5">
                          Max theoretical randomness is 2.0 bits (drops high-entropy noise)
                        </div>
                      </div>

                      <div>
                        <div className="flex justify-between text-[11px] mb-1">
                          <span className="text-[#8b949e]">Smoothing Window Size</span>
                          <span className="font-mono text-indigo-400 font-semibold">
                            {recipe.stat_col_window_size ?? 3} bp
                          </span>
                        </div>
                        <input
                          type="range"
                          min="1"
                          max="9"
                          step="2"
                          value={recipe.stat_col_window_size ?? 3}
                          onChange={(e) =>
                            update('stat_col_window_size', parseInt(e.target.value))
                          }
                          className="w-full accent-indigo-500"
                        />
                      </div>
                    </>
                  )}
                </>
              )}
            </div>
          )}
        </div>

        {/* Display 1: Sample Coverage */}
        <div
          className="shrink-0 border border-purple-500/30 rounded-lg bg-purple-500/[0.035] overflow-hidden"
          style={{ order: 1 }}
        >
          <button
            onClick={() => toggleSection('coverage')}
            className="w-full px-3 py-2 bg-purple-500/[0.08] flex items-center justify-between text-purple-100 font-medium"
          >
            <div className="flex items-center gap-2">
              <span className="w-5 h-5 rounded flex items-center justify-center text-purple-400"><Filter className="w-3.5 h-3.5" /></span>
              <span className="text-purple-200">Sample Filter Criteria</span>
            </div>
            {collapsedSections['coverage'] ? (
              <ChevronRight className="w-3.5 h-3.5 text-[#8b949e]" />
            ) : (
              <ChevronDown className="w-3.5 h-3.5 text-[#8b949e]" />
            )}
          </button>

          {!collapsedSections['coverage'] && (
            <div className="p-3 space-y-3">
              <div className="space-y-2.5">
                <label className="flex items-center justify-between cursor-pointer">
                  <span className="text-[#c9d1d9]">Drop samples below</span>
                  <input
                    type="checkbox"
                    checked={sampleOccupancyDraft > 0}
                    onChange={(e) => toggleSampleOccupancy(e.target.checked)}
                    className="rounded bg-[#1f242e] border-[#2d3545] text-purple-500 focus:ring-0"
                  />
                </label>

                {sampleOccupancyDraft > 0 && (
                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Minimum Locus Occupancy</span>
                      <span className="font-mono text-purple-400 font-semibold">
                        {sampleOccupancyDraft}%
                      </span>
                    </div>
                    <input
                      type="range"
                      min="1"
                      max="100"
                      step="1"
                      value={sampleOccupancyDraft}
                      onChange={(e) => setSampleOccupancyDraft(parseFloat(e.target.value))}
                      onPointerUp={(e) =>
                        commitSampleOccupancy(parseFloat(e.currentTarget.value))
                      }
                      onKeyUp={(e) =>
                        commitSampleOccupancy(parseFloat(e.currentTarget.value))
                      }
                      onBlur={(e) =>
                        commitSampleOccupancy(parseFloat(e.currentTarget.value))
                      }
                      className="w-full accent-purple-500"
                    />
                    <p className="mt-1 text-[10px] leading-snug text-[#8b949e]">
                      Below {sampleOccupancyDraft}% of loci → remove that sample from every alignment.
                    </p>
                  </div>
                )}
              </div>

              <div className="space-y-2.5">
                <label className="flex items-center justify-between cursor-pointer">
                  <span className="text-[#c9d1d9]">Filter Divergent Outliers</span>
                  <input
                    type="checkbox"
                    checked={recipe.trim_similarity}
                    onChange={(e) => update('trim_similarity', e.target.checked)}
                    className="rounded bg-[#1f242e] border-[#2d3545] text-purple-500 focus:ring-0"
                  />
                </label>

                {recipe.trim_similarity && (
                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Max Distance to Consensus</span>
                      <span className="font-mono text-purple-400 font-semibold">
                        {(recipe.similarity_threshold * 100).toFixed(0)}%
                      </span>
                    </div>
                    <input
                      type="range"
                      min="0.10"
                      max="0.80"
                      step="0.01"
                      value={recipe.similarity_threshold}
                      onChange={(e) =>
                        update('similarity_threshold', parseFloat(e.target.value))
                      }
                      className="w-full accent-purple-500"
                    />
                  </div>
                )}
              </div>

              <label className="flex items-center justify-between cursor-pointer">
                <span className="text-[#c9d1d9]">Drop low-coverage rows per locus</span>
                <input
                  type="checkbox"
                  checked={recipe.trim_coverage}
                  onChange={(e) => update('trim_coverage', e.target.checked)}
                  className="rounded bg-[#1f242e] border-[#2d3545] text-purple-500 focus:ring-0"
                />
              </label>

              {recipe.trim_coverage && (
                <>
                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Minimum Informative BP</span>
                      <span className="font-mono text-purple-400 font-semibold">
                        {recipe.min_coverage_bp} bp
                      </span>
                    </div>
                    <input
                      type="range"
                      min="0"
                      max="500"
                      step="10"
                      value={recipe.min_coverage_bp}
                      onChange={(e) => update('min_coverage_bp', parseInt(e.target.value))}
                      className="w-full accent-purple-500"
                    />
                  </div>

                  <div>
                    <div className="flex justify-between text-[11px] mb-1">
                      <span className="text-[#8b949e]">Min % Length Relative</span>
                      <span className="font-mono text-purple-400 font-semibold">
                        {recipe.min_coverage_percent}%
                      </span>
                    </div>
                    <input
                      type="range"
                      min="10"
                      max="90"
                      step="5"
                      value={recipe.min_coverage_percent}
                      onChange={(e) =>
                        update('min_coverage_percent', parseFloat(e.target.value))
                      }
                      className="w-full accent-purple-500"
                    />
                  </div>

                  <div className="pt-1">
                    <label className="text-[11px] text-[#8b949e] block mb-1">Relative To</label>
                    <select
                      value={recipe.relative_width}
                      onChange={(e) =>
                        update('relative_width', e.target.value as RelativeWidth)
                      }
                      className="w-full bg-[#1b2029] border border-[#2d3545] rounded px-2 py-1 text-xs text-[#c9d1d9] outline-none"
                    >
                      <option value="sample">Longest Sample in Locus</option>
                      <option value="alignment">Total Alignment Length</option>
                    </select>
                  </div>
                </>
              )}
            </div>
          )}
        </div>
      </div>
    </aside>
  );
};
