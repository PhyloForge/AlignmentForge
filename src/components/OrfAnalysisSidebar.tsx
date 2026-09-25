import React, { useRef } from 'react';
import { Settings2, ShieldCheck, Palette, Target, Play, Upload, Trash2, Dna, RotateCcw, ChevronDown, ChevronRight } from 'lucide-react';
import { 
  TrimmingRecipe, 
  GeneticCode, 
  AminoAcidViewerSettings, 
  AminoAcidColorScheme,
  OrfSearchMode
} from '../types';

interface OrfAnalysisSidebarProps {
  recipe: TrimmingRecipe;
  onChangeRecipe: (recipe: TrimmingRecipe) => void;
  onResetRecipe: () => void;
  aminoAcidViewerSettings: AminoAcidViewerSettings;
  onChangeAminoAcidViewerSettings: (settings: AminoAcidViewerSettings) => void;
}

/**
 * Upper bound for a count slider. A configuration file can hold a value above
 * the usual range, and a fixed `max` would snap it down as soon as the user
 * touched the control. Widening the track to fit keeps the loaded value intact.
 */
const sliderMax = (standardMax: number, currentValue: number | undefined): number =>
  Math.max(standardMax, Math.ceil(currentValue ?? 0));

export function OrfAnalysisSidebar({
  recipe,
  onChangeRecipe,
  onResetRecipe,
  aminoAcidViewerSettings,
  onChangeAminoAcidViewerSettings
}: OrfAnalysisSidebarProps) {
  // Sections start closed, matching the filtering sidebar.
  const [collapsedSections, setCollapsedSections] = React.useState<Record<string, boolean>>({
    analysis_modes: true,
    reference_guided_exons_introns: true,
    codon_qc: true,
    amino_acid_viewer: true,
  });

  const toggleSection = (key: string) =>
    setCollapsedSections((current) => ({ ...current, [key]: !current[key] }));

  const update = <K extends keyof TrimmingRecipe>(key: K, value: TrimmingRecipe[K]) => {
    onChangeRecipe({ ...recipe, [key]: value });
  };
  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleReferenceUpload = async (files: FileList | null) => {
    if (!files || files.length === 0) return;
    try {
      const newRefs: Record<string, string> = { ...recipe.orf_reference_sequences };
      for (let i = 0; i < files.length; i++) {
        const file = files[i];
        const text = await file.text();
        const lines = text.split('\n');
        let currentHeader = '';
        let currentSeq = '';

        for (const line of lines) {
          const trimmed = line.trim();
          if (trimmed.startsWith('>')) {
            if (currentHeader && currentSeq) {
              newRefs[currentHeader] = currentSeq;
            }
            currentHeader = trimmed.substring(1).trim().split(' ')[0];
            currentSeq = '';
          } else if (trimmed) {
            currentSeq += trimmed;
          }
        }
        if (currentHeader && currentSeq) {
          newRefs[currentHeader] = currentSeq;
        }
      }
      onChangeRecipe({
        ...recipe,
        orf_reference_sequences: newRefs,
        orf_use_references: true
      });
    } catch (err) {
      console.error('Failed to parse references', err);
    }
    
    if (fileInputRef.current) {
      fileInputRef.current.value = '';
    }
  };

  const refCount = Object.keys(recipe.orf_reference_sequences || {}).length;

  return (
    <div className="w-80 min-h-0 flex-none bg-[#14171d] border-l border-[#232833] flex flex-col h-full overflow-hidden select-none">
      <div className="flex flex-col h-full">
        
        {/* Sidebar Header */}
        <div className="flex-none p-3.5 border-b border-[#232833] bg-[#171b22]/70 flex items-center justify-between">
          <div className="flex items-center gap-2">
            <Dna className="w-4 h-4 text-pink-500" />
            <span className="font-semibold text-xs text-[#dce6ff] uppercase tracking-wider">
              ORF Configuration
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
        
        <div className="flex-1 overflow-y-auto p-4 space-y-6">
          {!recipe.enable_orf ? (
            <div className="shrink-0 border border-emerald-500/25 rounded-lg bg-emerald-500/5 p-4 text-center">
              <h3 className="mb-2 font-semibold text-emerald-100 text-sm">ORF Analysis is Off</h3>
              <p className="mb-4 text-[11px] text-[#8b949e]">
                Enable ORF Analysis to configure and find reading frames.
              </p>
              <button
                type="button"
                onClick={() => update('enable_orf', true)}
                className="inline-flex items-center justify-center gap-1.5 rounded bg-emerald-500/10 border border-emerald-500/20 px-3 py-1.5 text-xs font-medium text-emerald-400 transition-colors hover:bg-emerald-500/20"
              >
                <Play className="h-3.5 w-3.5" />
                Turn on ORF Analysis
              </button>
            </div>
          ) : (
            <div className="space-y-6">
              
              <div className="shrink-0 border border-purple-500/30 rounded-lg bg-purple-500/[0.035] overflow-hidden">
                <button
                  type="button"
                  onClick={() => toggleSection('analysis_modes')}
                  className="w-full px-3 py-2 bg-purple-500/[0.08] flex items-center justify-between text-purple-200 font-medium text-xs">
                  <div className="flex items-center gap-2">
                    <Settings2 className="w-4 h-4 text-purple-400" />
                    <span className="font-semibold text-purple-300">Analysis Modes</span>
                  </div>
                  {collapsedSections['analysis_modes'] ? (
                    <ChevronRight className="w-3.5 h-3.5 text-purple-400/70" />
                  ) : (
                    <ChevronDown className="w-3.5 h-3.5 text-purple-400/70" />
                  )}
                </button>
                {!collapsedSections['analysis_modes'] && (
                  <div className="p-3 space-y-3">
                  <div>
                    <label className="text-[11px] text-[#8b949e] block mb-1">ORF Search Mode</label>
                    <select
                      value={recipe.orf_search_mode}
                      onChange={(event) => update('orf_search_mode', event.target.value as OrfSearchMode)}
                      className="w-full bg-[#1b2029] border border-[#2d3545] rounded px-2 py-1 text-xs text-[#c9d1d9] outline-none"
                    >
                      <option value="continuouscds">Continuous CDS (Whole alignment)</option>
                      <option value="bestsharedsegment">Candidate ORF Extraction (Heuristic/Shared Segment)</option>
                      <option value="referenceguided" disabled={!recipe.orf_use_references || refCount === 0}>Reference-Guided CDS (Requires FASTA)</option>
                      <option value="referencecandidateorf" disabled={!recipe.orf_use_references || refCount === 0}>Reference-Guided Candidate ORF (Requires FASTA)</option>
                    </select>
                    <p className="mt-2 text-[10px] leading-snug text-[#8b949e]">
                      Sets the method for locating reading frames. Heuristic/Candidate extraction rapidly checks length and start/stop rules. Continuous CDS assumes the whole alignment is coding. Reference-guided uses your provided sequences.
                    </p>
                  </div>
                  
                </div>
                )}
              </div>

              <div className="shrink-0 border border-pink-500/30 rounded-lg bg-pink-500/[0.035] overflow-hidden">
                <button
                  type="button"
                  onClick={() => toggleSection('reference_guided_exons_introns')}
                  className="w-full px-3 py-2 bg-pink-500/[0.08] flex items-center justify-between text-pink-200 font-medium text-xs">
                  <div className="flex items-center gap-2">
                    <Target className="w-4 h-4 text-pink-400" />
                    <span className="font-semibold text-pink-300">Reference-Guided Exons & Introns</span>
                  </div>
                  {collapsedSections['reference_guided_exons_introns'] ? (
                    <ChevronRight className="w-3.5 h-3.5 text-pink-400/70" />
                  ) : (
                    <ChevronDown className="w-3.5 h-3.5 text-pink-400/70" />
                  )}
                </button>
                {!collapsedSections['reference_guided_exons_introns'] && (
                  <div className="p-3 space-y-3">
                  <input 
                    type="file" 
                    ref={fileInputRef}
                    multiple={true}
                    accept=".fa,.fasta,.fna" 
                    className="hidden" 
                    onChange={(e) => handleReferenceUpload(e.target.files)} 
                  />
                  <div className="grid grid-cols-2 gap-2">
                    <button
                      type="button"
                      onClick={() => fileInputRef.current?.click()}
                      className="inline-flex items-center justify-center gap-1.5 rounded border border-[#3a4250] bg-[#1f242e] px-2 py-1.5 text-[11px] font-medium text-[#c9d1d9] hover:bg-[#2a313d]"
                    >
                      <Upload className="w-3 h-3" /> Load FASTA
                    </button>
                    <button
                      type="button"
                      disabled={refCount === 0}
                      onClick={() => {
                        onChangeRecipe({
                          ...recipe,
                          orf_use_references: false,
                          orf_reference_sequences: {},
                          orf_search_mode: (recipe.orf_search_mode === 'referencecandidateorf' || recipe.orf_search_mode === 'referenceguided') ? 'continuouscds' as OrfSearchMode : recipe.orf_search_mode
                        });
                      }}
                      className="inline-flex items-center justify-center gap-1.5 rounded border border-rose-500/25 bg-rose-500/10 px-2 py-1.5 text-[11px] text-rose-300 disabled:opacity-35"
                    >
                      <Trash2 className="w-3 h-3" /> Clear
                    </button>
                  </div>
                  <div className="flex items-center justify-between rounded bg-[#101319] px-2 py-1.5 text-[11px]">
                    <span className="text-[#8b949e]">Loaded reference loci</span>
                    <span className="font-mono font-semibold text-pink-300">{refCount}</span>
                  </div>
                  <p className="text-[10px] leading-snug text-[#8b949e]">
                    Load reference FASTA sequences here to guide intron boundaries and frame translations.
                  </p>
                </div>
                )}
              </div>

              <div className="shrink-0 border border-teal-500/30 rounded-lg bg-teal-500/[0.035] overflow-hidden">
                <button
                  type="button"
                  onClick={() => toggleSection('codon_qc')}
                  className="w-full px-3 py-2 bg-teal-500/[0.08] flex items-center justify-between text-teal-200 font-medium text-xs">
                  <div className="flex items-center gap-2">
                    <ShieldCheck className="w-4 h-4 text-teal-400" />
                    <span className="font-semibold text-teal-300">Codon QC</span>
                  </div>
                  {collapsedSections['codon_qc'] ? (
                    <ChevronRight className="w-3.5 h-3.5 text-teal-400/70" />
                  ) : (
                    <ChevronDown className="w-3.5 h-3.5 text-teal-400/70" />
                  )}
                </button>
                {!collapsedSections['codon_qc'] && (
                  <div className="p-3 space-y-3">
                  
                  <div>
                    <label className="text-[11px] text-[#8b949e] block mb-1">Genetic code</label>
                    <select
                      value={recipe.genetic_code}
                      onChange={(event) => update('genetic_code', event.target.value as GeneticCode)}
                      className="w-full bg-[#1b2029] border border-[#2d3545] rounded px-2 py-1 text-xs text-[#c9d1d9] outline-none"
                    >
                      <option value="standard">Standard nuclear</option>
                      <option value="vertebratemitochondrial">Vertebrate mitochondrial</option>
                      <option value="invertebratemitochondrial">Invertebrate mitochondrial</option>
                    </select>
                  </div>

                  <div className="pt-2 mt-3 border-t border-[#2d3545]">
                    <label className="flex items-center justify-between cursor-pointer mb-3 mt-3">
                      <span className="text-[#c9d1d9] text-[11px]">Trim Terminal Frameshifts (!)</span>
                      <input
                        type="checkbox"
                        checked={recipe.macse_trim_terminal ?? true}
                        onChange={(e) => update('macse_trim_terminal', e.target.checked)}
                        className="rounded bg-[#1f242e] border-[#2d3545] text-teal-500 focus:ring-0"
                      />
                    </label>

                    <div className="space-y-4">
                      <div>
                        <div className="flex justify-between text-[11px] mb-1">
                          <span className="text-[#8b949e]">Max Frameshifts / Sample</span>
                          <span className="font-mono text-teal-400 font-semibold">{recipe.macse_max_internal_sample ?? 3}</span>
                        </div>
                        <input
                          type="range"
                          min="0"
                          max={sliderMax(20, recipe.macse_max_internal_sample ?? 3)}
                          step="1"
                          value={recipe.macse_max_internal_sample ?? 3}
                          onChange={(e) => update('macse_max_internal_sample', parseInt(e.target.value))}
                          className="w-full accent-teal-500"
                        />
                      </div>

                      <div>
                        <div className="flex justify-between text-[11px] mb-1">
                          <span className="text-[#8b949e]">Max Frameshift Cols / Locus</span>
                          <span className="font-mono text-teal-400 font-semibold">{recipe.macse_max_internal_locus ?? 10}</span>
                        </div>
                        <input
                          type="range"
                          min="0"
                          max={sliderMax(50, recipe.macse_max_internal_locus ?? 10)}
                          step="1"
                          value={recipe.macse_max_internal_locus ?? 10}
                          onChange={(e) => update('macse_max_internal_locus', parseInt(e.target.value))}
                          className="w-full accent-teal-500"
                        />
                      </div>

                      <div className="pt-2 border-t border-[#2d3545]/50">
                        <div className="mb-4 mt-2">
                          <label className="flex items-center justify-between cursor-pointer mb-1">
                            <span className="text-[#c9d1d9] text-[11px]">Mask Stop Codons</span>
                            <input
                              type="checkbox"
                              checked={recipe.stop_codon_action === 'maskcodon'}
                              onChange={(e) => update('stop_codon_action', e.target.checked ? 'maskcodon' : 'removesample')}
                              className="rounded bg-[#1f242e] border-[#2d3545] text-teal-500 focus:ring-0"
                            />
                          </label>

                        </div>
                        <div className="flex justify-between text-[11px] mb-1">
                          <span className="text-[#8b949e]">Max Stop Codons / Sample</span>
                          <span className="font-mono text-teal-400 font-semibold">{recipe.max_stop_codons_sample ?? 2}</span>
                        </div>
                        <input
                          type="range"
                          min="0"
                          max={sliderMax(10, recipe.max_stop_codons_sample ?? 2)}
                          step="1"
                          value={recipe.max_stop_codons_sample ?? 2}
                          onChange={(e) => update('max_stop_codons_sample', parseInt(e.target.value))}
                          className="w-full accent-teal-500"
                        />
                      </div>

                      <div>
                        <div className="flex justify-between text-[11px] mb-1">
                          <span className="text-[#8b949e]">Max Stop Codons / Locus</span>
                          <span className="font-mono text-teal-400 font-semibold">{recipe.max_stop_codons_locus ?? 5}</span>
                        </div>
                        <input
                          type="range"
                          min="0"
                          max={sliderMax(20, recipe.max_stop_codons_locus ?? 5)}
                          step="1"
                          value={recipe.max_stop_codons_locus ?? 5}
                          onChange={(e) => update('max_stop_codons_locus', parseInt(e.target.value))}
                          className="w-full accent-teal-500"
                        />
                      </div>
                    </div>
                  </div>
                </div>
                )}
              </div>

              <div className="shrink-0 border border-orange-500/30 rounded-lg bg-orange-500/[0.035] overflow-hidden">
                <button
                  type="button"
                  onClick={() => toggleSection('amino_acid_viewer')}
                  className="w-full px-3 py-2 bg-orange-500/[0.08] flex items-center justify-between text-orange-200 font-medium text-xs">
                  <div className="flex items-center gap-2">
                    <Palette className="w-4 h-4 text-orange-400" />
                    <span className="font-semibold text-orange-300">Amino Acid Viewer</span>
                  </div>
                  {collapsedSections['amino_acid_viewer'] ? (
                    <ChevronRight className="w-3.5 h-3.5 text-orange-400/70" />
                  ) : (
                    <ChevronDown className="w-3.5 h-3.5 text-orange-400/70" />
                  )}
                </button>
                {!collapsedSections['amino_acid_viewer'] && (
                  <div className="p-3 space-y-3">
                  <button
                    type="button"
                    onClick={() => onChangeAminoAcidViewerSettings({
                      ...aminoAcidViewerSettings,
                      enabled: !aminoAcidViewerSettings.enabled
                    })}
                    className={`w-full rounded border px-2.5 py-2 text-[11px] font-semibold transition-colors ${
                      aminoAcidViewerSettings.enabled
                        ? 'border-[#4a5361] bg-[#2a313d] text-[#dce6ff] hover:bg-[#333c4a]'
                        : 'border-[#3a4250] bg-[#1f242e] text-[#c9d1d9] hover:bg-[#2a313d]'
                    }`}
                  >
                    {aminoAcidViewerSettings.enabled ? 'Return Viewer to Nucleotides' : 'Convert Viewer to Amino Acids'}
                  </button>
                  
                  {aminoAcidViewerSettings.enabled && (
                    <div className="space-y-2.5 border-t border-orange-500/15 pt-2.5">
                      <div>
                        <label className="text-[11px] text-[#8b949e] block mb-1">Color palette</label>
                        <select
                          value={aminoAcidViewerSettings.colorScheme}
                          onChange={(event) => onChangeAminoAcidViewerSettings({
                            ...aminoAcidViewerSettings,
                            colorScheme: event.target.value as AminoAcidColorScheme
                          })}
                          className="w-full bg-[#1b2029] border border-[#2d3545] rounded px-2 py-1 text-xs text-[#c9d1d9] outline-none"
                        >
                          <option value="chemistry">Chemistry</option>
                          <option value="colorblind">Color-blind Safe</option>
                          <option value="zappo">Zappo</option>
                          <option value="taylor">Taylor</option>
                          <option value="hydrophobicity">Hydrophobicity</option>
                          <option value="monochrome">Monochrome</option>
                        </select>
                      </div>
                      <label className="flex items-center justify-between cursor-pointer">
                        <span className="text-[#c9d1d9] text-[11px]">Dim consensus matches</span>
                        <input
                          type="checkbox"
                          checked={aminoAcidViewerSettings.dimConsensusMatches}
                          onChange={(event) => onChangeAminoAcidViewerSettings({
                            ...aminoAcidViewerSettings,
                            dimConsensusMatches: event.target.checked,
                          })}
                          className="rounded bg-[#1f242e] border-[#2d3545] text-orange-500 focus:ring-0"
                        />
                      </label>
                    </div>
                  )}
                </div>
                )}
              </div>

            </div>
          )}
        </div>
      </div>
    </div>
  );
}
