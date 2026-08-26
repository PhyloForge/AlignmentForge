import React, { useState } from 'react';
import {
  X,
  Download,
  Layers,
  FolderOpen,
  CheckCircle,
  FileSpreadsheet,
  FileCode,
  FileText,
} from 'lucide-react';
import {
  AlignmentFormat,
  BatchExportConfig,
  BatchExportResult,
  ConcatenateConfig,
  ConcatenateResult,
  GroupedConcatenateConfig,
  GroupedConcatenateResult,
  TrimmingRecipe,
} from '../types';
import { openSaveDirectoryDialog, openFileDialog, runBatchExport, runConcatenate, runGroupedConcatenate } from '../tauriClient';
import { Sparkles, FileSearch } from 'lucide-react';

interface ExportModalProps {
  mode: 'batch' | 'concatenate' | 'group';
  onClose: () => void;
  selectedPaths: string[];
  allPaths: string[];
  recipe: TrimmingRecipe;
}

const joinDirectory = (parent: string, child: string) => {
  if (!parent) return child;
  return `${parent}${/[\\/]$/.test(parent) ? '' : '/'}${child}`;
};

const validateFolderName = (label: string, value: string) => {
  const name = value.trim();
  if (!name || name === '.' || name === '..' || /[\\/]/.test(name)) {
    throw new Error(`${label} must be one folder name without / or \\.`);
  }
  return name;
};

export const ExportModal: React.FC<ExportModalProps> = ({
  mode,
  onClose,
  selectedPaths,
  allPaths,
  recipe,
}) => {
  const [outputDir, setOutputDir] = useState<string>('output/trimmed_alignments');
  const [batchParentDir, setBatchParentDir] = useState<string>('output');
  const [batchOutputFolderName, setBatchOutputFolderName] = useState<string>('trimmed_alignments');
  const [generalAlignmentFolderName, setGeneralAlignmentFolderName] = useState<string>('all_alignments');
  const [orfAlignmentFolderName, setOrfAlignmentFolderName] = useState<string>('orf_alignments');
  const [intronFolderName, setIntronFolderName] = useState<string>('intron_alignments');
  const [outputPrefix, setOutputPrefix] = useState<string>('output/supermatrix');
  const [geneMappingPath, setGeneMappingPath] = useState<string>('');
  const [outputFormat, setOutputFormat] = useState<AlignmentFormat>('phylip');
  const [onlyPassing, setOnlyPassing] = useState<boolean>(true);
  const [exportGeneralAlignments, setExportGeneralAlignments] = useState<boolean>(true);
  const [exportOrfAlignments, setExportOrfAlignments] = useState<boolean>(recipe.enable_orf);
  const [exportScope, setExportScope] = useState<'all' | 'selected'>('all');

  const [saveSummaryCsv, setSaveSummaryCsv] = useState<boolean>(true);
  const [saveRecipeJson, setSaveRecipeJson] = useState<boolean>(true);
  const [exportIntrons, setExportIntrons] = useState<boolean>(
    recipe.enable_orf && recipe.orf_use_references && Object.keys(recipe.orf_reference_sequences ?? {}).length > 0
  );
  const [writeRaxmlPartitions, setWriteRaxmlPartitions] = useState<boolean>(true);
  const [writeNexusPartitions, setWriteNexusPartitions] = useState<boolean>(true);

  const [isRunning, setIsRunning] = useState<boolean>(false);
  const [batchResult, setBatchResult] = useState<BatchExportResult | null>(null);
  const [concatResult, setConcatResult] = useState<ConcatenateResult | null>(null);
  const [groupResult, setGroupResult] = useState<GroupedConcatenateResult | null>(null);
  const [errorMsg, setErrorMsg] = useState<string | null>(null);

  const activePaths =
    exportScope === 'selected' && selectedPaths.length > 0 ? selectedPaths : allPaths;

  const batchOutputDir = joinDirectory(batchParentDir, batchOutputFolderName.trim());

  const handlePickOutputDir = async () => {
    const picked = await openSaveDirectoryDialog();
    if (picked) {
      if (mode === 'batch') {
        setBatchParentDir(picked);
      } else {
        setOutputDir(picked);
        setOutputPrefix(`${picked}/supermatrix`);
      }
    }
  };

  const handlePickMappingFile = async () => {
    const picked = await openFileDialog();
    if (picked) {
      setGeneMappingPath(picked);
    }
  };

  const handleExecuteExport = async () => {
    setIsRunning(true);
    setErrorMsg(null);

    try {
      if (mode === 'batch') {
        const outputFolderName = validateFolderName('Output folder name', batchOutputFolderName);
        const generalFolderName = validateFolderName(
          'General alignment folder name',
          generalAlignmentFolderName
        );
        const orfFolderName = validateFolderName('ORF alignment folder name', orfAlignmentFolderName);
        const intronName = validateFolderName('Intron folder name', intronFolderName);
        const activeFolderNames = [
          exportGeneralAlignments ? generalFolderName : null,
          exportOrfAlignments && recipe.enable_orf ? orfFolderName : null,
          exportIntrons ? intronName : null,
        ].filter((name): name is string => name !== null);
        if (new Set(activeFolderNames.map((name) => name.toLocaleLowerCase())).size !== activeFolderNames.length) {
          throw new Error('Exported alignment folders must have different names.');
        }
        const config: BatchExportConfig = {
          input_paths: activePaths,
          output_directory: joinDirectory(batchParentDir, outputFolderName),
          general_alignment_directory_name: generalFolderName,
          orf_alignment_directory_name: orfFolderName,
          intron_directory_name: intronName,
          output_format: outputFormat,
          only_passing: onlyPassing,
          export_general_alignments: exportGeneralAlignments,
          export_orf_alignments: exportOrfAlignments,
          save_recipe_json: saveRecipeJson,
          save_summary_csv: saveSummaryCsv,
          export_introns: exportIntrons,
        };
        const res = await runBatchExport(config, recipe);
        setBatchResult(res);
      } else if (mode === 'concatenate') {
        const config: ConcatenateConfig = {
          input_paths: activePaths,
          output_file_prefix: outputPrefix,
          output_format: outputFormat,
          only_passing: onlyPassing,
          write_raxml_partitions: writeRaxmlPartitions,
          write_nexus_partitions: writeNexusPartitions,
        };
        const res = await runConcatenate(config, recipe);
        setConcatResult(res);
      } else if (mode === 'group') {
        if (!geneMappingPath) {
          throw new Error('Please select a gene mapping file.');
        }
        const config: GroupedConcatenateConfig = {
          input_paths: activePaths,
          output_directory: outputDir,
          gene_mapping_csv_path: geneMappingPath,
          output_format: outputFormat,
          only_passing: onlyPassing,
          write_raxml_partitions: writeRaxmlPartitions,
          write_nexus_partitions: writeNexusPartitions,
        };
        const res = await runGroupedConcatenate(config, recipe);
        setGroupResult(res);
      }
    } catch (e: any) {
      setErrorMsg(e.toString());
    } finally {
      setIsRunning(false);
    }
  };

  return (
    <div className="fixed inset-0 z-50 bg-black/70 backdrop-blur-sm flex items-center justify-center p-4 select-none">
      <div className="bg-[#14171d] border border-[#232833] rounded-xl w-[540px] max-h-[90vh] overflow-hidden shadow-2xl flex flex-col text-[#c9d1d9]">
        {/* Header */}
        <div className="px-5 py-3.5 border-b border-[#232833] flex items-center justify-between bg-[#171b22]">
          <div className="flex items-center gap-2 font-semibold text-sm text-[#dce6ff]">
            {mode === 'batch' ? (
              <>
                <Download className="w-4 h-4 text-blue-400" />
                <span>Batch Export Trimmed Alignments</span>
              </>
            ) : mode === 'concatenate' ? (
              <>
                <Layers className="w-4 h-4 text-emerald-400" />
                <span>Concatenate into Supermatrix & Partitions</span>
              </>
            ) : (
              <>
                <Sparkles className="w-4 h-4 text-fuchsia-400" />
                <span>Batch Concatenate Exons by Gene</span>
              </>
            )}
          </div>
          <button
            onClick={onClose}
            className="p-1 rounded-md text-[#8b949e] hover:text-[#c9d1d9] hover:bg-[#232833]"
          >
            <X className="w-4 h-4" />
          </button>
        </div>

        {/* Content Body */}
        <div className="p-5 space-y-4 text-xs overflow-y-auto">
          {batchResult || concatResult || groupResult ? (
            /* Success Feedback */
            <div className="py-4 text-center space-y-3">
              <div className="w-12 h-12 rounded-full bg-emerald-500/20 text-emerald-400 flex items-center justify-center mx-auto">
                <CheckCircle className="w-6 h-6" />
              </div>
              <h3 className="font-semibold text-sm text-[#dce6ff]">
                {mode === 'batch' 
                  ? 'Batch Export Complete!' 
                  : mode === 'concatenate' 
                  ? 'Supermatrix Assembly Complete!' 
                  : 'Gene Batch Complete!'}
              </h3>
              {batchResult && (
                <div className="font-mono text-xs text-[#8b949e] bg-[#0e1014] p-3 rounded border border-[#232833] text-left space-y-1">
                  <div>General alignments exported: <b className="text-emerald-400">{batchResult.total_exported}</b></div>
                  <div>Catalog QC failures: <b className="text-rose-400">{batchResult.total_discarded}</b></div>
                  {(batchResult.total_orfs_exported ?? 0) > 0 && (
                    <div>Accepted ORF alignments exported: <b className="text-violet-300">{batchResult.total_orfs_exported}</b></div>
                  )}
                  {(batchResult.total_introns_exported ?? 0) > 0 && (
                    <div>Exported intron alignments: <b className="text-blue-300">{batchResult.total_introns_exported}</b></div>
                  )}
                  {batchResult.intron_directory_path && (
                    <div className="text-[11px] text-blue-400 pt-1">Intron folder: {batchResult.intron_directory_path}</div>
                  )}
                  {batchResult.alignment_directory_path && (
                    <div className="text-[11px] text-emerald-400 pt-1">General alignment folder: {batchResult.alignment_directory_path}</div>
                  )}
                  {batchResult.orf_directory_path && (
                    <div className="text-[11px] text-violet-400 pt-1">ORF alignment folder: {batchResult.orf_directory_path}</div>
                  )}
                  {batchResult.summary_csv_path && (
                    <div className="text-[11px] text-blue-400 pt-1">
                      📄 Summary CSV: {batchResult.summary_csv_path}
                    </div>
                  )}
                </div>
              )}
              {concatResult && (
                <div className="font-mono text-xs text-[#8b949e] bg-[#0e1014] p-3 rounded border border-[#232833] text-left space-y-1">
                  <div>Total Taxa: <b className="text-emerald-400">{concatResult.total_taxa}</b></div>
                  <div>Total Length: <b className="text-emerald-400">{concatResult.total_length.toLocaleString()} bp</b></div>
                  <div>Concatenated Loci: <b className="text-emerald-400">{concatResult.total_loci}</b></div>
                  <div className="text-[11px] text-cyan-400 pt-1 truncate">
                    💾 Supermatrix: {concatResult.supermatrix_path}
                  </div>
                </div>
              )}
              {groupResult && (
                <div className="font-mono text-xs text-[#8b949e] bg-[#0e1014] p-3 rounded border border-[#232833] text-left space-y-1">
                  <div>Generated Gene Supermatrices: <b className="text-emerald-400">{groupResult.total_genes}</b></div>
                  <div>Total Exons Assembled: <b className="text-emerald-400">{groupResult.total_exons_processed}</b></div>
                  <div className="text-[11px] text-cyan-400 pt-1 truncate">
                    💾 Output Directory: {groupResult.output_directory}
                  </div>
                </div>
              )}
              <div className="pt-2">
                <button
                  onClick={onClose}
                  className="px-5 py-2 bg-blue-600 hover:bg-blue-500 text-white rounded-md font-medium text-xs shadow-md transition-colors"
                >
                  Done
                </button>
              </div>
            </div>
          ) : (
            <>
              {/* Metadata File (Gene Batch only) */}
              {mode === 'group' && (
                <div className="mb-4">
                  <label className="text-[11px] font-semibold text-[#8b949e] uppercase block mb-1.5">
                    Exon to Gene Mapping File (CSV/TSV/TXT)
                  </label>
                  <div className="flex items-center gap-2">
                    <input
                      type="text"
                      readOnly
                      value={geneMappingPath}
                      placeholder="Select metadata document..."
                      className="flex-1 bg-[#1b2029] border border-[#2d3545] rounded-md px-3 py-1.5 font-mono text-xs text-[#c9d1d9] outline-none"
                    />
                    <button
                      onClick={handlePickMappingFile}
                      className="px-3 py-1.5 bg-[#1f242e] hover:bg-[#28303d] border border-[#2d3545] rounded-md text-xs text-[#c9d1d9] flex items-center gap-1.5 whitespace-nowrap"
                    >
                      <FileSearch className="w-3.5 h-3.5 text-fuchsia-400" />
                      <span>Browse File</span>
                    </button>
                  </div>
                </div>
              )}

              {/* Output Directory / Path */}
              <div>
                <label className="text-[11px] font-semibold text-[#8b949e] uppercase block mb-1.5">
                  {mode === 'batch'
                    ? 'Parent Output Destination'
                    : mode === 'group'
                      ? 'Output Folder Destination'
                      : 'Supermatrix Output Path Prefix'}
                </label>
                <div className="flex items-center gap-2">
                  <input
                    type="text"
                    value={mode === 'batch' ? batchParentDir : mode === 'group' ? outputDir : outputPrefix}
                    onChange={(e) =>
                      mode === 'batch'
                        ? setBatchParentDir(e.target.value)
                        : mode === 'group'
                          ? setOutputDir(e.target.value)
                          : setOutputPrefix(e.target.value)
                    }
                    className="flex-1 bg-[#1b2029] border border-[#2d3545] rounded-md px-3 py-1.5 font-mono text-xs text-[#c9d1d9] outline-none"
                  />
                  <button
                    onClick={handlePickOutputDir}
                    className="px-3 py-1.5 bg-[#1f242e] hover:bg-[#28303d] border border-[#2d3545] rounded-md text-xs text-[#c9d1d9] flex items-center gap-1.5"
                  >
                    <FolderOpen className="w-3.5 h-3.5 text-blue-400" />
                    <span>Browse</span>
                  </button>
                </div>
              </div>

              {mode === 'batch' && (
                <div className="rounded-md border border-[#232833] bg-[#11141a] p-3 space-y-2.5">
                  <label className="text-[11px] font-semibold text-[#8b949e] uppercase block">
                    Folder Names
                  </label>
                  <div className="grid grid-cols-2 gap-2.5">
                    <label className="space-y-1">
                      <span className="text-[10px] text-[#8b949e]">Output</span>
                      <input
                        type="text"
                        value={batchOutputFolderName}
                        onChange={(e) => setBatchOutputFolderName(e.target.value)}
                        className="w-full bg-[#1b2029] border border-[#2d3545] rounded-md px-2.5 py-1.5 font-mono text-xs text-[#c9d1d9] outline-none"
                      />
                    </label>
                    <label className="space-y-1">
                      <span className="text-[10px] text-[#8b949e]">General alignments</span>
                      <input
                        type="text"
                        value={generalAlignmentFolderName}
                        onChange={(e) => setGeneralAlignmentFolderName(e.target.value)}
                        className="w-full bg-[#1b2029] border border-[#2d3545] rounded-md px-2.5 py-1.5 font-mono text-xs text-[#c9d1d9] outline-none"
                      />
                    </label>
                    <label className="space-y-1">
                      <span className="text-[10px] text-[#8b949e]">ORF alignments</span>
                      <input
                        type="text"
                        value={orfAlignmentFolderName}
                        onChange={(e) => setOrfAlignmentFolderName(e.target.value)}
                        className="w-full bg-[#1b2029] border border-[#2d3545] rounded-md px-2.5 py-1.5 font-mono text-xs text-[#c9d1d9] outline-none"
                      />
                    </label>
                    <label className="space-y-1">
                      <span className="text-[10px] text-[#8b949e]">Intron alignments</span>
                      <input
                        type="text"
                        value={intronFolderName}
                        onChange={(e) => setIntronFolderName(e.target.value)}
                        className="w-full bg-[#1b2029] border border-[#2d3545] rounded-md px-2.5 py-1.5 font-mono text-xs text-[#c9d1d9] outline-none"
                      />
                    </label>
                  </div>
                  <div className="text-[10px] font-mono text-[#6e7681] break-all">
                    Output: {batchOutputDir || '—'}
                  </div>
                </div>
              )}

              {/* Output Format */}
              <div className="grid grid-cols-2 gap-3">
                <div>
                  <label className="text-[11px] font-semibold text-[#8b949e] uppercase block mb-1.5">
                    Output Format
                  </label>
                  <select
                    value={outputFormat}
                    onChange={(e) => setOutputFormat(e.target.value as AlignmentFormat)}
                    className="w-full bg-[#1b2029] border border-[#2d3545] rounded-md px-2.5 py-1.5 text-xs text-[#c9d1d9] outline-none"
                  >
                    <option value="phylip">Relaxed PHYLIP (.phy)</option>
                    <option value="fasta">FASTA (.fa)</option>
                    <option value="nexus">NEXUS (.nex)</option>
                  </select>
                </div>

                <div>
                  <label className="text-[11px] font-semibold text-[#8b949e] uppercase block mb-1.5">
                    Dataset Scope
                  </label>
                  <select
                    value={exportScope}
                    onChange={(e) => setExportScope(e.target.value as 'all' | 'selected')}
                    className="w-full bg-[#1b2029] border border-[#2d3545] rounded-md px-2.5 py-1.5 text-xs text-[#c9d1d9] outline-none"
                  >
                    <option value="all">All Loci ({allPaths.length})</option>
                    <option value="selected" disabled={selectedPaths.length === 0}>
                      Selected Loci Only ({selectedPaths.length})
                    </option>
                  </select>
                </div>
              </div>

              {/* Checkboxes & Options */}
              <div className="space-y-2 pt-1 border-t border-[#232833]">
                <label className="flex items-center gap-2 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={onlyPassing}
                    onChange={(e) => setOnlyPassing(e.target.checked)}
                    className="rounded bg-[#1f242e] border-[#2d3545] text-blue-500"
                  />
                  <span className="text-[#c9d1d9]">
                    {mode === 'batch'
                      ? `Apply Catalog quality gates to general alignments and introns (${recipe.min_taxa} taxa, ${recipe.min_length} bp, ${recipe.max_gap_percent}% gap). ORF acceptance remains independent.`
                      : `Export only loci passing quality gates (${recipe.min_taxa} taxa, ${recipe.min_length} bp, ${recipe.max_gap_percent}% gap)`}
                  </span>
                </label>

                {mode === 'batch' ? (
                  <>
                    <label className="flex items-center gap-2 cursor-pointer">
                      <input
                        type="checkbox"
                        checked={exportGeneralAlignments}
                        onChange={(e) => setExportGeneralAlignments(e.target.checked)}
                        className="rounded bg-[#1f242e] border-[#2d3545] text-emerald-500"
                      />
                      <span className="text-[#c9d1d9]">
                        Export general alignments to <code>{generalAlignmentFolderName.trim() || '—'}/</code> using Catalog Pass/Fail
                      </span>
                    </label>
                    {recipe.enable_orf && (
                      <label className="flex items-center gap-2 cursor-pointer">
                        <input
                          type="checkbox"
                          checked={exportOrfAlignments}
                          onChange={(e) => setExportOrfAlignments(e.target.checked)}
                          className="rounded bg-[#1f242e] border-[#2d3545] text-violet-500"
                        />
                        <span className="text-[#c9d1d9]">
                          Export accepted ORF alignments independently to <code>{orfAlignmentFolderName.trim() || '—'}/</code>
                        </span>
                      </label>
                    )}
                    {recipe.enable_orf && recipe.orf_use_references && Object.keys(recipe.orf_reference_sequences ?? {}).length > 0 && (
                      <label className="flex items-center gap-2 cursor-pointer">
                        <input
                          type="checkbox"
                          checked={exportIntrons}
                          onChange={(e) => setExportIntrons(e.target.checked)}
                          className="rounded bg-[#1f242e] border-[#2d3545] text-blue-500"
                        />
                        <span className="text-[#c9d1d9]">
                          Export reference-trimmed introns independently to <code>{intronFolderName.trim() || '—'}/</code>
                        </span>
                      </label>
                    )}
                    <label className="flex items-center gap-2 cursor-pointer">
                      <input
                        type="checkbox"
                        checked={saveSummaryCsv}
                        onChange={(e) => setSaveSummaryCsv(e.target.checked)}
                        className="rounded bg-[#1f242e] border-[#2d3545] text-blue-500"
                      />
                      <span className="text-[#c9d1d9]">
                        Write summary CSV (<code>alignment-trimming_summary.csv</code>)
                      </span>
                    </label>

                    <label className="flex items-center gap-2 cursor-pointer">
                      <input
                        type="checkbox"
                        checked={saveRecipeJson}
                        onChange={(e) => setSaveRecipeJson(e.target.checked)}
                        className="rounded bg-[#1f242e] border-[#2d3545] text-blue-500"
                      />
                      <span className="text-[#c9d1d9]">
                        Save reproducible recipe (<code>recipe.json</code>)
                      </span>
                    </label>
                  </>
                ) : (
                  <>
                    <label className="flex items-center gap-2 cursor-pointer">
                      <input
                        type="checkbox"
                        checked={writeRaxmlPartitions}
                        onChange={(e) => setWriteRaxmlPartitions(e.target.checked)}
                        className={mode === 'group' ? "rounded bg-[#1f242e] border-[#2d3545] text-fuchsia-500" : "rounded bg-[#1f242e] border-[#2d3545] text-emerald-500"}
                      />
                      <span className="text-[#c9d1d9]">
                        Generate RAxML / RAxML-NG partition file (<code>_partitions.txt</code>)
                      </span>
                    </label>

                    <label className="flex items-center gap-2 cursor-pointer">
                      <input
                        type="checkbox"
                        checked={writeNexusPartitions}
                        onChange={(e) => setWriteNexusPartitions(e.target.checked)}
                        className={mode === 'group' ? "rounded bg-[#1f242e] border-[#2d3545] text-fuchsia-500" : "rounded bg-[#1f242e] border-[#2d3545] text-emerald-500"}
                      />
                      <span className="text-[#c9d1d9]">
                        Generate IQ-TREE / NEXUS partition block (<code>_partitions.nex</code>)
                      </span>
                    </label>
                  </>
                )}
              </div>

              {errorMsg && (
                <div className="p-2.5 rounded bg-rose-500/10 border border-rose-500/30 text-rose-300 text-xs">
                  {errorMsg}
                </div>
              )}

              {/* Action Buttons */}
              <div className="flex items-center justify-end gap-2 pt-3 border-t border-[#232833]">
                <button
                  onClick={onClose}
                  className="px-3.5 py-1.5 rounded-md bg-[#1b2029] hover:bg-[#232833] text-[#8b949e] hover:text-[#c9d1d9] text-xs font-medium transition-colors"
                >
                  Cancel
                </button>
                <button
                  onClick={handleExecuteExport}
                  disabled={isRunning}
                  className="px-4 py-1.5 rounded-md bg-blue-600 hover:bg-blue-500 text-white text-xs font-semibold shadow-md transition-colors flex items-center gap-1.5 disabled:opacity-50"
                >
                  {isRunning ? (
                    <span>Processing...</span>
                  ) : mode === 'batch' ? (
                    <>
                      <Download className="w-3.5 h-3.5" />
                      <span>Start Batch Export</span>
                    </>
                  ) : mode === 'concatenate' ? (
                    <>
                      <Layers className="w-3.5 h-3.5" />
                      <span>Assemble Supermatrix</span>
                    </>
                  ) : (
                    <>
                      <Sparkles className="w-3.5 h-3.5" />
                      <span>Batch Concatenate Exons</span>
                    </>
                  )}
                </button>
              </div>
            </>
          )}
        </div>
      </div>
    </div>
  );
};
