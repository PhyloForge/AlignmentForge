# AlignmentForge code review and implementation plan

Review date: 8 September 2026. This review covers the current working files, including the existing uncommitted changes. It focuses on the Rust filtering pipeline, parsers, export paths, browser workers, and the interface that controls them.

Implement the numbered items in order, except where a dependency is stated. P1 means that the issue can change sequence data, retain the wrong samples, or produce an incorrect export. P2 means a material correctness, reliability, or performance improvement. Optional items need not block the corrections.

Keep the existing Rust modules and React components. Use ordinary loops and small functions. Do not change biological assumptions as part of a performance edit. Do not commit, push, or build an installer.

## Checks completed

- The existing Rust library suite passes: 70 tests, 0 failures.
- TypeScript checking passes with `tsc --noEmit`.
- The display helper check passes.
- ESLint reports 0 errors and 23 warnings. Most warnings concern unused names; two concern effect dependencies. These warnings are not the main review findings.
- A temporary Rust program imported the current engine, parser, and export modules and ran synthetic cases. It confirmed the PHYLIP, NEXUS, coverage, consensus, statistics, annotation, empty-output, ORF-frame, and export issues described below. It did not read user sequence data into the conversation.
- Small readability edits were made in `parsers/mod.rs`, `pipeline/recipe.rs`, and the `extend_state` signature in `algorithms/reference.rs`. These edits expand compressed statements. They do not change filtering or parsing behavior.
- No handoff file is tracked. The global handoff ignore rule applies.

The review did not run a packaged build, a browser interaction test, a scientific comparison against external trimming programs, or a performance benchmark. Speed improvements below are supported by code inspection; no speedup factor is claimed. Passing the current tests does not establish that all filter combinations are correct.

## Numbered recommendations

1. **P1 — Correct unnamed continuation rows in interleaved PHYLIP input.**

   Location: [parsers/phylip.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/parsers/phylip.rs:109), especially selection of the first sequence shorter than `expected_length`.

   **Evidence:** A two-taxon, three-block file with declared length 12 parses successfully but assigns sequence chunks to the wrong taxa. The first block is `a AAAA`, `b CCCC`; the second block is `GGGG`, `TTTT`; the third block is `ACGT`, `TGCA`. Expected sequences are `AAAAGGGGACGT` and `CCCCTTTTTGCA`. Both parsed sequences differ from these expectations. Shape checks cannot detect this because both wrong sequences have length 12.

   **Implementation:** Preserve block boundaries long enough to distinguish sequential records from interleaved records. Once a first interleaved block defines all taxa, advance an explicit row index through every continuation block. Do not select the first incomplete sequence for that layout. Keep the existing named-continuation path. Use separate, short sequential and interleaved loops if needed; do not add a general parser framework. Reject incomplete blocks with a useful error.

   **Validation:** Add the exact three-block case above. Also test two blocks, four blocks, named later blocks, sequential wrapped sequences, strict ten-column names, and a truncated last block. Assert taxon-to-sequence identity, not only lengths. This correction changes parsed sequence data for affected files.

2. **P1 — Read NEXUS matrix metadata and taxon names correctly.**

   Location: [parsers/nexus.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/parsers/nexus.rs:24).

   **Evidence:** `FORMAT DATATYPE=DNA MATCHCHAR=.;` with rows `a ATGC` and `b ....` leaves four literal dots in the second sequence. Those dots then enter filter calculations as observed characters. Changing `NCHAR=4` to `NCHAR=40` is also accepted. The parser currently supplies no declared width to `validate_alignment_shape`. It splits quoted names on whitespace, and it does not handle a matrix terminator attached to the last sequence token.

   **Implementation:** Read `NTAX`, `NCHAR`, `GAP`, `MISSING`, and `MATCHCHAR` from the supported DATA/CHARACTERS block. Normalize declared gap and missing symbols to the engine's existing symbols. Resolve match characters against the first sequence after assembling interleaved chunks. Check declared dimensions. Use a small tokenizer for quoted taxon names, comments, and the matrix semicolon. Reject unsupported syntax explicitly if it cannot be parsed safely. Quote and escape names in `write_nexus` so exported names with spaces round-trip.

   **Validation:** Test match characters, custom missing/gap symbols, quoted names with spaces and escaped quotes, inline comments, an attached semicolon, declared taxon/width mismatches, and a missing matrix terminator. Confirm that equivalent FASTA and NEXUS inputs produce identical normalized matrices and filter outcomes. Do not silently interpret protein data as nucleotide data; see item 8.

3. **P1 — Pass the same dataset context to catalog, preview, and every export.**

   Locations: [commands.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/commands.rs:234), [export/batch.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/export/batch.rs:138), [export/concatenate.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/export/concatenate.rs:53), and [export/group.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/export/group.rs:100).

   **Evidence:** Exports call `apply_recipe(..., 0)`. The engine then uses each input locus's taxon count as the occupancy denominator. A four-taxon locus from an eight-taxon dataset fails a 75% catalog occupancy threshold but is exported by batch export with `only_passing=true`. A second problem is the use of `excluded_taxa.is_empty()` to decide whether dataset exclusions have been computed. An empty list is a valid computed result. Selected-locus and grouped exports can therefore recompute occupancy on a smaller subset and remove samples that the catalog retained.

   **Implementation:** Resolve the original dataset's unique-taxon count and runtime recipe once in the command layer. Pass that count explicitly to each export function and each `apply_recipe` call. Treat the supplied runtime recipe as resolved even when its exclusion list is empty. For direct export callers without a loaded dataset, resolve context once from all requested inputs before selecting gene groups. Do not recompute it inside a group. Keep the documented independence of general-alignment assessment and ORF acceptance.

   **Validation:** Use two loci with a total union of eight taxa, one containing four taxa. At 75% locus occupancy, catalog, preview, general export, and concatenation must agree. Add a selected-export case where a taxon passes dataset sample occupancy but fails occupancy within the selected subset. Test both an empty exclusion list and a nonempty list. Check grouped export against the full requested dataset, not one gene at a time.

4. **P1 — Apply absolute and percentage coverage thresholds without bypasses.**

   Locations: [algorithms/coverage.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/algorithms/coverage.rs:14) and [pipeline/engine.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/pipeline/engine.rs:223).

   **Evidence:** For `AAAA`, `AA--`, and `----`, a minimum of 4 bp retains all three rows; a minimum of 5 bp also retains all three. The early return when `min_coverage_bp >= alignment_len` bypasses both thresholds. A two-row alignment with `AAAA` and `----` also retains both rows with a 3 bp minimum because the engine and helper both skip coverage filtering at two taxa.

   **Implementation:** Keep the empty-input guard. Remove the width-based bypass and apply the existing `count >= min_bp && percent >= min_percent` rule for all nonempty taxon sets. Remove the corresponding engine condition that requires more than two rows. Coverage is an individual-row measurement and does not need a consensus. Preserve `RelativeWidth::Sample` as the longest observed sample before this filtering pass. Do not remove low-sample guards from consensus-based filters as an unrelated change.

   **Validation:** With the three rows above, 4 bp must retain only `AAAA`; 5 bp must retain none. Test one and two rows, zero observed bases, lowercase `n`, `?`, threshold equality, percentage-only filtering, and both relative-width options. Add a pipeline case where manual exclusions reduce a locus to two rows before coverage runs.

5. **P1 — Use an exact consensus for filters that remove samples or mask sequence.**

   Location: [algorithms/stats.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/algorithms/stats.rs:56); callers include `similarity.rs` and `segments.rs`.

   **Evidence:** The consensus samples rows at a fixed stride. In a 100-row fixture with 60 `AAAA` rows and 40 `CCCC` rows, placing the C rows at the first 40 even indices makes the filter retain all 40 C rows and remove all 60 A rows. Reordering the same sequences makes it retain all 60 A rows and remove all 40 C rows. No biological data changed.

   **Implementation:** Count every row for the consensus used by filtering. Change the count array from `u16` to `usize` when removing sampling. Preserve the existing character treatment and deterministic tie order in this correction. If approximate divergence is retained for a display metric, give that approximation an explicit name and description; do not use it to make removal decisions.

   **Validation:** Add the 100-row fixture and verify the retained taxon set after several deterministic permutations. Also test 49, 50, 99, and 100 rows, majority ties, and missing-only columns. Apply the same invariance test to segment masking. Measure runtime after the correctness fix; do not trade row-order independence for speed.

6. **P1 — Calculate final statistics after MACSE normalization.**

   Locations: [pipeline/engine.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/pipeline/engine.rs:674), [algorithms/informative.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/algorithms/informative.rs:106), and `pipeline/catalog.rs`.

   **Evidence:** With `AA!`, `AA!`, `AAA`, and `AAA`, other sequence filters disabled, the diff reports 0% gaps and one informative site. The returned sequences contain `N` in place of `!`, so their actual values are 16.7% gaps and zero informative sites. The locus passes a 10% maximum gap threshold and a minimum of one informative site despite failing both in the output.

   **Implementation:** Keep `!` available until MACSE QC has inspected it. Convert surviving markers to `N` before calculating final gap/site statistics and applying assessment. Make catalog metrics use the same final matrix. In nucleotide site statistics, treat `!` as missing rather than as a resolved biological state. Define the same treatment for raw statistics without altering the stored raw sequences. Apply that definition consistently to observed-base totals and coverage where these are intended to count biological observations.

   **Validation:** Add the four-row case above. Compare every reported final statistic with a direct calculation from `transformed_alignment.sequences`. Test ORF enabled and disabled, retained internal markers, terminal marker masking, and ordinary sequences without markers. Keep the existing MACSE rejection thresholds unchanged.

7. **P1 — Report export failures and prevent output-name collisions.**

   Locations: [export/batch.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/export/batch.rs:132), [export/concatenate.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/export/concatenate.rs:47), [export/group.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/export/group.rs:89), and `components/ExportModal.tsx`.

   **Evidence:** Batch export returns `Ok` when its only input path does not exist. It also returns `Ok` when the intended output file path is an existing directory. Parsing uses `.ok()`, alignment writing uses `.is_ok()`, and several sidecar writes discard errors. Group export can increment its success count after a failed write. Batch outputs are named by file stem, so two input directories containing the same stem can write the same destination in parallel.

   **Implementation:** Precompute destination paths and detect duplicate destinations before writing. Include case-insensitive collisions where relevant to the output filesystem. Use a simple per-file error list for batch export, with input path, output path, and error text. Separate failed, filtered, and successfully exported counts. For concatenation, stop with an error if an input cannot be read; do not silently create a smaller matrix. Propagate summary, recipe, and partition write failures. Return a saved path only after the write succeeds. Write alignment files to sibling temporary files and rename them after successful completion so an interrupted write is not presented as a complete alignment.

   **Implementation detail:** Use the existing `csv` dependency to write the summary CSV, rather than comma interpolation. Gene names used as output filenames must be single filename components; reject separators, `.` and `..`. Use that same existing CSV reader for quoted mapping fields. These are ordinary filename and table cases, not a new validation framework.

   **Validation:** Test a missing input, an output path that is a directory, duplicate stems from separate input directories, a failed sidecar write, a locus name containing a comma, and a gene name containing a path separator. Confirm that failure counts and UI messages identify the affected file and that incomplete files are not counted as exported.

8. **P1 — Reject duplicate sample identifiers and basic malformed sequence records.**

   Locations: [parsers/mod.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/parsers/mod.rs:13), [parsers/fasta.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/parsers/fasta.rs:39), and the other parser entry points.

   **Evidence:** Two FASTA records named `same` with sequences `AAAA` and `CCCC` are accepted. Later maps keyed by taxon name can keep only one sequence, while row counts still include both. FASTA also collects sequence bytes before the first header. Valid UTF-8 is accepted as sequence text even though downstream algorithms use ASCII byte coordinates, and some masking code indexes character vectors with those byte positions.

   **Implementation:** Extend the existing shared validation function with nonempty, unique taxon-name checks and a basic sequence-alphabet check. Preserve case. Permit the nucleotide ambiguity and MACSE symbols the application intentionally supports. Reject sequence content before the first FASTA header, empty names, and non-ASCII sequence characters with the record name and position. Keep Unicode in taxon names where writers can preserve it. Do not silently merge duplicate samples or change their names.

   **Scope decision:** The scanner accepts `.faa`, but the filters and NEXUS writer assume nucleotides. Make the supported input alphabet explicit. The small correction is to reject unsupported protein input with a clear message; adding full protein filtering is a separate feature. Also handle taxon names with whitespace safely when exporting PHYLIP: reject unsupported names or use a documented mapping, rather than emitting an ambiguous file.

   **Validation:** Test duplicate and empty names, sequence text before a header, a non-ASCII sequence character, empty records, valid lowercase/IUPAC/MACSE input, and names with spaces. A duplicate record must fail before statistics, occupancy, or concatenation runs.

9. **P2 — Convert every mask annotation to raw alignment coordinates.**

   Locations: [pipeline/engine.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/pipeline/engine.rs:365), the corresponding segment/fallback/ORF append operations, and [MsaViewer.tsx](/Users/chutter/Bioinformatics/Github/AlignmentForge/src/components/MsaViewer.tsx:685).

   **Evidence:** For two `--AAAAAAAA` rows and one `--AAAACCCC` row, removing the two initial gap columns and masking with a four-base window produces a mask at `[4, 8)`. Its actual raw position is `[6, 10)`. The viewer compares the reported positions directly with raw columns. Reference slicing, reverse complementation, and ORF boundary trimming introduce additional offsets.

   **Implementation:** Define `TrimmingDiff.masked_segments` as zero-based, half-open raw intervals. Add one small helper that maps a stage-local interval through that stage's column map and emits contiguous raw runs. Split intervals across removed columns rather than shading the intervening raw positions. Map HMM and window masks immediately when generated. For ORF masks, account for both orientation and the columns removed by ORF trimming before translating the returned positions. Preserve taxon names and the existing raw-column reason map.

   **Validation:** Add the exact prefix-deletion fixture. Also test an internal deleted column, a reference-anchored region, a negative reading frame, and a masked codon after a frame offset. Assert exact raw intervals numerically. No screenshot inspection is required. This corrects annotations; it must not change the sequence output.

10. **P2 — Make empty-output and threshold-endpoint behavior explicit.**

   Locations: [pipeline/engine.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/pipeline/engine.rs:687), `pipeline/catalog.rs::apply_assessment`, and `algorithms/assess.rs`.

   **Evidence:** Three `----` rows become zero-width sequences after gap-only cleanup, but `diff.pass` is true when optional assessment is disabled. Export subsequently skips the locus because its length is zero. Also, a valid configuration with `max_gap_percent=0` disables the gap gate: a 50%-gap locus passes. The configuration currently accepts zero without explaining that sentinel behavior.

   **Implementation:** Treat zero taxa and zero output width as unconditional pipeline failures, with distinct reasons, in both direct and cached assessment paths. For the gap maximum, explicitly choose and document the zero-value convention. Recommended behavior is that a maximum of 0 means no gaps are permitted, while 100 permits any gap percentage; update the comparison accordingly if that convention is adopted. If zero must remain a legacy disable value, describe it in configuration comments and the UI. Do not silently reinterpret existing saved recipes without documenting the choice.

   **Validation:** Test zero rows, zero columns with rows retained, all columns removed by different filters, assessment on/off, and gap percentages just below, equal to, and above the configured limit. Check agreement among summary status, preview status, and export eligibility.

11. **P2 — Keep ORF candidate coordinates, output frame, and QC frame consistent.**

   Location: [algorithms/orf.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/algorithms/orf.rs:930).

   **Evidence:** With three `AATGGCTGCTTAA` rows, ORF enabled, automatic frame shifting disabled, and stop action `Keep`, the result selects frame +2, keeps the original 13 bases, and reports `found_valid_orf=true`. Codon QC still starts at column zero. Candidate measurements and QC therefore refer to different triplets.

   **Implementation:** Preserve candidate discovery and its reported frame. Before QC, decide whether the selected interval can be produced under the enabled transformations. If the candidate needs a forbidden shift, preserve the unmodified sequence and report that the candidate cannot be applied; do not validate it using another frame. Define the treatment of incomplete terminal codons when shifting is disabled. Make final stop detection, ORF validity, and ORF export acceptance use the same applied-frame outcome.

   **Validation:** Test the +2 fixture above, +3, a negative frame with reverse flipping disabled, and a complete frame +1 alignment. Test each stop action. Retain the deliberate independence of Catalog QC and ORF acceptance; an unavailable ORF need not invalidate the general-alignment branch.

12. **P2 — Use the same masking and cleanup steps in the reference fallback path.**

   Location: [pipeline/engine.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/pipeline/engine.rs:458).

   **Code inspection:** The main path runs HMM/window masking followed by missing-only column cleanup. The hybrid fallback restores the pre-reference matrix and reruns masking, but omits that cleanup before ORF search. Candidate coordinates and frame discovery can therefore depend on whether a failed reference attempt occurred. Existing fallback tests do not enable masking.

   **Implementation:** Extract only the duplicated masking-plus-cleanup sequence into a small local helper, or explicitly add the missing cleanup with the same column-map updates. Reuse it in both paths. Preserve the fallback restore point and discard annotations from the abandoned attempt.

   **Scientific decision:** Before changing cleanup for coding data, specify whether a column of `N`/`?` represents an unknown biological position or a removable alignment column. Removing one or two positions can alter downstream triplets. A codon-preserving mode needs an explicit rule, such as retaining partial missing triplets. Record this decision and test it; do not change the general missing-column rule solely to make a test pass.

   **Validation:** Extend the existing guided-failure fixture with masking that makes an internal region missing-only. Compare fallback processing with candidate processing from the same restored matrix. Check retained raw coordinates, sequence width, candidate frame, and mask annotations. Add coding cases with one, two, and three missing positions.

13. **P2 — Describe statistical and profile filters according to their implemented behavior.**

   Locations: [algorithms/statistical_columns.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/algorithms/statistical_columns.rs:154), `algorithms/hmm.rs`, `components/FilterSidebar.tsx`, and filter configuration comments.

   **Code inspection:** The `Gappyout (Gap Elbow)` UI option uses an average-similarity multiplier and fixed gap cutoffs; there is no elbow calculation in that branch. The block filter accepts `_max_nonconserved` but does not use it. The HMM-labelled function calculates leave-one-out emission ratios and smooths them; it has no hidden-state transitions or forward/backward computation. These observations do not establish equivalence to external programs with similar names.

   **Implementation:** Keep current numerical behavior while labeling these as AlignmentForge similarity, conserved-block, and profile-confidence heuristics. State the actual rule in help text and method metadata. Preserve serialized enum values for compatibility if labels change. Document `stat_col_max_nonconserved` as unsupported/deprecated rather than implying it is active. If bridging nonconserved runs or exact external-tool compatibility is wanted, implement it as a separate, specified feature after confirming scientific requirements.

   **Validation:** Add small fixed matrices showing each heuristic's threshold and gap behavior. Add a conserved-block case with an internal nonconserved run to document current behavior. Check that the UI, saved recipe comments, and export report name the same method. Do not claim external-tool parity without reference-output tests.

14. **P2 — Settle browser requests on cancellation and expose calculation failures.**

   Locations: [enginePool.ts](/Users/chutter/Bioinformatics/Github/AlignmentForge/src/engine/enginePool.ts:48), [engineWorker.ts](/Users/chutter/Bioinformatics/Github/AlignmentForge/src/engine/engineWorker.ts:89), and [App.tsx](/Users/chutter/Bioinformatics/Github/AlignmentForge/src/App.tsx:427).

   **Code inspection:** Worker termination clears pending requests without rejecting their promises. A worker error rejects existing requests but leaves the failed worker available for later requests. Catalog recalculation errors are written to the console; the displayed summaries can remain from another recipe. The chunk loop posts progress but does not yield to incoming worker messages, so obsolete runs can queue up during filter edits. The UI generation check prevents stale result assignment, but it does not cancel the computation.

   **Implementation:** Reject all pending promises with a clear cancellation error before termination. Mark failed workers unavailable and require a clean reload or explicit recovery. Show a catalog calculation error near the controls and distinguish the last successfully applied recipe from pending/failed settings. Coalesce queued recalculations to the newest recipe. If adding yields between chunks, add a generation check before each chunk and before finish, because the WASM session has one mutable pending buffer. Keep progress callbacks associated with request IDs.

   **Validation:** With a small fake Worker, test termination during a request, a worker error followed by a new request, rapid recipe A/B/C edits, and a dataset switch. Every promise must settle. Only the current generation may update summaries or progress. Add one real worker smoke test when the browser harness is available; no image capture is needed.

15. **P2 — Remove avoidable matrix copies before changing algorithms for speed.**

   Locations: [pipeline/engine.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/pipeline/engine.rs:281), [pipeline/catalog.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/pipeline/catalog.rs:548), and [wasm_api.rs](/Users/chutter/Bioinformatics/Github/AlignmentForge/src-tauri/src/wasm_api.rs:306).

   **Code inspection:** The engine builds raw ORF overlay sequences and fallback snapshots even when ORF processing is off. Catalog summarization clones the complete trimmed matrix and diff when the ORF branch does not differ. Browser assessment clones large taxon-retention maps and immediately empties them before serialization. Recipes can also carry a whole reference bank, which is cloned repeatedly during per-locus summarization.

   **Implementation:** First allocate ORF overlay data only when ORF analysis requires it, and allocate fallback snapshots only for the hybrid reference mode. Borrow the existing alignment/diff when the two summary branches are identical. Provide a lightweight summary path that does not construct then discard retention fields. Measure reference-bank copying separately; if material, borrow the matching reference sequence during per-locus processing instead of copying the full bank. Keep this change within existing engine functions and simple data types.

   **Validation:** Compare full outputs, masks, reasons, and summaries before and after each optimization on fixed synthetic nucleotide datasets. Measure release runtime and peak memory for ORF off, ORF on without references, and reference-guided processing. Launch each benchmark with output redirected to a file and inspect it after completion. Report dataset dimensions and build mode. Keep only changes with a measurable benefit; no target speedup is assumed.

16. **P2 — Reuse measured summaries for desktop assessment-only edits.**

   Locations: `commands.rs::recalculate_catalog`, `state.rs::AlignmentCache`, `pipeline/catalog.rs::apply_assessment`, and `App.tsx` assessment dependencies.

   **Code inspection:** The browser caches unassessed measurements and reapplies pass thresholds. The desktop recalculation path clones cached alignments and reruns sequence processing even when only a taxon, length, gap, or informative-site gate changes. Reference-guided processing makes that repeated work more expensive.

   **Implementation:** Keep one cached desktop run of unassessed summaries with its sequence-processing recipe identity and dataset generation. For assessment-only changes, call the existing `apply_assessment` function on that run. Invalidate on dataset replacement, manual exclusions, reference changes, or a processing-field change. Include reference contents or an explicit reference revision: `orf_reference_sequences` is skipped by normal recipe serialization, so a plain serialized recipe is not a sufficient cache key. Start with one cache entry; do not add a generic cache framework.

   **Validation:** Verify identical summaries from fresh processing and cached assessment at threshold boundaries. Instrument the test so a gate-only edit demonstrably skips sequence processing. Check invalidation after replacing a reference, changing a filter, and loading a new dataset. Benchmark this separately from item 15.

17. **Optional — Export an analysis record that can reproduce the result.**

   Locations: `export/batch.rs`, `pipeline/recipe.rs::orf_reference_sequences`, and `components/ExportModal.tsx`.

   **Evidence and benefit:** Batch export labels `recipe.json` as reproducible, but reference sequences are omitted by serialization. A reference-guided recipe alone cannot reproduce its output. Dataset scope, source changes after loading, and the independent Catalog/ORF outcomes also need to be clear to a later reader.

   **Implementation:** Keep ordinary filter configuration small. For an export run, optionally save the reference sequences as a separate FASTA file, the effective recipe, input paths and identifying file metadata, engine version, original dataset taxon count, exclusion lists, and per-locus outcome/error reasons. Record whether export rereads source files or uses the loaded snapshot; if sources have changed since loading, require a reload or report the mismatch. Add a compact per-stage count report for removed taxa, removed columns, and masked bases if it can be gathered from existing stage results. Use plain JSON/CSV files and the existing writers.

   **Validation:** Reload a saved run with its references and compare output matrices and statuses. Include a reference-guided locus, a failed locus, and a manually discarded sample. Confirm that a missing reference sidecar produces a visible message. This feature must not change filter rules.

18. **P2 — Turn the confirmed cases into permanent regression coverage.**

   Locations: existing `#[cfg(test)]` modules, `scripts/check-display-helpers.ts`, `package.json`, and `.github/workflows/ci.yml`.

   **Implementation:** Add each regression beside the code it checks as its correction is implemented. Use short in-memory matrices and temporary files generated by the test. Replace core export tests that run assertions only when `../test_data/...` exists with self-contained fixtures. Keep example-data tests separate. Add a thin native/WASM comparison for parser and recipe entry points; sharing Rust source does not test serialization, dataset context, or worker behavior. The comment in `defaultRecipe.ts` refers to `npm run parity`, but that script is absent: implement a real preset/entry-point parity check or correct the comment.

   **Required comparisons:** Exact retained taxa and sequences; raw/final coordinate maps; gap, variable-site, and informative-site statistics; Catalog and ORF status in their intended separate branches; exported alignment contents; partition lengths and one-based boundaries; surfaced file errors. For row-order tests, compare taxon-to-sequence mappings rather than output row order.

   **Validation commands:** Run `cargo test --manifest-path src-tauri/Cargo.toml --lib`, `npm run lint`, `npm run check:display`, and `tsc --noEmit` after the relevant fixes. Build WASM and run the new parity check when changing the WASM boundary or shared recipe contract. Redirect lengthy output to logs. Do not build a DMG. Avoid repository-wide formatting and unrelated warning cleanup.

## Implementation boundaries and completion criteria

The current plan contains proposed behavior changes; those changes have not been applied during this review. Only the small readability edits listed above are complete.

Items 1–8 address the highest-impact data and export problems. Items 9–14 complete coordinate, ORF, and runtime reliability work. Items 15–16 are performance work and must preserve validated outputs. Item 17 is an optional feature. Item 18 is implemented alongside the corrections, not postponed until the end.

Before calling the work complete, the confirmed synthetic failures must have permanent tests with the expected results stated here. Export must agree with the chosen dataset context and must report partial failure. Document the scientific decisions in items 10–13 before changing their semantics. Do not claim that an ORF is biologically correct solely because it has no internal stop codons, or that the local heuristic filters reproduce an external program without a comparison.
