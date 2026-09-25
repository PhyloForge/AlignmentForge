use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use crate::models::Alignment;
use crate::pipeline::catalog::UnassessedSummary;

/// Per-locus taxon lists, without the sequences.
///
/// The dataset-wide occupancy filter only needs to know which samples appear in
/// which locus. Cloning every sequence to answer that question costs the whole
/// dataset in memory on each alignment the user opens.
pub type TaxonPresence = Vec<Vec<String>>;

/// Loci measured with one set of processing settings, before the thresholds.
pub struct MeasuredRun {
    /// Identifies the processing settings (`processing_key` in commands.rs).
    pub key: String,
    pub loci: Vec<UnassessedSummary>,
}

/// In-memory cache of parsed alignments, keyed by file path.
/// Stored in `tauri::State` so all commands can borrow from memory
/// instead of re-reading files from disk on every recipe change.
#[derive(Clone)]
pub struct AlignmentCache {
    inner: Arc<Mutex<HashMap<String, Arc<Alignment>>>>,
    /// The latest catalog run, so a change of pass thresholds only assesses
    /// these loci again.
    measured: Arc<Mutex<Option<Arc<MeasuredRun>>>>,
    catalog_generation: Arc<AtomicU64>,
}

/// Takes a lock, recovering from a poisoned mutex.
///
/// Poisoning only means some other thread panicked while holding the lock.
/// The data itself stays consistent, and losing the whole cache to that is
/// worse than continuing.
fn lock<T>(mutex: &Mutex<T>) -> std::sync::MutexGuard<'_, T> {
    mutex.lock().unwrap_or_else(|poisoned| poisoned.into_inner())
}

impl AlignmentCache {
    pub fn new() -> Self {
        Self {
            inner: Arc::new(Mutex::new(HashMap::new())),
            measured: Arc::new(Mutex::new(None)),
            catalog_generation: Arc::new(AtomicU64::new(0)),
        }
    }

    /// Marks a new catalog recalculation as current and invalidates older jobs.
    pub fn begin_catalog_job(&self) -> u64 {
        self.catalog_generation.fetch_add(1, Ordering::SeqCst) + 1
    }

    pub fn is_catalog_job_current(&self, generation: u64) -> bool {
        self.catalog_generation.load(Ordering::SeqCst) == generation
    }

    /// Replace the entire cache with a new set of alignments.
    pub fn store(&self, alignments: Vec<Alignment>) {
        {
            let mut map = lock(&self.inner);
            map.clear();
            map.reserve(alignments.len());
            for align in alignments {
                map.insert(align.file_path.clone(), Arc::new(align));
            }
        }
        *lock(&self.measured) = None;
    }

    /// Get a single alignment by file path.
    pub fn get(&self, file_path: &str) -> Option<Arc<Alignment>> {
        lock(&self.inner).get(file_path).cloned()
    }

    /// Taxon names per locus, without copying any sequence data.
    pub fn taxon_presence(&self) -> TaxonPresence {
        let map = lock(&self.inner);
        map.values().map(|align| align.taxa.clone()).collect()
    }

    pub fn unique_taxon_count(&self) -> usize {
        let map = lock(&self.inner);
        map.values()
            .flat_map(|alignment| alignment.taxa.iter())
            .collect::<HashSet<_>>()
            .len()
    }

    /// True when an alignment with this path is cached.
    pub fn contains(&self, file_path: &str) -> bool {
        lock(&self.inner).contains_key(file_path)
    }

    /// Get the alignments matching the given paths (preserving order).
    pub fn get_by_paths(&self, paths: &[String]) -> Vec<Arc<Alignment>> {
        let map = lock(&self.inner);
        paths.iter().filter_map(|p| map.get(p).cloned()).collect()
    }

    /// The latest measured run, when it used these settings on these paths.
    pub fn measured_run(&self, key: &str, paths: &[String]) -> Option<Arc<MeasuredRun>> {
        let measured = lock(&self.measured);
        let run = measured.as_ref()?;
        let same_loci = run.loci.len() == paths.len()
            && run
                .loci
                .iter()
                .zip(paths)
                .all(|(locus, path)| &locus.summary.file_path == path);
        (run.key == key && same_loci).then(|| run.clone())
    }

    /// The latest measured run, for its per-sample lists.
    pub fn latest_measured_run(&self) -> Option<Arc<MeasuredRun>> {
        lock(&self.measured).clone()
    }

    /// Keeps a measured run, unless a newer job has started since `generation`.
    /// The check and the store happen under one lock, so an older job cannot
    /// replace the run of a newer one.
    pub fn store_measured_run(
        &self,
        generation: u64,
        run: MeasuredRun,
    ) -> Option<Arc<MeasuredRun>> {
        let mut measured = lock(&self.measured);
        if !self.is_catalog_job_current(generation) {
            return None;
        }
        let run = Arc::new(run);
        *measured = Some(run.clone());
        Some(run)
    }
}
