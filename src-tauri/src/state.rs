use std::collections::HashMap;
use std::collections::HashSet;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use crate::models::Alignment;

/// Per-locus taxon lists, without the sequences.
///
/// The dataset-wide occupancy filter only needs to know which samples appear in
/// which locus. Cloning every sequence to answer that question costs the whole
/// dataset in memory on each alignment the user opens.
pub type TaxonPresence = Vec<Vec<String>>;

/// In-memory cache of parsed alignments, keyed by file path.
/// Stored in `tauri::State` so all commands can borrow from memory
/// instead of re-reading files from disk on every recipe change.
#[derive(Clone)]
pub struct AlignmentCache {
    inner: Arc<Mutex<HashMap<String, Alignment>>>,
    catalog_generation: Arc<AtomicU64>,
}

impl AlignmentCache {
    pub fn new() -> Self {
        Self {
            inner: Arc::new(Mutex::new(HashMap::new())),
            catalog_generation: Arc::new(AtomicU64::new(0)),
        }
    }

    /// Takes the cache lock, recovering from a poisoned mutex.
    ///
    /// Poisoning only means some other thread panicked while holding the lock.
    /// The map itself stays consistent, and losing the whole cache to that is
    /// worse than continuing.
    fn lock(&self) -> std::sync::MutexGuard<'_, HashMap<String, Alignment>> {
        self.inner.lock().unwrap_or_else(|poisoned| poisoned.into_inner())
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
        let mut map = self.lock();
        map.clear();
        map.reserve(alignments.len());
        for align in alignments {
            map.insert(align.file_path.clone(), align);
        }
    }

    /// Get a clone of a single alignment by file path.
    pub fn get(&self, file_path: &str) -> Option<Alignment> {
        let map = self.lock();
        map.get(file_path).cloned()
    }

    /// Get clones of all cached alignments.
    ///
    /// This copies every sequence in the dataset. Prefer `taxon_presence` when
    /// only the sample names are needed.
    pub fn get_all(&self) -> Vec<Alignment> {
        let map = self.lock();
        map.values().cloned().collect()
    }

    /// Taxon names per locus, without copying any sequence data.
    pub fn taxon_presence(&self) -> TaxonPresence {
        let map = self.lock();
        map.values().map(|align| align.taxa.clone()).collect()
    }

    pub fn unique_taxon_count(&self) -> usize {
        let map = self.lock();
        map.values()
            .flat_map(|alignment| alignment.taxa.iter())
            .collect::<HashSet<_>>()
            .len()
    }

    /// Number of cached alignments.
    pub fn len(&self) -> usize {
        self.lock().len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// True when an alignment with this path is cached.
    pub fn contains(&self, file_path: &str) -> bool {
        self.lock().contains_key(file_path)
    }

    /// Get clones of alignments matching the given paths (preserving order).
    pub fn get_by_paths(&self, paths: &[String]) -> Vec<Alignment> {
        let map = self.lock();
        paths
            .iter()
            .filter_map(|p| map.get(p).cloned())
            .collect()
    }
}
