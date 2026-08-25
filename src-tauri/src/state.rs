use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use crate::models::Alignment;

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

    /// Marks a new catalog recalculation as current and invalidates older jobs.
    pub fn begin_catalog_job(&self) -> u64 {
        self.catalog_generation.fetch_add(1, Ordering::SeqCst) + 1
    }

    pub fn is_catalog_job_current(&self, generation: u64) -> bool {
        self.catalog_generation.load(Ordering::SeqCst) == generation
    }

    /// Replace the entire cache with a new set of alignments.
    pub fn store(&self, alignments: Vec<Alignment>) {
        let mut map = self.inner.lock().unwrap();
        map.clear();
        map.reserve(alignments.len());
        for align in alignments {
            map.insert(align.file_path.clone(), align);
        }
    }

    /// Get a clone of a single alignment by file path.
    pub fn get(&self, file_path: &str) -> Option<Alignment> {
        let map = self.inner.lock().unwrap();
        map.get(file_path).cloned()
    }

    /// Get clones of all cached alignments.
    pub fn get_all(&self) -> Vec<Alignment> {
        let map = self.inner.lock().unwrap();
        map.values().cloned().collect()
    }

    /// Get clones of alignments matching the given paths (preserving order).
    pub fn get_by_paths(&self, paths: &[String]) -> Vec<Alignment> {
        let map = self.inner.lock().unwrap();
        paths
            .iter()
            .filter_map(|p| map.get(p).cloned())
            .collect()
    }
}
