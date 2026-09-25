pub mod algorithms;
pub mod models;
pub mod parsers;
pub mod pipeline;

// Desktop-only surface. The browser build compiles the same engine to
// WebAssembly, where there is no filesystem, no Tauri shell, and no thread
// pool, so these modules are left out of that target.
#[cfg(not(target_arch = "wasm32"))]
pub mod commands;
#[cfg(not(target_arch = "wasm32"))]
pub mod export;
#[cfg(not(target_arch = "wasm32"))]
pub mod filter_config;
#[cfg(not(target_arch = "wasm32"))]
pub mod state;

#[cfg(target_arch = "wasm32")]
pub mod wasm_api;

#[cfg(not(target_arch = "wasm32"))]
use commands::*;
#[cfg(not(target_arch = "wasm32"))]
use state::AlignmentCache;

#[cfg(not(target_arch = "wasm32"))]
#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .manage(AlignmentCache::new())
        .plugin(tauri_plugin_dialog::init())
        .invoke_handler(tauri::generate_handler![
            scan_directory,
            get_alignment,
            recalculate_catalog,
            run_batch_export,
            run_concatenate,
            run_grouped_concatenate,
            save_filter_config,
            load_filter_config,
            save_alignment_stats_csv,
            get_presets,
        ])
        .run(tauri::generate_context!())
        .expect("error while running AlignmentForge application");
}
