pub mod algorithms;
pub mod commands;
pub mod export;
pub mod filter_config;
pub mod models;
pub mod parsers;
pub mod pipeline;
pub mod state;

use commands::*;
use state::AlignmentCache;

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        .manage(AlignmentCache::new())
        .plugin(tauri_plugin_opener::init())
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_fs::init())
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
