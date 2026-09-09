fn main() {
    // A build script runs on the host, so the target has to come from Cargo's
    // environment rather than a cfg attribute. The browser target compiles the
    // engine only; Tauri's build step generates desktop shell artifacts and
    // must not run for WebAssembly.
    let target = std::env::var("TARGET").unwrap_or_default();
    if target.contains("wasm32") {
        return;
    }
    tauri_build::build();
}
