pub mod batch;
pub mod concatenate;
pub mod group;

use std::fs;
use std::path::Path;

use crate::models::AlignmentFormat;
use crate::parsers::write_alignment;

/// Writes a file under a temporary name in the same folder, then renames it to
/// `path`. An export that stops part way never leaves a partial file under the
/// final name, and a failed write keeps the file of an earlier export.
pub fn write_atomic<F>(path: &Path, write: F) -> Result<(), String>
where
    F: FnOnce(&Path) -> Result<(), String>,
{
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| format!("Invalid output path: {}", path.display()))?;
    let temporary_path = path.with_file_name(format!(".{file_name}.tmp"));
    if let Err(error) = write(&temporary_path) {
        let _ = fs::remove_file(&temporary_path);
        return Err(error);
    }
    fs::rename(&temporary_path, path).map_err(|error| {
        let _ = fs::remove_file(&temporary_path);
        format!("Failed to finalize {}: {error}", path.display())
    })
}

pub fn write_alignment_atomic(
    path: &Path,
    taxa: &[String],
    sequences: &[String],
    format: AlignmentFormat,
) -> Result<(), String> {
    write_atomic(path, |temporary_path| {
        write_alignment(temporary_path, taxa, sequences, format)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn failed_write_leaves_no_partial_file() {
        let root = std::env::temp_dir().join("alignmentforge_atomic_write_test");
        let _ = fs::remove_dir_all(&root);
        fs::create_dir_all(&root).unwrap();
        let path = root.join("supermatrix.fa");
        fs::write(&path, "earlier export\n").unwrap();

        let result = write_atomic(&path, |temporary_path| {
            fs::write(temporary_path, ">a\nAC").unwrap();
            Err("disk full".to_string())
        });

        assert_eq!(result.unwrap_err(), "disk full");
        assert_eq!(fs::read_to_string(&path).unwrap(), "earlier export\n");
        assert_eq!(fs::read_dir(&root).unwrap().count(), 1, "the temporary file remains");
        let _ = fs::remove_dir_all(root);
    }
}
