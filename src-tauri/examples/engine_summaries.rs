//! Prints the desktop engine's catalog summaries as JSON, for the browser
//! parity check in `scripts/check-engine-parity.mjs`.
//!
//! Reads `{ "files": [...], "recipes": [...] }` from standard input and prints
//! one list of summaries for each recipe.

use std::collections::HashSet;
use std::io::Read;

use alignment_forge_lib::parsers::parse_alignment;
use alignment_forge_lib::pipeline::catalog::evaluate_recipe_on_alignments_with_progress;
use alignment_forge_lib::pipeline::recipe::TrimmingRecipe;

#[derive(serde::Deserialize)]
struct Request {
    files: Vec<String>,
    recipes: Vec<TrimmingRecipe>,
}

fn main() -> Result<(), String> {
    let mut input = String::new();
    std::io::stdin()
        .read_to_string(&mut input)
        .map_err(|error| format!("Failed to read the request: {error}"))?;
    let request: Request =
        serde_json::from_str(&input).map_err(|error| format!("Invalid request: {error}"))?;

    let alignments = request
        .files
        .iter()
        .map(parse_alignment)
        .collect::<Result<Vec<_>, _>>()?;
    // The desktop app passes the number of unique samples from the folder scan.
    let total_unique_taxa = alignments
        .iter()
        .flat_map(|alignment| alignment.taxa.iter())
        .collect::<HashSet<_>>()
        .len();

    let runs: Vec<_> = request
        .recipes
        .iter()
        .map(|recipe| {
            evaluate_recipe_on_alignments_with_progress(
                &alignments,
                recipe,
                total_unique_taxa,
                None::<fn(usize, usize, &str)>,
            )
            .0
        })
        .collect();
    let output = serde_json::to_string(&runs)
        .map_err(|error| format!("Failed to serialize the summaries: {error}"))?;
    println!("{output}");
    Ok(())
}
