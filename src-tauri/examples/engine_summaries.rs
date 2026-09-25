//! Prints the desktop engine's catalog summaries and viewer data as JSON, for
//! the browser parity check in `scripts/check-engine-parity.mjs`.
//!
//! Reads `{ "files": [...], "recipes": [...] }` from standard input. For each
//! recipe, prints one list of summaries and one viewer entry for each file.

use std::collections::HashSet;
use std::io::Read;

use alignment_forge_lib::algorithms::informative::calculate_parsimony_informative_sites;
use alignment_forge_lib::algorithms::stats::compute_majority_consensus;
use alignment_forge_lib::parsers::parse_alignment;
use alignment_forge_lib::pipeline::catalog::{
    evaluate_recipe_on_alignments_with_progress, recipe_with_taxon_presence,
};
use alignment_forge_lib::pipeline::engine::apply_recipe;
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

    let summaries: Vec<_> = request
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

    // The viewer data, as the `get_alignment` command builds it.
    let dataset_taxa: Vec<&[String]> = alignments
        .iter()
        .map(|alignment| alignment.taxa.as_slice())
        .collect();
    let views: Vec<Vec<serde_json::Value>> = request
        .recipes
        .iter()
        .map(|recipe| {
            let runtime_recipe = recipe_with_taxon_presence(recipe, &dataset_taxa);
            alignments
                .iter()
                .map(|alignment| {
                    let (trimmed_alignment, diff) =
                        apply_recipe(alignment, &runtime_recipe, total_unique_taxa);
                    let (_, _, pis_mask) =
                        calculate_parsimony_informative_sites(&alignment.sequences, true);
                    serde_json::json!({
                        "trimmed_alignment": trimmed_alignment,
                        "diff": diff,
                        "pis_mask": pis_mask,
                        "majority_consensus": compute_majority_consensus(&alignment.sequences, false),
                    })
                })
                .collect()
        })
        .collect();

    let output =
        serde_json::to_string(&serde_json::json!({ "summaries": summaries, "views": views }))
            .map_err(|error| format!("Failed to serialize the summaries: {error}"))?;
    println!("{output}");
    Ok(())
}
