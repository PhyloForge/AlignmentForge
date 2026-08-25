use std::collections::BTreeMap;

/// Port of PhyloProcessR `collapseAlignmentSpecies.R`
/// Collapses per-individual sequences to species-level sequences by merging all taxa that share
/// the same Genus_species prefix (or specified delimiter) into a single combined consensus sequence.
/// Informative bases take priority over gaps/missing data.
pub fn collapse_species(
    taxa: &[String],
    sequences: &[String],
    species_delimiter: &str,
    max_missing_percent: f64,
) -> (Vec<String>, Vec<String>) {
    if sequences.is_empty() || taxa.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let length = sequences[0].len();

    // Group sequence indices by species prefix
    let mut species_groups: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (i, taxon) in taxa.iter().enumerate() {
        let species_key = if let Some(idx) = taxon.find(species_delimiter) {
            taxon[..idx].to_string()
        } else {
            // If no delimiter found, take first two underscore parts (e.g. Genus_species)
            let parts: Vec<&str> = taxon.split('_').collect();
            if parts.len() >= 2 {
                format!("{}_{}", parts[0], parts[1])
            } else {
                taxon.clone()
            }
        };
        species_groups.entry(species_key).or_default().push(i);
    }

    let mut collapsed_taxa = Vec::new();
    let mut collapsed_seqs = Vec::new();

    for (species, indices) in species_groups {
        let mut merged_seq = Vec::with_capacity(length);
        let mut gap_count = 0usize;

        for col in 0..length {
            let mut resolved_base = b'-';

            // Find first informative base among individuals
            for &idx in &indices {
                if let Some(&b) = sequences[idx].as_bytes().get(col) {
                    let u = b.to_ascii_uppercase();
                    if u != b'-' && u != b'?' && u != b'N' {
                        resolved_base = u;
                        break;
                    }
                }
            }

            if resolved_base == b'-' {
                gap_count += 1;
            }
            merged_seq.push(resolved_base);
        }

        let missing_pct = (gap_count as f64 / length as f64) * 100.0;
        if missing_pct <= max_missing_percent {
            collapsed_taxa.push(species);
            collapsed_seqs.push(String::from_utf8(merged_seq).unwrap_or_default());
        }
    }

    (collapsed_taxa, collapsed_seqs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_collapse_species() {
        let taxa = vec![
            "Rana_temporaria_Sample1".to_string(),
            "Rana_temporaria_Sample2".to_string(),
            "Bufo_bufo_Sample1".to_string(),
        ];
        let seqs = vec![
            "A-G-T".to_string(),
            "-T--T".to_string(),
            "CCCCC".to_string(),
        ];

        let (c_taxa, c_seqs) = collapse_species(&taxa, &seqs, "_Sample", 100.0);
        assert_eq!(c_taxa, vec!["Bufo_bufo".to_string(), "Rana_temporaria".to_string()]);
        assert_eq!(c_seqs[0], "CCCCC");
        assert_eq!(c_seqs[1], "ATG-T"); // Merged A from sample 1, T from sample 2, G from sample 1
    }
}
