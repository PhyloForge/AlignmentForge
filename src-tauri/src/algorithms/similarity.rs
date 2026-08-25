use crate::algorithms::stats::{compute_majority_consensus, pairwise_distance_to_target};

/// Port of PhyloProcessR `trimSampleSimilarity.R`
/// Identifies and removes samples whose pairwise divergence from the majority consensus
/// meets or exceeds `similarity_threshold` (e.g. 0.40).
/// This eliminates off-target captures, paralogs, and reverse-complemented sequences.
pub fn filter_sample_similarity(
    taxa: &[String],
    sequences: &[String],
    similarity_threshold: f64,
) -> (Vec<String>, Vec<String>, Vec<String>) {
    if sequences.len() <= 2 {
        return (taxa.to_vec(), sequences.to_vec(), Vec::new());
    }

    let consensus = compute_majority_consensus(sequences, false);
    let distances = pairwise_distance_to_target(sequences, &consensus);

    let mut kept_taxa = Vec::new();
    let mut kept_seqs = Vec::new();
    let mut dropped_taxa = Vec::new();

    for i in 0..taxa.len() {
        let dist = distances.get(i).copied().unwrap_or(0.0);
        if dist < similarity_threshold {
            kept_taxa.push(taxa[i].clone());
            kept_seqs.push(sequences[i].clone());
        } else {
            dropped_taxa.push(taxa[i].clone());
        }
    }

    (kept_taxa, kept_seqs, dropped_taxa)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_similarity_trimming() {
        let taxa = vec![
            "Species_A".to_string(),
            "Species_B".to_string(),
            "Species_C".to_string(),
            "Paralog_Rogue".to_string(),
        ];
        let seqs = vec![
            "ATGCATGCATGC".to_string(),
            "ATGCATGCATGG".to_string(),
            "ATGCATGCATGA".to_string(),
            "GGGGGGGGGGGG".to_string(), // High divergence from ATGC consensus
        ];

        let (kept_t, kept_s, dropped) = filter_sample_similarity(&taxa, &seqs, 0.40);
        assert_eq!(kept_t.len(), 3);
        assert_eq!(dropped, vec!["Paralog_Rogue".to_string()]);
        assert_eq!(kept_s.len(), 3);
    }
}
