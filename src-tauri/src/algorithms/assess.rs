/// Port of PhyloProcessR `alignmentAssess.R`
/// Evaluates whether an alignment satisfies minimum quality and completeness gates:
/// 1. `min_taxa`: Minimum absolute number of surviving sequences required.
/// 2. `min_taxa_occupancy_percent`: Minimum taxon occupancy percentage across the entire dataset.
/// 3. `min_length`: Minimum length (bp) required.
/// 4. `max_gap_percent`: Maximum overall gap percentage permitted across the locus.
/// 5. `min_pis_count`: Minimum number of parsimony-informative sites required.
/// 6. `min_pis_percent`: Minimum proportion of parsimony-informative sites required.
/// 7. `min_variable_count`: Minimum number of variable sites required.
/// 8. `min_variable_percent`: Minimum proportion of variable sites required.
pub fn assess_alignment(
    num_taxa: usize,
    total_dataset_taxa: usize,
    length: usize,
    gap_percent: f64,
    min_taxa: usize,
    min_taxa_occupancy_percent: f64,
    min_length: usize,
    max_gap_percent: f64,
    pis_count: usize,
    pis_percent: f64,
    min_pis_count: usize,
    min_pis_percent: f64,
    variable_count: usize,
    variable_percent: f64,
    min_variable_count: usize,
    min_variable_percent: f64,
) -> (bool, Vec<String>) {
    let mut fail_reasons = Vec::new();

    if num_taxa == 0 {
        fail_reasons.push("0 surviving taxa (all samples pruned)".to_string());
    } else {
        if min_taxa > 0 && num_taxa < min_taxa {
            fail_reasons.push(format!("Taxa count ({} < min {})", num_taxa, min_taxa));
        }

        if total_dataset_taxa > 0 && min_taxa_occupancy_percent > 0.0 {
            let occupancy_pct = (num_taxa as f64 / total_dataset_taxa as f64) * 100.0;
            if occupancy_pct < min_taxa_occupancy_percent {
                fail_reasons.push(format!(
                    "Taxon occupancy ({:.1}% < min {:.1}%, {}/{} taxa)",
                    occupancy_pct, min_taxa_occupancy_percent, num_taxa, total_dataset_taxa
                ));
            }
        }
    }

    if min_length > 0 && length < min_length {
        fail_reasons.push(format!("Length ({} bp < min {} bp)", length, min_length));
    }

    if max_gap_percent > 0.0 && gap_percent > max_gap_percent {
        fail_reasons.push(format!(
            "Gap percentage ({:.1}% > max {:.1}%)",
            gap_percent, max_gap_percent
        ));
    }

    if min_pis_count > 0 && pis_count < min_pis_count {
        fail_reasons.push(format!(
            "Parsimony-informative sites ({} < min {})",
            pis_count, min_pis_count
        ));
    }

    if min_pis_percent > 0.0 && pis_percent < min_pis_percent {
        fail_reasons.push(format!(
            "Parsimony-informative proportion ({:.1}% < min {:.1}%)",
            pis_percent, min_pis_percent
        ));
    }

    if min_variable_count > 0 && variable_count < min_variable_count {
        fail_reasons.push(format!(
            "Variable sites ({} < min {})",
            variable_count, min_variable_count
        ));
    }

    if min_variable_percent > 0.0 && variable_percent < min_variable_percent {
        fail_reasons.push(format!(
            "Variable-site proportion ({:.1}% < min {:.1}%)",
            variable_percent, min_variable_percent
        ));
    }

    let pass = fail_reasons.is_empty();
    (pass, fail_reasons)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_assess_pass_fail() {
        // 10 taxa present out of 20 total dataset taxa (50% occupancy)
        let (pass1, reasons1) = assess_alignment(
            10, 20, 500, 15.0, 4, 50.0, 100, 50.0, 40, 8.0, 0, 0.0, 60, 12.0, 0, 0.0,
        );
        assert!(pass1);
        assert!(reasons1.is_empty());

        // 8 taxa present out of 20 total dataset taxa (40% occupancy < 50% min required)
        let (pass2, reasons2) = assess_alignment(
            8, 20, 500, 15.0, 4, 50.0, 100, 50.0, 40, 8.0, 0, 0.0, 60, 12.0, 0, 0.0,
        );
        assert!(!pass2);
        assert_eq!(reasons2.len(), 1);
        assert!(reasons2[0].contains("Taxon occupancy"));

        // Low length and high gap
        let (pass3, reasons3) = assess_alignment(
            10, 20, 50, 60.0, 4, 0.0, 100, 50.0, 1, 2.0, 0, 0.0, 2, 4.0, 0, 0.0,
        );
        assert!(!pass3);
        assert_eq!(reasons3.len(), 2);

        let (pass4, reasons4) = assess_alignment(
            10, 20, 500, 10.0, 4, 0.0, 100, 50.0, 12, 2.4, 20, 5.0, 30, 6.0, 0, 0.0,
        );
        assert!(!pass4);
        assert_eq!(reasons4.len(), 2);
        assert!(reasons4.iter().any(|reason| reason.contains("sites")));
        assert!(reasons4.iter().any(|reason| reason.contains("proportion")));

        let (pass5, reasons5) = assess_alignment(
            10, 20, 500, 10.0, 4, 0.0, 100, 50.0, 30, 6.0, 0, 0.0, 12, 2.4, 20, 5.0,
        );
        assert!(!pass5);
        assert_eq!(reasons5.len(), 2);
        assert!(reasons5
            .iter()
            .any(|reason| reason.contains("Variable sites")));
        assert!(reasons5
            .iter()
            .any(|reason| reason.contains("Variable-site proportion")));
    }
}
