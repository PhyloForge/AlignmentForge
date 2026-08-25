use crate::algorithms::orf::{translate_codon, GeneticCode};

pub fn apply_codon_qc(
    taxa: &[String],
    sequences: &[String],
    genetic_code: GeneticCode,
    trim_terminal_frameshifts: bool,
    max_internal_macse_sample: usize,
    max_internal_macse_locus: usize,
    max_stop_codons_sample: usize,
    max_stop_codons_locus: usize,
) -> (Vec<String>, Vec<String>, Vec<String>) {
    if sequences.is_empty() {
        return (taxa.to_vec(), sequences.to_vec(), vec![]);
    }

    let mut current_seqs: Vec<Vec<u8>> = sequences.iter().map(|s| s.as_bytes().to_vec()).collect();
    let length = current_seqs[0].len();
    let codon_count = length / 3;

    let mut dropped_taxa = Vec::new();
    let mut kept_taxa = Vec::new();
    let mut kept_seqs = Vec::new();

    let mut total_locus_macse_cols = vec![false; length];
    let mut total_locus_stop_codons = 0;

    for (t_idx, seq) in current_seqs.iter_mut().enumerate() {
        let taxon = &taxa[t_idx];
        
        // 1. Trim terminal frameshifts (!)
        if trim_terminal_frameshifts {
            // Trim 5'
            for c_idx in 0..codon_count {
                let start = c_idx * 3;
                let codon = &seq[start..start + 3];
                if codon.contains(&b'!') {
                    seq[start] = b'-';
                    seq[start + 1] = b'-';
                    seq[start + 2] = b'-';
                } else {
                    break;
                }
            }
            // Trim 3'
            for c_idx in (0..codon_count).rev() {
                let start = c_idx * 3;
                let codon = &seq[start..start + 3];
                if codon.contains(&b'!') {
                    seq[start] = b'-';
                    seq[start + 1] = b'-';
                    seq[start + 2] = b'-';
                } else {
                    break;
                }
            }
        }

        // 2. Count internal frameshifts & stops
        let mut sample_internal_macse = 0;
        let mut sample_stop_codons = 0;

        for c_idx in 0..codon_count {
            let start = c_idx * 3;
            let codon = &seq[start..start + 3];
            
            // Macse
            if codon.contains(&b'!') {
                sample_internal_macse += codon.iter().filter(|&&b| b == b'!').count();
                for i in 0..3 {
                    if codon[i] == b'!' {
                        total_locus_macse_cols[start + i] = true;
                    }
                }
            }
            
            // Stop Codons
            if !codon.contains(&b'-') && !codon.contains(&b'N') && !codon.contains(&b'?') && !codon.contains(&b'!') {
                let aa = translate_codon(codon, genetic_code);
                if aa == '*' && c_idx < codon_count - 1 {
                    sample_stop_codons += 1;
                    total_locus_stop_codons += 1;
                }
            }
        }

        if sample_internal_macse > max_internal_macse_sample || sample_stop_codons > max_stop_codons_sample {
            dropped_taxa.push(taxon.clone());
        } else {
            kept_taxa.push(taxon.clone());
            kept_seqs.push(seq.clone());
        }
    }

    let locus_macse_cols = total_locus_macse_cols.iter().filter(|&&b| b).count();
    
    // Check locus limits
    if locus_macse_cols > max_internal_macse_locus || total_locus_stop_codons > max_stop_codons_locus {
        // Drop the entire locus
        return (vec![], vec![], taxa.to_vec());
    }

    let out_seqs = kept_seqs.into_iter().map(|s| String::from_utf8(s).unwrap()).collect();
    (kept_taxa, out_seqs, dropped_taxa)
}

pub fn convert_macse_to_n(sequences: &mut [String]) {
    for seq in sequences.iter_mut() {
        unsafe {
            let bytes = seq.as_bytes_mut();
            for b in bytes.iter_mut() {
                if *b == b'!' {
                    *b = b'N';
                }
            }
        }
    }
}
