pub fn convert_macse_to_n(sequences: &mut [String]) {
    for seq in sequences.iter_mut() {
        if seq.as_bytes().contains(&b'!') {
            *seq = seq.replace('!', "N");
        }
    }
}
