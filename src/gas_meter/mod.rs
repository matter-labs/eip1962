const FIELD_MODULUS_COEFF: u64 = 1;
const GROUP_SIZE_COEFF: u64 = 2;
const BLS12_X_HAMMING_COEFF: u64 = 3;
const BLS12_X_POTENTIALLY_NEGATIVE_COEFF: u64 = 4;

pub(crate) fn calculate_hamming_weight(representation: &[u64]) -> u32 {
    let mut weight = 0;
    for el in representation.iter() {
        weight += el.count_ones();
    }

    weight
}

