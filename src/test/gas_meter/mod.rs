pub(crate) mod bls12;
pub(crate) mod bn;
pub(crate) mod mnt4;
pub(crate) mod mnt6;
pub(crate) mod arithmetic_ops;

mod monte_carlo;

use crate::field::biguint_to_u64_vec;

use num_bigint::BigUint;
use num_traits::Zero;

pub(crate) fn make_x_bit_length_and_hamming_weight(bit_length: usize, hamming_weight: usize) -> BigUint {
    assert!(bit_length > 0);
    assert!(hamming_weight > 0);
    assert!(bit_length >= hamming_weight);
    if bit_length == hamming_weight {
        let mut x = BigUint::from(1u64);
        x <<= bit_length;
        x -= BigUint::from(1u64);
        assert!(!x.is_zero(), "made zero for {} bits and {} hamming", bit_length, hamming_weight);
        assert!(x.bits() == bit_length);
        return x;
    }

    let mut x = BigUint::from(1u64);
    x <<= bit_length - 1;
    for i in 1..hamming_weight {
        let mut tmp = BigUint::from(1u64);
        tmp <<= i - 1;
        x += tmp;
    }

    assert!(!x.is_zero(), "made zero for {} bits and {} hamming", bit_length, hamming_weight);
    assert!(x.bits() == bit_length);

    x
}

pub(crate) fn six_u_plus_two(u: &BigUint, u_is_positive: bool) -> (BigUint, usize, usize) {
    let r = if u_is_positive { 
        BigUint::from(6u64) * u + BigUint::from(2u64)
    } else {
        BigUint::from(6u64) * u - BigUint::from(2u64)
    };

    let r_vec = biguint_to_u64_vec(r.clone());
    let num_bits = r.bits();
    let mut hamming = 0;
    for r in r_vec.into_iter() {
        hamming += r.count_ones();
    }

    let hamming = hamming as usize;

    assert!(hamming <= num_bits);

    (r, num_bits, hamming)
}