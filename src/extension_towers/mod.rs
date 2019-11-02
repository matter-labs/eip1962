pub mod fp2;
pub mod fp3;
pub mod fp4_as_2_over_2;
pub mod fp6_as_2_over_3;
pub mod fp6_as_3_over_2;
pub mod fp12_as_2_over3_over_2;

use crate::fp::Fp;
use crate::field::{SizedPrimeField};
use crate::traits::FieldElement;
use crate::representation::{ElementRepr};
use crate::traits::ZeroAndOne;
use crate::constants::*;

pub(crate) fn is_non_nth_root<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>
(
    element: & Fp<'a, FE, F>,
    modulus: &MaxFieldUint,
    n: u64
) -> bool {
    if element.is_zero() {
        return false;
    }
    let one = MaxFieldUint::from(1u64);
    let mut power = *modulus;
    power -= one;

    let divisor = MaxFieldUint::from(n);
    let (power, rem) = power.div_mod(divisor);
    if !rem.is_zero() {
        return false;
    }
    let l = element.pow(power.as_ref());
    let one = Fp::one(element.field);
    if l == one {
        return false;
    } else {
        return true;
    }
}

pub(crate) fn is_non_nth_root_fp2<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>
(
    element: & self::fp2::Fp2<'a, FE, F>,
    modulus: &MaxFieldUint,
    n: u64
) -> bool {
    if element.is_zero() {
        return false;
    }
    let mut power = MaxFieldSquaredUint::from(modulus.as_ref());
    power *= power;
    power -= MaxFieldSquaredUint::from(1u64);

    let divisor = MaxFieldSquaredUint::from(n);
    let (power, rem) = power.div_mod(divisor);
    if !rem.is_zero() {
        return false;
    }
    let l = element.pow(power.as_ref());
    let one = self::fp2::Fp2::one(element.extension_field);
    if l == one {
        return false;
    } else {
        return true;
    }
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use num_traits::Num;
    use num_integer::Integer;
    use num_traits::Zero;
    use crate::field::*;
    use crate::fp::Fp;
    use crate::extension_towers::fp2::{Extension2, Fp2};
    use crate::extension_towers::fp6_as_3_over_2::{Extension3Over2, Fp6};
    use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};
    use crate::traits::FieldElement;
    // use crate::pairings::*;
    use super::*;
    use crate::traits::ZeroAndOne;

    // #[test]
    // fn test_bn254_extension() {
    //     let modulus = BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    //     let base_field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    //     let mut fp_non_residue = Fp::one(&base_field);
    //     fp_non_residue.negate(); // non-residue is -1
    //     let modulus_minus_one_by_2 = (modulus.clone() - BigUint::from(1u64)) >> 1;
    //     let not_a_square = is_non_nth_root(&fp_non_residue, modulus.clone(), 2u64);
    //     assert!(not_a_square);

    //     let mut extension_2 = Extension2::new(fp_non_residue);
    //     extension_2.calculate_frobenius_coeffs(modulus.clone());

    //     let one = Fp::one(&base_field);

    //     // non-residue is u+9
    //     let mut fp2_non_residue = Fp2::zero(&extension_2);
    //     let fp_9_repr = U256Repr::from(9u64);
    //     let fp_9 = Fp::from_repr(&base_field, fp_9_repr).unwrap(); 
    //     fp2_non_residue.c0 = fp_9.clone();
    //     fp2_non_residue.c1 = one.clone();

    //     let not_a_cube = is_non_nth_root_fp2(&fp2_non_residue, modulus.clone(), 3u64);
    //     assert!(not_a_cube);

    //     let not_a_square = is_non_nth_root_fp2(&fp2_non_residue, modulus.clone(), 2u64);
    //     assert!(not_a_square);

    //     let not_a_6th_root = is_non_nth_root_fp2(&fp2_non_residue, modulus.clone(), 6u64);
    //     assert!(not_a_6th_root);

    //     let mut extension_6 = Extension3Over2::new(fp2_non_residue.clone());

    //     let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6_as_3_over_2(modulus.clone(), &extension_6).unwrap();

    //     extension_6.frobenius_coeffs_c1 = coeffs_c1;
    //     extension_6.frobenius_coeffs_c2 = coeffs_c2;
    //     extension_6.frobenius_coeffs_are_calculated = true;

    //     let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(&extension_6));
    // }

}