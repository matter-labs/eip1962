extern crate test as rust_test;

use num_bigint::BigUint;
use num_traits::FromPrimitive;
use num_integer::Integer;
use num_traits::Zero;
use crate::field::{U384Repr, U832Repr, new_field};
use crate::fp::Fp;
use crate::traits::{FieldElement};
use crate::traits::ZeroAndOne;
use crate::extension_towers::fp2::{Fp2, Extension2};
use num_traits::Num;
use crate::pairings::*;
use crate::integers::MaxFieldUint;
use crate::square_root::*;

use rust_test::Bencher;

#[bench]
fn bench_bls12_381_sqrt_in_fp(b: &mut Bencher) {
    use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
    use crate::extension_towers::Fp6Fp12FrobeniusBaseElements;

    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1
 
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    let mut extension_2 = Extension2::new(fp_non_residue);
    extension_2.calculate_frobenius_coeffs(&modulus).unwrap();

    let mut el = Fp::one(&base_field);
    el.double();
    el.negate();

    b.iter(|| {
        let _ = sqrt_for_three_mod_four(&el).unwrap();
    });
}

#[bench]
fn bench_bls12_381_sqrt_in_fp_with_modulus_check(b: &mut Bencher) {
    use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
    use crate::extension_towers::Fp6Fp12FrobeniusBaseElements;

    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1
 
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    let mut extension_2 = Extension2::new(fp_non_residue);
    extension_2.calculate_frobenius_coeffs(&modulus).unwrap();

    let mut el = Fp::one(&base_field);
    el.double();
    el.negate();

    b.iter(|| {
        let _ = sqrt(&el, None).unwrap();
    });
}

#[bench]
fn bench_bls12_381_sqrt_in_fp2(b: &mut Bencher) {
    use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
    use crate::extension_towers::Fp6Fp12FrobeniusBaseElements;

    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1
 
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    let mut extension_2 = Extension2::new(fp_non_residue);
    extension_2.calculate_frobenius_coeffs(&modulus).unwrap();

    let mut four = Fp::one(&base_field);
    four.double();
    four.double();

    let mut fp2_el = Fp2::zero(&extension_2);
    fp2_el.c0 = four.clone();
    fp2_el.c1 = four.clone();

    fp2_el.square();

    b.iter(|| {
        let _ = sqrt_for_three_mod_four_ext2(&fp2_el).unwrap();
    });
}


#[bench]
fn bench_bls12_381_sqrt_in_fp_exp_baseline(b: &mut Bencher) {
    use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
    use crate::extension_towers::Fp6Fp12FrobeniusBaseElements;

    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1
 
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    let mut extension_2 = Extension2::new(fp_non_residue);
    extension_2.calculate_frobenius_coeffs(&modulus).unwrap();

    let mut el = Fp::one(&base_field);
    el.double();
    el.negate();

    use crate::field::SizedPrimeField;
    use crate::representation::ElementRepr;

    let mut modulus_minus_three_by_four = *el.field.modulus();
    modulus_minus_three_by_four.shr(2);

    b.iter(|| {
        let _ = el.pow(&modulus_minus_three_by_four.as_ref());
    });
}
