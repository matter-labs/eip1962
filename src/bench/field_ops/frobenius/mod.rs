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
use crate::extension_towers::fp3::{Fp3, Extension3};
use crate::extension_towers::fp6_as_2_over_3;
use crate::extension_towers::fp6_as_3_over_2;
use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};
use num_traits::Num;
use crate::pairings::*;
use crate::integers::MaxFieldUint;

use rust_test::Bencher;

#[bench]
fn bench_cp6_frobenius_optimized(b: &mut Bencher) {
    let modulus = BigUint::from_str_radix("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
    let base_field = new_field::<U832Repr>("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
    let nonres_repr = U832Repr::from(13);
    let mut fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    let mut extension_3 = Extension3::new(fp_non_residue.clone());
    extension_3.calculate_frobenius_coeffs_optimized(&modulus).unwrap();

    let one = Fp::one(&base_field);

    let mut fp3_non_residue = Fp3::zero(&extension_3); // non-residue is 13 + 0*u + 0*u^2
    fp3_non_residue.c0 = fp_non_residue;

    b.iter(|| {
        let mut extension_6 = fp6_as_2_over_3::Extension2Over3::new(fp3_non_residue.clone());
        extension_6.calculate_frobenius_coeffs_optimized(&modulus).unwrap()
    });
}

// #[bench]
// fn bench_bls12_381_frob_fp12_using_sliding(b: &mut Bencher) {
//     use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
//     let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
//     let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
//     let mut fp_non_residue = Fp::one(&base_field);
//     fp_non_residue.negate(); // non-residue is -1
 
//     let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

//     let mut extension_2 = Extension2::new(fp_non_residue);
//     extension_2.calculate_frobenius_coeffs(&modulus).unwrap();

//     let one = Fp::one(&base_field);

//     let mut fp2_non_residue = Fp2::zero(&extension_2);
//     fp2_non_residue.c0 = one.clone();
//     fp2_non_residue.c1 = one.clone();

//     let base = WindowExpBase::new(&fp2_non_residue, Fp2::one(&extension_2), 15, 7);

//     let mut extension_6 = fp6_as_3_over_2::Extension3Over2::new(fp2_non_residue);
//     extension_6.calculate_frobenius_coeffs(&modulus, &base).unwrap();

//     b.iter(|| {
//         let mut extension_12 = Extension2Over3Over2::new(fp6_as_3_over_2::Fp6::zero(&extension_6));
//         extension_12.calculate_frobenius_coeffs(&modulus, &base).unwrap()
//     });
// }

#[bench]
fn bench_bls12_381_frob_fp12_using_optimized_arith(b: &mut Bencher) {
    use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1
 
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    let mut extension_2 = Extension2::new(fp_non_residue);
    extension_2.calculate_frobenius_coeffs(&modulus).unwrap();

    let one = Fp::one(&base_field);

    let mut fp2_non_residue = Fp2::zero(&extension_2);
    fp2_non_residue.c0 = one.clone();
    fp2_non_residue.c1 = one.clone();

    let mut extension_6 = fp6_as_3_over_2::Extension3Over2::new(fp2_non_residue);
    extension_6.calculate_frobenius_coeffs_optimized(&modulus).unwrap();

    b.iter(|| {
        let mut extension_12 = Extension2Over3Over2::new(fp6_as_3_over_2::Fp6::zero(&extension_6));
        extension_12.calculate_frobenius_coeffs_optimized(&modulus).unwrap()
    });
}

// #[bench]
// fn bench_bls12_381_frob_fp6_using_sliding(b: &mut Bencher) {
//     use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
//     let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
//     let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
//     let mut fp_non_residue = Fp::one(&base_field);
//     fp_non_residue.negate(); // non-residue is -1
 
//     let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

//     let mut extension_2 = Extension2::new(fp_non_residue);
//     extension_2.calculate_frobenius_coeffs(&modulus).unwrap();

//     let one = Fp::one(&base_field);

//     let mut fp2_non_residue = Fp2::zero(&extension_2);
//     fp2_non_residue.c0 = one.clone();
//     fp2_non_residue.c1 = one.clone();

//     let base = WindowExpBase::new(&fp2_non_residue, Fp2::one(&extension_2), 15, 7);

//     b.iter(|| {
//         let mut extension_6 = fp6_as_3_over_2::Extension3Over2::new(fp2_non_residue.clone());
//         extension_6.calculate_frobenius_coeffs(&modulus, &base).unwrap();
//     });
// }

#[bench]
fn bench_bls12_381_frob_fp6_using_optimized_arith(b: &mut Bencher) {
    use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1
 
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    let mut extension_2 = Extension2::new(fp_non_residue);
    extension_2.calculate_frobenius_coeffs(&modulus).unwrap();

    let one = Fp::one(&base_field);

    let mut fp2_non_residue = Fp2::zero(&extension_2);
    fp2_non_residue.c0 = one.clone();
    fp2_non_residue.c1 = one.clone();

    b.iter(|| {
        let mut extension_6 = fp6_as_3_over_2::Extension3Over2::new(fp2_non_residue.clone());
        extension_6.calculate_frobenius_coeffs_optimized(&modulus).unwrap();
    });
}

#[bench]
fn bench_bls12_381_frob_fp6_using_base_precomp(b: &mut Bencher) {
    use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
    use crate::extension_towers::Fp6Fp12FrobeniusBaseElements;

    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1
 
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    let mut extension_2 = Extension2::new(fp_non_residue);
    extension_2.calculate_frobenius_coeffs(&modulus).unwrap();

    let one = Fp::one(&base_field);

    let mut fp2_non_residue = Fp2::zero(&extension_2);
    fp2_non_residue.c0 = one.clone();
    fp2_non_residue.c1 = one.clone();

    let precomp_base = Fp6Fp12FrobeniusBaseElements::construct(&modulus, &fp2_non_residue).unwrap();

    b.iter(|| {
        let mut extension_6 = fp6_as_3_over_2::Extension3Over2::new(fp2_non_residue.clone());
        extension_6.calculate_frobenius_coeffs_with_precomp(&precomp_base).unwrap();
    });
}

#[bench]
fn bench_bls12_381_frob_fp12_using_base_precomp(b: &mut Bencher) {
    use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
    use crate::extension_towers::Fp6Fp12FrobeniusBaseElements;

    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1
 
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    let mut extension_2 = Extension2::new(fp_non_residue);
    extension_2.calculate_frobenius_coeffs(&modulus).unwrap();

    let one = Fp::one(&base_field);

    let mut fp2_non_residue = Fp2::zero(&extension_2);
    fp2_non_residue.c0 = one.clone();
    fp2_non_residue.c1 = one.clone();

    let precomp_base = Fp6Fp12FrobeniusBaseElements::construct(&modulus, &fp2_non_residue).unwrap();

    let mut extension_6 = fp6_as_3_over_2::Extension3Over2::new(fp2_non_residue.clone());
    extension_6.calculate_frobenius_coeffs_with_precomp(&precomp_base).unwrap();

    b.iter(|| {
        let mut extension_12 = Extension2Over3Over2::new(fp6_as_3_over_2::Fp6::zero(&extension_6));
        extension_12.calculate_frobenius_coeffs_with_precomp(&precomp_base).unwrap();
    });
}

#[bench]
fn bench_bls12_381_frob_base_precomp_time(b: &mut Bencher) {
    use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
    use crate::extension_towers::Fp6Fp12FrobeniusBaseElements;

    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1
 
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    let mut extension_2 = Extension2::new(fp_non_residue);
    extension_2.calculate_frobenius_coeffs(&modulus).unwrap();

    let one = Fp::one(&base_field);

    let mut fp2_non_residue = Fp2::zero(&extension_2);
    fp2_non_residue.c0 = one.clone();
    fp2_non_residue.c1 = one.clone();

    b.iter(|| {
        let precomp_base = Fp6Fp12FrobeniusBaseElements::construct(&modulus, &fp2_non_residue).unwrap();
    });
}