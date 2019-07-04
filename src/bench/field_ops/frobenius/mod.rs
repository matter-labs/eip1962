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

use rust_test::Bencher;

#[bench]
fn bench_cp6_frobenius(b: &mut Bencher) {
    let modulus = BigUint::from_str_radix("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
    let base_field = new_field::<U832Repr>("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
    let nonres_repr = U832Repr::from(13);
    let mut fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

    let mut extension_3 = Extension3::new(fp_non_residue.clone());
    let (coeffs_1, coeffs_2) = frobenius_calculator_fp3(modulus.clone(), &extension_3).unwrap();
    extension_3.frobenius_coeffs_c1 = coeffs_1;
    extension_3.frobenius_coeffs_c2 = coeffs_2;
    extension_3.frobenius_coeffs_are_calculated = true;

    let one = Fp::one(&base_field);

    let mut fp3_non_residue = Fp3::zero(&extension_3); // non-residue is 13 + 0*u + 0*u^2
    fp3_non_residue.c0 = fp_non_residue;

    let mut extension_6 = fp6_as_2_over_3::Extension2Over3::new(fp3_non_residue);

    b.iter(|| {
        frobenius_calculator_fp6_as_2_over_3(modulus.clone(), &extension_6).unwrap()
    });
}

#[bench]
fn bench_bls12_381_frob_fp12(b: &mut Bencher) {
    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1

    let mut extension_2 = Extension2::new(fp_non_residue);

    let coeffs = frobenius_calculator_fp2(&extension_2).unwrap();
    extension_2.frobenius_coeffs_c1 = coeffs;
    extension_2.frobenius_coeffs_are_calculated = true;

    let one = Fp::one(&base_field);

    let mut fp2_non_residue = Fp2::zero(&extension_2);
    fp2_non_residue.c0 = one.clone();
    fp2_non_residue.c1 = one.clone();

    let mut extension_6 = fp6_as_3_over_2::Extension3Over2::new(fp2_non_residue);

    let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6_as_3_over_2(modulus.clone(), &extension_6).unwrap();

    extension_6.frobenius_coeffs_c1 = coeffs_c1;
    extension_6.frobenius_coeffs_c2 = coeffs_c2;
    extension_6.frobenius_coeffs_are_calculated = true;

    let mut fp2_non_residue = Fp2::zero(&extension_2);

    let mut extension_12 = Extension2Over3Over2::new(fp6_as_3_over_2::Fp6::zero(&extension_6));

    b.iter(|| frobenius_calculator_fp12(modulus.clone(), &extension_12).unwrap());
}

#[bench]
fn bench_bls12_381_frob_fp12_using_sliding(b: &mut Bencher) {
    use crate::sliding_window_exp::{WindowExpBase, IntoWindows};
    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1

    let mut extension_2 = Extension2::new(fp_non_residue);

    let coeffs = frobenius_calculator_fp2(&extension_2).unwrap();
    extension_2.frobenius_coeffs_c1 = coeffs;
    extension_2.frobenius_coeffs_are_calculated = true;

    let one = Fp::one(&base_field);

    let mut fp2_non_residue = Fp2::zero(&extension_2);
    fp2_non_residue.c0 = one.clone();
    fp2_non_residue.c1 = one.clone();

    let base = WindowExpBase::new(&fp2_non_residue, Fp2::one(&extension_2), 15, 7);

    let mut extension_6 = fp6_as_3_over_2::Extension3Over2::new(fp2_non_residue);

    let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6_as_3_over_2_using_sliding_window(modulus.clone(), &base, &extension_6).unwrap();

    extension_6.frobenius_coeffs_c1 = coeffs_c1;
    extension_6.frobenius_coeffs_c2 = coeffs_c2;
    extension_6.frobenius_coeffs_are_calculated = true;

    let mut fp2_non_residue = Fp2::zero(&extension_2);

    let mut extension_12 = Extension2Over3Over2::new(fp6_as_3_over_2::Fp6::zero(&extension_6));

    b.iter(|| frobenius_calculator_fp12_using_sliding_window(modulus.clone(), &base, &extension_12).unwrap());
}