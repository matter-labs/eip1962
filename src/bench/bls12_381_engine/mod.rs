extern crate test as rust_test;
use self::rust_test::Bencher;

use num_bigint::BigUint;
use num_traits::FromPrimitive;
use num_integer::Integer;
use crate::traits::{FieldElement};
use crate::traits::ZeroAndOne;
use num_traits::Num;
use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
use crate::weierstrass::{Group, CurveOverFpParameters, CurveOverFp2Parameters};
use crate::pairings::{PairingEngine, TwistType};
use crate::engines::bls12_381::*;

#[bench]
fn bench_bls12_engine_sugroup_g1_double_and_add(b: &mut Bencher) {
    let point = &BLS12_381_G1_GENERATOR;
    let subgroup = &BLS12_381_SUBGROUP_ORDER;
    b.iter(|| {
        assert!(point.mul(&subgroup[..]).is_zero());
    });
}

#[bench]
fn bench_bls12_engine_g1_double_and_add_worst_case(b: &mut Bencher) {
    let point = &BLS12_381_G1_GENERATOR;
    let worst_case_scalar = [std::u64::MAX; 4];
    b.iter(|| {
        assert!(!point.mul(&worst_case_scalar[..]).is_zero());
    });
}

#[bench]
fn bench_bls12_engine_sugroup_g2_double_and_add(b: &mut Bencher) {
    let point = &BLS12_381_G2_GENERATOR;
    let subgroup = &BLS12_381_SUBGROUP_ORDER;
    b.iter(|| {
        assert!(point.mul(&subgroup[..]).is_zero());
    });
}

#[bench]
fn bench_bls12_engine_g2_double_and_add_worst_case(b: &mut Bencher) {
    let point = &BLS12_381_G2_GENERATOR;
    let worst_case_scalar = [std::u64::MAX; 4];
    b.iter(|| {
        assert!(!point.mul(&worst_case_scalar[..]).is_zero());
    });
}

#[bench]
fn bench_bls12_engine_sugroup_g2_window_3(b: &mut Bencher) {
    let point = &BLS12_381_G2_GENERATOR;
    let subgroup = &BLS12_381_SUBGROUP_ORDER;
    b.iter(|| {
        assert!(point.wnaf_mul_with_window_size(&subgroup[..], 3).is_zero());
    });
}

#[bench]
fn bench_bls12_engine_sugroup_g2_window_4(b: &mut Bencher) {
    let point = &BLS12_381_G2_GENERATOR;
    let subgroup = &BLS12_381_SUBGROUP_ORDER;
    b.iter(|| {
        assert!(point.wnaf_mul_with_window_size(&subgroup[..], 4).is_zero());
    });
}

#[bench]
fn bench_bls12_engine_sugroup_g2_window_5(b: &mut Bencher) {
    let point = &BLS12_381_G2_GENERATOR;
    let subgroup = &BLS12_381_SUBGROUP_ORDER;
    b.iter(|| {
        assert!(point.wnaf_mul_with_window_size(&subgroup[..], 5).is_zero());
    });
}

#[bench]
fn bench_bls12_engine_sugroup_g2_window_6(b: &mut Bencher) {
    let point = &BLS12_381_G2_GENERATOR;
    let subgroup = &BLS12_381_SUBGROUP_ORDER;
    b.iter(|| {
        assert!(point.wnaf_mul_with_window_size(&subgroup[..], 6).is_zero());
    });
}

#[bench]
fn bench_bls12_engine_sugroup_g2_window_7(b: &mut Bencher) {
    let point = &BLS12_381_G2_GENERATOR;
    let subgroup = &BLS12_381_SUBGROUP_ORDER;
    b.iter(|| {
        assert!(point.wnaf_mul_with_window_size(&subgroup[..], 7).is_zero());
    });
}


#[bench]
fn bench_bls12_engine_pair_2(b: &mut Bencher) {
    let g1_point = BLS12_381_G1_GENERATOR.clone();
    let g2_point = BLS12_381_G2_GENERATOR.clone();
    let g1s = vec![g1_point; 2];
    let g2s = vec![g2_point; 2];
    let pairs = 
    b.iter(|| {
        assert!(BLS12_381_PAIRING_ENGINE.pair(&g1s, &g2s).is_some());
    });
}

#[bench]
fn bench_bls12_engine_pair_4(b: &mut Bencher) {
    let g1_point = BLS12_381_G1_GENERATOR.clone();
    let g2_point = BLS12_381_G2_GENERATOR.clone();
    let g1s = vec![g1_point; 4];
    let g2s = vec![g2_point; 4];
    let pairs = 
    b.iter(|| {
        assert!(BLS12_381_PAIRING_ENGINE.pair(&g1s, &g2s).is_some());
    });
}

#[cfg(feature = "mappings")]
#[bench]
fn bench_bls12_engine_map_fp_to_g1(b: &mut Bencher) {
    let mut x = BLS12_381_FP_ONE.clone();
    x.double();
    x.double();
    x.square();

    b.iter(|| {
        crate::engines::bls12_381::mapping::fp_to_g1(&x).unwrap();
    });
}

#[cfg(feature = "mappings")]
#[bench]
fn bench_bls12_engine_map_fp2_to_g2(b: &mut Bencher) {
    let mut x = BLS12_381_FP2_ONE.clone();
    x.double();
    x.double();
    x.square();

    b.iter(|| {
        crate::engines::bls12_381::mapping::fp2_to_g2(&x).unwrap();
    });
}