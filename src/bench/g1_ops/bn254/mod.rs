extern crate test as rust_test;

use crate::field::*;
use crate::fp::Fp;
use crate::weierstrass::curve::*;
use crate::traits::FieldElement;
use rust_test::Bencher;
use crate::multiexp::{ben_coster, ben_coster_wnaf, peppinger};
use crate::weierstrass::Group;

const MULTIEXP_NUM_POINTS: usize = 100;

#[bench]
fn bench_doubling_bn254(b: &mut Bencher) {
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let one = Fp::one(&field);
    let a_coeff = Fp::zero(&field);
    let mut b_coeff = one.clone();
    b_coeff.double();
    b_coeff.add_assign(&one);

    let curve = WeierstrassCurve::new(
        &group, 
        a_coeff, 
        b_coeff);

    let mut two = one.clone();
    two.double();

    let point = CurvePoint::point_from_xy(
        &curve, 
        one, 
        two);

    b.iter(|| point.clone().double());
}

#[bench]
fn bench_addition_bn254(b: &mut Bencher) {
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let one = Fp::one(&field);
    let a_coeff = Fp::zero(&field);
    let mut b_coeff = one.clone();
    b_coeff.double();
    b_coeff.add_assign(&one);

    let curve = WeierstrassCurve::new(
        &group, 
        a_coeff, 
        b_coeff);

    let mut two = one.clone();
    two.double();

    let point = CurvePoint::point_from_xy(
        &curve, 
        one, 
        two);

    let mut another = point.clone();
    another.double();

    b.iter(|| point.clone().add_assign(&another));
}

#[bench]
fn bench_multiplication_bn254(b: &mut Bencher) {
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let one = Fp::one(&field);
    let a_coeff = Fp::zero(&field);
    let mut b_coeff = one.clone();
    b_coeff.double();
    b_coeff.add_assign(&one);

    let curve = WeierstrassCurve::new(
        &group, 
        a_coeff, 
        b_coeff);

    let mut two = one.clone();
    two.double();

    let point = CurvePoint::point_from_xy(
        &curve, 
        one, 
        two);

    // scalar is order - 1
    let scalar = [0x43e1f593f0000000,
                0x2833e84879b97091,
                0xb85045b68181585d,
                0x30644e72e131a029];
    
    b.iter(|| point.mul(&scalar));
}

#[bench]
fn bench_multiplication_bn254_into_affine(b: &mut Bencher) {
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let one = Fp::one(&field);
    let a_coeff = Fp::zero(&field);
    let mut b_coeff = one.clone();
    b_coeff.double();
    b_coeff.add_assign(&one);

    let curve = WeierstrassCurve::new(
        &group, 
        a_coeff, 
        b_coeff);

    let mut two = one.clone();
    two.double();

    let point = CurvePoint::point_from_xy(
        &curve, 
        one, 
        two);

    // scalar is order - 1
    let scalar = [0x43e1f593f0000000,
                0x2833e84879b97091,
                0xb85045b68181585d,
                0x30644e72e131a029];
    
    b.iter(|| point.mul(&scalar).into_xy());
}

#[bench]
fn bench_multiplication_bn254_into_affine_wnaf(b: &mut Bencher) {
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let one = Fp::one(&field);
    let a_coeff = Fp::zero(&field);
    let mut b_coeff = one.clone();
    b_coeff.double();
    b_coeff.add_assign(&one);

    let curve = WeierstrassCurve::new(
        &group, 
        a_coeff, 
        b_coeff);

    let mut two = one.clone();
    two.double();

    let point = CurvePoint::point_from_xy(
        &curve, 
        one, 
        two);

    // scalar is order - 1
    let scalar = U256Repr([0x43e1f593f0000000,
                0x2833e84879b97091,
                0xb85045b68181585d,
                0x30644e72e131a029]);
    
    b.iter(|| point.wnaf_mul_impl(scalar).into_xy());
}

#[bench]
fn bench_multiplication_bn254_g2_into_affine_wnaf(b: &mut Bencher) {
    use num_bigint::BigUint;
    use num_traits::Num;
    use crate::extension_towers::fp2::{Fp2, Extension2};
    use crate::weierstrass::twist::{WeierstrassCurveTwist, TwistPoint};
    let base_field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let scalar_field = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1

    let mut extension_2 = Extension2 {
        field: &base_field,
        non_residue: fp_non_residue,
        frobenius_coeffs_c1: [Fp::zero(&base_field), Fp::zero(&base_field)]
    };

    let one = Fp::one(&base_field);

    // non-residue is u+9
    let mut fp2_non_residue = Fp2::zero(&extension_2);
    let fp_9_repr = U256Repr::from(9u64);
    let fp_9 = Fp::from_repr(&base_field, fp_9_repr).unwrap(); 
    fp2_non_residue.c0 = fp_9.clone();
    fp2_non_residue.c1 = one.clone();

    let b_fp = Fp::from_repr(&base_field, U256Repr::from(3)).unwrap();
    // here it's b/(u+9)
    let mut b_fp2 = fp2_non_residue.inverse().unwrap();
    b_fp2.mul_by_fp(&b_fp);

    let a_fp2 = Fp2::zero(&extension_2);

    let twist = WeierstrassCurveTwist::new(&scalar_field, &extension_2, a_fp2, b_fp2);

    let q_x_0 = BigUint::from_str_radix("10857046999023057135944570762232829481370756359578518086990519993285655852781", 10).unwrap().to_bytes_be();
    let q_x_1 = BigUint::from_str_radix("11559732032986387107991004021392285783925812861821192530917403151452391805634", 10).unwrap().to_bytes_be();
    let q_y_0 = BigUint::from_str_radix("8495653923123431417604973247489272438418190587263600148770280649306958101930", 10).unwrap().to_bytes_be();
    let q_y_1 = BigUint::from_str_radix("4082367875863433681332203403145435568316851327593401208105741076214120093531", 10).unwrap().to_bytes_be();

    let q_x_0 = Fp::from_be_bytes(&base_field, &q_x_0, true).unwrap();
    let q_x_1 = Fp::from_be_bytes(&base_field, &q_x_1, true).unwrap();
    let q_y_0 = Fp::from_be_bytes(&base_field, &q_y_0, true).unwrap();
    let q_y_1 = Fp::from_be_bytes(&base_field, &q_y_1, true).unwrap();

    let mut q_x = Fp2::zero(&extension_2);
    q_x.c0 = q_x_0;
    q_x.c1 = q_x_1;

    let mut q_y = Fp2::zero(&extension_2);
    q_y.c0 = q_y_0;
    q_y.c1 = q_y_1;

    let point = TwistPoint::point_from_xy(&twist, q_x, q_y);

    // scalar is order - 1
    let scalar = U256Repr([0x43e1f593f0000000,
                0x2833e84879b97091,
                0xb85045b68181585d,
                0x30644e72e131a029]);
    
    b.iter(|| point.wnaf_mul_impl(scalar).into_xy());
}

#[bench]
fn bench_field_inverse(b: &mut Bencher) {
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let mut be_repr = vec![0u8; 32];
    be_repr[31] = 7u8;
    let element = Fp::from_be_bytes(&field, &be_repr[..], false).unwrap();
    
    b.iter(|| element.inverse().unwrap());
}

#[bench]
fn bench_field_mont_inverse(b: &mut Bencher) {
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let mut be_repr = vec![0u8; 32];
    be_repr[31] = 7u8;
    let element = Fp::from_be_bytes(&field, &be_repr[..], false).unwrap();
    
    b.iter(|| element.mont_inverse().unwrap());
}

#[bench]
fn bench_ben_coster_bn254(b: &mut Bencher) {
    use crate::representation::ElementRepr;
    use rand::{RngCore, SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let one = Fp::one(&field);
    let a_coeff = Fp::zero(&field);
    let mut b_coeff = one.clone();
    b_coeff.double();
    b_coeff.add_assign(&one);

    let curve = WeierstrassCurve::new(
        &group, 
        a_coeff, 
        b_coeff);

    let mut two = one.clone();
    two.double();

    let point = CurvePoint::point_from_xy(
        &curve, 
        one, 
        two);

    let pairs: Vec<_> = (0..MULTIEXP_NUM_POINTS).map(|_| {
        let mut scalar = U256Repr::default();
        let mut bytes = vec![0u8; 32];
        rng.fill_bytes(&mut bytes[1..]);
        scalar.read_be(& bytes[..]).unwrap();

        (point.clone(), scalar)
    }).collect();

    b.iter(move || ben_coster(pairs.clone()));
}

#[bench]
fn bench_ben_coster_bn254_using_wnaf(b: &mut Bencher) {
    use crate::representation::ElementRepr;
    use rand::{RngCore, SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let one = Fp::one(&field);
    let a_coeff = Fp::zero(&field);
    let mut b_coeff = one.clone();
    b_coeff.double();
    b_coeff.add_assign(&one);

    let curve = WeierstrassCurve::new(
        &group, 
        a_coeff, 
        b_coeff);

    let mut two = one.clone();
    two.double();

    let point = CurvePoint::point_from_xy(
        &curve, 
        one, 
        two);

    let pairs: Vec<_> = (0..MULTIEXP_NUM_POINTS).map(|_| {
        let mut scalar = U256Repr::default();
        let mut bytes = vec![0u8; 32];
        rng.fill_bytes(&mut bytes[1..]);
        scalar.read_be(& bytes[..]).unwrap();

        (point.clone(), scalar)
    }).collect();

    b.iter(move || ben_coster_wnaf(pairs.clone()));
}

#[bench]
fn bench_peppinger_bn254(b: &mut Bencher) {
    use crate::representation::ElementRepr;
    use rand::{RngCore, SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let one = Fp::one(&field);
    let a_coeff = Fp::zero(&field);
    let mut b_coeff = one.clone();
    b_coeff.double();
    b_coeff.add_assign(&one);

    let curve = WeierstrassCurve::new(
        &group, 
        a_coeff, 
        b_coeff);

    let mut two = one.clone();
    two.double();

    let point = CurvePoint::point_from_xy(
        &curve, 
        one, 
        two);

    let pairs: Vec<_> = (0..MULTIEXP_NUM_POINTS).map(|_| {
        let mut scalar = U256Repr::default();
        let mut bytes = vec![0u8; 32];
        rng.fill_bytes(&mut bytes[1..]);
        scalar.read_be(& bytes[..]).unwrap();

        (point.clone(), scalar)
    }).collect();

    b.iter(move || peppinger(pairs.clone()));
}

#[bench]
fn bench_naive_multiexp_bn254(b: &mut Bencher) {
    use crate::representation::ElementRepr;
    use rand::{RngCore, SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let one = Fp::one(&field);
    let a_coeff = Fp::zero(&field);
    let mut b_coeff = one.clone();
    b_coeff.double();
    b_coeff.add_assign(&one);

    let curve = WeierstrassCurve::new(
        &group, 
        a_coeff, 
        b_coeff);

    let mut two = one.clone();
    two.double();

    let point = CurvePoint::point_from_xy(
        &curve, 
        one, 
        two);

    let pairs: Vec<_> = (0..MULTIEXP_NUM_POINTS).map(|_| {
        let mut scalar = U256Repr::default();
        let mut bytes = vec![0u8; 32];
        rng.fill_bytes(&mut bytes[1..]);
        scalar.read_be(& bytes[..]).unwrap();

        (point.clone(), scalar)
    }).collect();


    b.iter(move || {
        let mut pairs: Vec<_> = pairs.iter().map(|el| el.0.mul(el.1)).collect();
        let mut acc = pairs.pop().unwrap();
        while let Some(p) = pairs.pop() {
            acc.add_assign(&p);
        }
    });
}

#[bench]
fn bench_wnaf_multiexp_bn254(b: &mut Bencher) {
    use crate::representation::ElementRepr;
    use rand::{RngCore, SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let one = Fp::one(&field);
    let a_coeff = Fp::zero(&field);
    let mut b_coeff = one.clone();
    b_coeff.double();
    b_coeff.add_assign(&one);

    let curve = WeierstrassCurve::new(
        &group, 
        a_coeff, 
        b_coeff);

    let mut two = one.clone();
    two.double();

    let point = CurvePoint::point_from_xy(
        &curve, 
        one, 
        two);

    let pairs: Vec<_> = (0..MULTIEXP_NUM_POINTS).map(|_| {
        let mut scalar = U256Repr::default();
        let mut bytes = vec![0u8; 32];
        rng.fill_bytes(&mut bytes[1..]);
        scalar.read_be(& bytes[..]).unwrap();

        (point.clone(), scalar)
    }).collect();


    b.iter(move || {
        let mut pairs: Vec<_> = pairs.iter().map(|el| el.0.wnaf_mul_impl(el.1)).collect();
        let mut acc = pairs.pop().unwrap();
        while let Some(p) = pairs.pop() {
            acc.add_assign(&p);
        }
    });
}