extern crate test as rust_test;

use num_bigint::BigUint;
use num_traits::FromPrimitive;
use num_integer::Integer;
use num_traits::Zero;
use crate::field::{U384Repr, U256Repr, new_field};
use crate::fp::Fp;
use crate::traits::{FieldElement};
use crate::extension_towers::fp2::{Fp2, Extension2};
use crate::extension_towers::fp6_as_3_over_2::{Fp6, Extension3Over2};
use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};
use num_traits::Num;
use crate::pairings::{frobenius_calculator_fp2, frobenius_calculator_fp6_as_3_over_2, frobenius_calculator_fp12};
use crate::weierstrass::{Group};
use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
use crate::weierstrass::twist::{TwistPoint, WeierstrassCurveTwist};
use crate::pairings::{PairingEngine};
use crate::representation::ElementRepr;
use rust_test::Bencher;

use crate::pairings::bn::{TwistType, BnInstance};

#[bench]
fn bench_bn254_pairing(b: &mut Bencher) {
    let modulus = BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let base_field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let scalar_field = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let mut fp_non_residue = Fp::one(&base_field);
    fp_non_residue.negate(); // non-residue is -1

    let mut extension_2 = Extension2 {
        field: &base_field,
        non_residue: fp_non_residue,
        frobenius_coeffs_c1: [Fp::zero(&base_field), Fp::zero(&base_field)]
    };

    let coeffs = frobenius_calculator_fp2(&extension_2).unwrap();
    extension_2.frobenius_coeffs_c1 = coeffs;

    let one = Fp::one(&base_field);

    // non-residue is u+9
    let mut fp2_non_residue = Fp2::zero(&extension_2);
    let fp_9_repr = U256Repr::from(9u64);
    let fp_9 = Fp::from_repr(&base_field, fp_9_repr).unwrap(); 
    fp2_non_residue.c0 = fp_9.clone();
    fp2_non_residue.c1 = one.clone();

    let f_c1 = [Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2)];

    let mut extension_6 = Extension3Over2 {
        non_residue: fp2_non_residue.clone(),
        field: &extension_2,
        frobenius_coeffs_c1: f_c1.clone(),
        frobenius_coeffs_c2: f_c1,
    };

    let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6_as_3_over_2(modulus.clone(), &extension_6).unwrap();

    extension_6.frobenius_coeffs_c1 = coeffs_c1;
    extension_6.frobenius_coeffs_c2 = coeffs_c2;

    let f_c1 = [Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2)];

    let mut extension_12 = Extension2Over3Over2 {
        non_residue: Fp6::zero(&extension_6),
        field: &extension_6,
        frobenius_coeffs_c1: f_c1,
    };

    let coeffs = frobenius_calculator_fp12(modulus.clone(), &extension_12).unwrap();
    extension_12.frobenius_coeffs_c1 = coeffs;

    let b_fp = Fp::from_repr(&base_field, U256Repr::from(3)).unwrap();
    // here it's b/(u+9)
    let mut b_fp2 = fp2_non_residue.inverse().unwrap();
    b_fp2.mul_by_fp(&b_fp);

    let a_fp = Fp::zero(&base_field);
    let a_fp2 = Fp2::zero(&extension_2);

    let curve = WeierstrassCurve::new(&scalar_field, a_fp, b_fp);
    let twist = WeierstrassCurveTwist::new(&scalar_field, &extension_2, a_fp2, b_fp2);

    let p_x = BigUint::from_str_radix("1", 10).unwrap().to_bytes_be();
    let p_y = BigUint::from_str_radix("2", 10).unwrap().to_bytes_be();

    let q_x_0 = BigUint::from_str_radix("10857046999023057135944570762232829481370756359578518086990519993285655852781", 10).unwrap().to_bytes_be();
    let q_x_1 = BigUint::from_str_radix("11559732032986387107991004021392285783925812861821192530917403151452391805634", 10).unwrap().to_bytes_be();
    let q_y_0 = BigUint::from_str_radix("8495653923123431417604973247489272438418190587263600148770280649306958101930", 10).unwrap().to_bytes_be();
    let q_y_1 = BigUint::from_str_radix("4082367875863433681332203403145435568316851327593401208105741076214120093531", 10).unwrap().to_bytes_be();

    let p_x = Fp::from_be_bytes(&base_field, &p_x, true).unwrap();
    let p_y = Fp::from_be_bytes(&base_field, &p_y, true).unwrap();

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

    let p = CurvePoint::point_from_xy(&curve, p_x, p_y);
    // println!("P.x = {}", p.x.into_repr());
    let q = TwistPoint::point_from_xy(&twist, q_x, q_y);
    // println!("Q.x = {}", q.x.c0.repr);
    // println!("Q.y = {}", q.y.c0.repr);

    // let x = BigUint::from_str_radix("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507", 10).unwrap();
    // println!("x = {}", x);
    // println!("x = {:x}", x);

    assert!(p.check_on_curve());
    assert!(q.check_on_curve());

    let mut minus_one_over_2 = Fp::one(&base_field);
    minus_one_over_2.negate();
    let mut two = Fp::one(&base_field);
    two.double();
    let two_inv = two.inverse().unwrap();
    minus_one_over_2.mul_assign(&two_inv);

    let non_residue_in_p_minus_one_over_2 = fp2_non_residue.pow(&minus_one_over_2.into_repr());

    let u = U256Repr::from(4965661367192848881);
    let mut six_u_plus_2 = u;
    six_u_plus_2.mul2();
    let two_u = six_u_plus_2;
    six_u_plus_2.mul2();
    six_u_plus_2.add_nocarry(&two_u);
    let mut two = U256Repr::from(1);
    two.mul2();
    six_u_plus_2.add_nocarry(&two);

    // println!("Expected coeff = {:x}", BigUint::from_str_radix("827617134098165717451808940080463277390770457691666780560712143809003953598", 10).unwrap());
    // println!("Expected coeff = {:x}", BigUint::from_str_radix("987776078024262725561041258416387561158070255475504730561661362421251696401", 10).unwrap());
    // println!("Expected coeff = {:x}", BigUint::from_str_radix("2813988028633040066320201189843971639620433430176492766961373503539074898364", 10).unwrap());

    let engine = BnInstance {
        u: vec![4965661367192848881],
        u_is_negative: false,
        six_u_plus_2: six_u_plus_2.0[..2].to_vec(),
        twist_type: TwistType::D,
        base_field: &base_field,
        curve: &curve,
        curve_twist: &twist,
        fp2_extension: &extension_2,
        fp6_extension: &extension_6,
        fp12_extension: &extension_12,
        non_residue_in_p_minus_one_over_2: non_residue_in_p_minus_one_over_2
    };

    b.iter(|| {
        engine.pair(&[p.clone()], &[q.clone()]).unwrap();
    });
}