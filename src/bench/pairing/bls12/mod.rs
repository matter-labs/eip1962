extern crate test as rust_test;
use self::rust_test::Bencher;

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
use crate::pairings::{PairingEngine, TwistType};
use crate::pairings::bls12::{Bls12Instance};

#[bench]
fn bench_bls12_381_pairing(b: &mut Bencher) {
    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let scalar_field = new_field::<U256Repr>("52435875175126190479447740508185965837690552500527637822603658699938581184513", 10).unwrap();
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

    let mut fp2_non_residue = Fp2::zero(&extension_2);
    fp2_non_residue.c0 = one.clone();
    fp2_non_residue.c1 = one.clone();

    let f_c1 = [Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2)];

    let mut extension_6 = Extension3Over2 {
        non_residue: fp2_non_residue,
        field: &extension_2,
        frobenius_coeffs_c1: f_c1.clone(),
        frobenius_coeffs_c2: f_c1,
    };

    let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6_as_3_over_2(modulus.clone(), &extension_6).unwrap();

    extension_6.frobenius_coeffs_c1 = coeffs_c1;
    extension_6.frobenius_coeffs_c2 = coeffs_c2;

    let mut fp2_non_residue = Fp2::zero(&extension_2);

        let f_c1 = [Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2)];

    let mut extension_12 = Extension2Over3Over2 {
        non_residue: Fp6::zero(&extension_6),
        field: &extension_6,
        frobenius_coeffs_c1: f_c1,
    };


    let coeffs = frobenius_calculator_fp12(modulus, &extension_12).unwrap();
    extension_12.frobenius_coeffs_c1 = coeffs;

    let b_fp = Fp::from_repr(&base_field, U384Repr::from(4)).unwrap();
    let mut b_fp2 = Fp2::zero(&extension_2);
    b_fp2.c0 = b_fp.clone();
    b_fp2.c1 = b_fp.clone();

    let a_fp = Fp::zero(&base_field);
    let a_fp2 = Fp2::zero(&extension_2);

    let curve = WeierstrassCurve::new(&scalar_field, a_fp, b_fp);
    let twist = WeierstrassCurveTwist::new(&scalar_field, &extension_2, a_fp2, b_fp2);

    let p_x = BigUint::from_str_radix("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507", 10).unwrap().to_bytes_be();
    let p_y = BigUint::from_str_radix("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569", 10).unwrap().to_bytes_be();

    let q_x_0 = BigUint::from_str_radix("352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160", 10).unwrap().to_bytes_be();
    let q_x_1 = BigUint::from_str_radix("3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758", 10).unwrap().to_bytes_be();
    let q_y_0 = BigUint::from_str_radix("1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905", 10).unwrap().to_bytes_be();
    let q_y_1 = BigUint::from_str_radix("927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582", 10).unwrap().to_bytes_be();

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

    // let x = BigUint::from_str_radix("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507", 10).unwrap();
    // println!("x = {}", x);
    // println!("x = {:x}", x);

    assert!(p.check_on_curve());
    assert!(q.check_on_curve());

    let bls12_engine = Bls12Instance {
        x: vec![0xd201000000010000],
        x_is_negative: true,
        twist_type: TwistType::M,
        base_field: &base_field,
        curve: &curve,
        curve_twist: &twist,
        fp2_extension: &extension_2,
        fp6_extension: &extension_6,
        fp12_extension: &extension_12,
    };

    b.iter(|| {
        bls12_engine.pair(&[p.clone()], &[q.clone()]).unwrap();
    });
}