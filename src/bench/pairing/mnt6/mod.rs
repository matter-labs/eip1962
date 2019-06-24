use num_bigint::BigUint;
use num_traits::FromPrimitive;
use num_integer::Integer;
use num_traits::Zero;
use crate::field::{U320Repr, new_field, biguint_to_u64_vec};
use crate::fp::Fp;
use crate::traits::{FieldElement};
use crate::extension_towers::fp3::{Fp3, Extension3};
use crate::extension_towers::fp6_as_2_over_3::{Fp6, Extension2Over3};
use num_traits::Num;
use crate::pairings::{frobenius_calculator_fp3, frobenius_calculator_fp6_as_2_over_3};
use crate::weierstrass::{Group};
use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
use crate::weierstrass::cubic_twist::{TwistPoint, WeierstrassCurveTwist};
use crate::pairings::{PairingEngine};
use rust_test::Bencher;

#[bench]
fn bench_mnt6_pairing(b: &mut Bencher) {
    let modulus = BigUint::from_str_radix("475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137", 10).unwrap();
    let base_field = new_field::<U320Repr>("475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137", 10).unwrap();
    let nonres_repr = U320Repr::from(5);
    let mut fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

    let mut extension_3 = Extension3 {
        field: &base_field,
        non_residue: fp_non_residue.clone(),
        frobenius_coeffs_c1: [Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field)],
        frobenius_coeffs_c2: [Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field)]
    };

    let (coeffs_1, coeffs_2) = frobenius_calculator_fp3(modulus.clone(), &extension_3).unwrap();
    extension_3.frobenius_coeffs_c1 = coeffs_1;
    extension_3.frobenius_coeffs_c2 = coeffs_2;

    let one = Fp::one(&base_field);

    let mut fp3_non_residue = Fp3::zero(&extension_3); // non-residue is 13 + 0*u + 0*u^2
    fp3_non_residue.c0 = fp_non_residue;

    let f_c1 = [Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field),
                Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field)];

    let mut extension_6 = Extension2Over3 {
        non_residue: fp3_non_residue,
        field: &extension_3,
        frobenius_coeffs_c1: f_c1
    };

    let coeffs = frobenius_calculator_fp6_as_2_over_3(modulus, &extension_6).unwrap();
    extension_6.frobenius_coeffs_c1 = coeffs;

    let b_fp = BigUint::from_str_radix("106700080510851735677967319632585352256454251201367587890185989362936000262606668469523074", 10).unwrap().to_bytes_be();
    let b_fp = Fp::from_be_bytes(&base_field, &b_fp, true).unwrap();

    let a_fp = Fp::from_repr(&base_field, U320Repr::from(11)).unwrap();

    let mut twist = Fp3::zero(&extension_3);
    twist.c1 = one.clone();

    let mut twist_squared = twist.clone();
    twist_squared.square();

    let mut twist_cubed = twist_squared.clone();
    twist_cubed.mul_assign(&twist);

    let mut a_fp3 = twist_squared.clone();
    a_fp3.mul_by_fp(&a_fp);

    let mut b_fp3 = twist_cubed.clone();
    b_fp3.mul_by_fp(&b_fp);

    let scalar_field = new_field::<U320Repr>("475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081", 10).unwrap();

    let curve = WeierstrassCurve::new(&scalar_field, a_fp, b_fp);
    let curve_twist = WeierstrassCurveTwist::new(&scalar_field, &extension_3, a_fp3, b_fp3);

    let p_x = BigUint::from_str_radix("336685752883082228109289846353937104185698209371404178342968838739115829740084426881123453", 10).unwrap().to_bytes_be();
    let p_y = BigUint::from_str_radix("402596290139780989709332707716568920777622032073762749862342374583908837063963736098549800", 10).unwrap().to_bytes_be();

    let q_x_0 = BigUint::from_str_radix("421456435772811846256826561593908322288509115489119907560382401870203318738334702321297427", 10).unwrap().to_bytes_be();
    let q_x_1 = BigUint::from_str_radix("103072927438548502463527009961344915021167584706439945404959058962657261178393635706405114", 10).unwrap().to_bytes_be();
    let q_x_2 = BigUint::from_str_radix("143029172143731852627002926324735183809768363301149009204849580478324784395590388826052558", 10).unwrap().to_bytes_be();
    
    let q_y_0 = BigUint::from_str_radix("464673596668689463130099227575639512541218133445388869383893594087634649237515554342751377", 10).unwrap().to_bytes_be();
    let q_y_1 = BigUint::from_str_radix("100642907501977375184575075967118071807821117960152743335603284583254620685343989304941678", 10).unwrap().to_bytes_be();
    let q_y_2 = BigUint::from_str_radix("123019855502969896026940545715841181300275180157288044663051565390506010149881373807142903", 10).unwrap().to_bytes_be();

    let p_x = Fp::from_be_bytes(&base_field, &p_x, true).unwrap();
    let p_y = Fp::from_be_bytes(&base_field, &p_y, true).unwrap();

    let q_x_0 = Fp::from_be_bytes(&base_field, &q_x_0, true).unwrap();
    let q_x_1 = Fp::from_be_bytes(&base_field, &q_x_1, true).unwrap();
    let q_x_2 = Fp::from_be_bytes(&base_field, &q_x_2, true).unwrap();

    let q_y_0 = Fp::from_be_bytes(&base_field, &q_y_0, true).unwrap();
    let q_y_1 = Fp::from_be_bytes(&base_field, &q_y_1, true).unwrap();
    let q_y_2 = Fp::from_be_bytes(&base_field, &q_y_2, true).unwrap();

    let mut q_x = Fp3::zero(&extension_3);
    q_x.c0 = q_x_0;
    q_x.c1 = q_x_1;
    q_x.c2 = q_x_2;

    let mut q_y = Fp3::zero(&extension_3);
    q_y.c0 = q_y_0;
    q_y.c1 = q_y_1;
    q_y.c2 = q_y_2;

    let p = CurvePoint::point_from_xy(&curve, p_x, p_y);
    let q = TwistPoint::point_from_xy(&curve_twist, q_x, q_y);

    let x: Vec<u64> = vec![
        0xdc9a1b671660000, 0x46609756bec2a33f, 0x1eef55
    ];

    assert!(p.check_on_curve());
    assert!(q.check_on_curve());

    let engine = super::MNT6Instance {
        x: x,
        x_is_negative: true,
        exp_w0: vec![0xdc9a1b671660000, 0x46609756bec2a33f, 0x1eef55],
        exp_w1: vec![1u64],
        exp_w0_is_negative: true,
        base_field: &base_field,
        curve: &curve,
        curve_twist: &curve_twist,
        twist: twist,
        fp3_extension: &extension_3,
        fp6_extension: &extension_6,
    };

    b.iter(|| {
        engine.pair(&[p.clone()], &[q.clone()]).unwrap();
    });
}