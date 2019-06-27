use num_bigint::BigUint;
use num_traits::FromPrimitive;
use num_integer::Integer;
use num_traits::Zero;
use crate::field::{U320Repr, new_field, biguint_to_u64_vec};
use crate::fp::Fp;
use crate::traits::{FieldElement};
use crate::extension_towers::fp2::{Fp2, Extension2};
use crate::extension_towers::fp4_as_2_over_2::{Fp4, Extension2Over2};
use num_traits::Num;
use crate::pairings::{frobenius_calculator_fp2, frobenius_calculator_fp4_as_2_over_2};
use crate::weierstrass::{Group};
use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
use crate::weierstrass::twist::{TwistPoint, WeierstrassCurveTwist};
use crate::pairings::{PairingEngine};
use rust_test::Bencher;

#[bench]
fn bench_mnt4_pairing(b: &mut Bencher) {
    let modulus = BigUint::from_str_radix("475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081", 10).unwrap();
    let base_field = new_field::<U320Repr>("475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081", 10).unwrap();
    let nonres_repr = U320Repr::from(17);
    let fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

    let mut extension_2 = Extension2 {
        field: &base_field,
        non_residue: fp_non_residue.clone(),
        frobenius_coeffs_c1: [Fp::zero(&base_field), Fp::zero(&base_field)]
    };

    let coeffs = frobenius_calculator_fp2(&extension_2).unwrap();
    extension_2.frobenius_coeffs_c1 = coeffs;

    let one = Fp::one(&base_field);

    let mut fp3_non_residue = Fp2::zero(&extension_2); // non-residue is 13 + 0*u + 0*u^2
    fp3_non_residue.c0 = fp_non_residue;

    let f_c1 = [Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field),
                Fp::zero(&base_field)];

    let mut extension_4 = Extension2Over2 {
        non_residue: fp3_non_residue,
        field: &extension_2,
        frobenius_coeffs_c1: f_c1
    };

    let coeffs = frobenius_calculator_fp4_as_2_over_2(modulus, &extension_4).unwrap();
    extension_4.frobenius_coeffs_c1 = coeffs;

    let b_fp = BigUint::from_str_radix("423894536526684178289416011533888240029318103673896002803341544124054745019340795360841685", 10).unwrap().to_bytes_be();
    let b_fp = Fp::from_be_bytes(&base_field, &b_fp, true).unwrap();

    let a_fp = Fp::from_repr(&base_field, U320Repr::from(2)).unwrap();

    let mut twist = Fp2::zero(&extension_2);
    twist.c1 = one.clone();

    let mut twist_squared = twist.clone();
    twist_squared.square();

    let mut twist_cubed = twist_squared.clone();
    twist_cubed.mul_assign(&twist);

    let mut a_fp3 = twist_squared.clone();
    a_fp3.mul_by_fp(&a_fp);

    let mut b_fp3 = twist_cubed.clone();
    b_fp3.mul_by_fp(&b_fp);

    // let scalar_field = new_field::<U320Repr>("475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137", 10).unwrap();
    let group_order = BigUint::from_str_radix("475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137", 10).unwrap();
    let group_order = biguint_to_u64_vec(group_order);
    let curve = WeierstrassCurve::new(group_order.clone(), a_fp, b_fp);
    let curve_twist = WeierstrassCurveTwist::new(group_order.clone(), &extension_2, a_fp3, b_fp3);

    let p_x = BigUint::from_str_radix("60760244141852568949126569781626075788424196370144486719385562369396875346601926534016838", 10).unwrap().to_bytes_be();
    let p_y = BigUint::from_str_radix("363732850702582978263902770815145784459747722357071843971107674179038674942891694705904306", 10).unwrap().to_bytes_be();

    let q_x_0 = BigUint::from_str_radix("438374926219350099854919100077809681842783509163790991847867546339851681564223481322252708", 10).unwrap().to_bytes_be();
    let q_x_1 = BigUint::from_str_radix("37620953615500480110935514360923278605464476459712393277679280819942849043649216370485641", 10).unwrap().to_bytes_be();

    let q_y_0 = BigUint::from_str_radix("37437409008528968268352521034936931842973546441370663118543015118291998305624025037512482", 10).unwrap().to_bytes_be();
    let q_y_1 = BigUint::from_str_radix("424621479598893882672393190337420680597584695892317197646113820787463109735345923009077489", 10).unwrap().to_bytes_be();

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
    let q = TwistPoint::point_from_xy(&curve_twist, q_x, q_y);

    let x: Vec<u64> = vec![
        0xdc9a1b671660000, 0x46609756bec2a33f, 0x1eef55
    ];

    assert!(p.check_on_curve());
    assert!(q.check_on_curve());

    let engine = super::MNT4Instance {
        x: biguint_to_u64_vec(BigUint::from_str_radix("689871209842287392837045615510547309923794944", 10).unwrap()),
        x_is_negative: false,
        exp_w0: biguint_to_u64_vec(BigUint::from_str_radix("689871209842287392837045615510547309923794945", 10).unwrap()),
        exp_w1: vec![1u64],
        exp_w0_is_negative: false,
        base_field: &base_field,
        curve: &curve,
        curve_twist: &curve_twist,
        twist: twist,
        fp2_extension: &extension_2,
        fp4_extension: &extension_4,
    };

    b.iter(|| {
        engine.pair(&[p.clone()], &[q.clone()]).unwrap();
    });
}