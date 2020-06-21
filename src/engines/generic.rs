use crate::traits::*;
use crate::weierstrass::*;
use crate::weierstrass::curve::*;
use crate::pairings::{PairingEngine};

pub struct Engine<'a,
    F: FieldElement + ZeroAndOne + 'a,
    FEXT: FieldElement + ZeroAndOne + 'a,
    GT: FieldElement + 'a,
    C: CurveParameters<BaseFieldElement = F> + 'a,
    CTW: CurveParameters<BaseFieldElement = FEXT> + 'a,
    P: PairingEngine<PairingResult = GT, G1 = CurvePoint<'a, C>, G2 = CurvePoint<'a, CTW>>,
    AUX: Sized + 'a
> {
    pub curve: &'a WeierstrassCurve<'a, C>,
    pub twist: &'a WeierstrassCurve<'a, CTW>,
    pub g1_generator: P::G1,
    pub g2_generator: P::G2,
    pub pairing_engine: &'a P,
    
    pub(crate) base_params: &'a C,
    pub(crate) twist_params: &'a CTW,
    pub(crate) aux: AUX,
}

impl<'a,
F: FieldElement + ZeroAndOne + 'a,
FEXT: FieldElement + ZeroAndOne + 'a,
GT: FieldElement + 'a,
C: CurveParameters<BaseFieldElement = F> + 'a,
CTW: CurveParameters<BaseFieldElement = FEXT> + 'a,
P: PairingEngine<PairingResult = GT, G1 = CurvePoint<'a, C>, G2 = CurvePoint<'a, CTW>>,
AUX: Sized + 'a
> Engine<'a, F, FEXT, GT, C, CTW, P, AUX> {
    pub fn get_base_field_params(&self) -> <F as ZeroAndOne>::Params {
        self.base_params.params()
    }

    pub fn get_extension_field_params(&self) -> <FEXT as ZeroAndOne>::Params {
        self.twist_params.params()
    }
}

// #[cfg(test)]
// mod test {
//     use super::*;
//     use num_bigint::BigUint;
//     use num_traits::Num;
//     use crate::test::{biguint_to_u64_vec};

//     fn construct_field() -> BoxRef<PrimeField<U256Repr>> {
//         let base_field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
//         let base_field = BoxRef::new(Box::from(base_field));

//         base_field
//     }

//     fn construct_ext<'a>(base_field: &'a PrimeField<U256Repr>) -> BoxRef<Extension2<'a, U256Repr, PrimeField<U256Repr>>> {
//         let mut fp_non_residue = Fp::one(base_field);
//         fp_non_residue.negate(); // non-residue is -1

//         let modulus = BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
//         let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

//         let mut extension_2 = Extension2::new(fp_non_residue);
//         extension_2.calculate_frobenius_coeffs(&modulus).expect("must work");

//         let extension_2 = BoxRef::new(Box::from(extension_2));

//         extension_2
//     }

//     // // fn create_bn254_engine() -> Bn254Engine<'static> {
//     // fn create_bn254_engine<'a>() -> (BoxRef<PrimeField<U256Repr>>, Fp<'a, U256Repr, PrimeField<U256Repr>>) {
//     //     let modulus = BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
//     //     let base_field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
//     //     let base_field = BoxRef::new(Box::from(base_field));
//     //     let group_order = BigUint::from_str_radix("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
//     //     let group_order_vec = biguint_to_u64_vec(group_order);
//     //     let mut fp_non_residue = Fp::one(base_field.as_ref());
//     //     fp_non_residue.negate(); // non-residue is -1

//     //     let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

//     //     // let mut extension_2 = Extension2::new(fp_non_residue);
//     //     // extension_2.calculate_frobenius_coeffs(&modulus).expect("must work");

//     //     // let extension_2 = BoxRef::new(Box::from(extension_2));

//     //     let one = Fp::one(&*base_field);

//     //     (base_field, one)



//     //     // // non-residue is u+9
//     //     // let mut fp2_non_residue = Fp2::zero(&*extension_2);
//     //     // let fp_9_repr = U256Repr::from(9u64);
//     //     // let fp_9 = Fp::from_repr(&*base_field, fp_9_repr).unwrap(); 
//     //     // fp2_non_residue.c0 = fp_9.clone();
//     //     // fp2_non_residue.c1 = one.clone();

//     //     // let mut extension_6 = Extension3Over2::new(fp2_non_residue.clone());
//     //     // extension_6.calculate_frobenius_coeffs_optimized(&modulus).expect("must work");
//     //     // let extension_6 = Arc::from(extension_6);

//     //     // let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(extension_6.as_ref()));
//     //     // extension_12.calculate_frobenius_coeffs_optimized(&modulus).expect("must work");
//     //     // let extension_12 = Arc::from(extension_12);

//     //     // let b_fp = Fp::from_repr(base_field.as_ref(), U256Repr::from(3)).unwrap();
//     //     // // here it's b/(u+9)
//     //     // let mut b_fp2 = fp2_non_residue.inverse().unwrap();
//     //     // b_fp2.mul_by_fp(&b_fp);

//     //     // let a_fp = Fp::zero(base_field.as_ref());
//     //     // let a_fp2 = Fp2::zero(extension_2.as_ref());

//     //     // let fp_params = CurveOverFpParameters::new(base_field.as_ref());
//     //     // let fp_params = Arc::from(fp_params);
//     //     // let fp2_params = CurveOverFp2Parameters::new(extension_2.as_ref());
//     //     // let fp2_params = Arc::from(fp2_params);

//     //     // let mut group_order = [0u64; 4];
//     //     // group_order.copy_from_slice(&group_order_vec[..4]);
//     //     // let group_order = Arc::from(group_order);

//     //     // let curve = WeierstrassCurve::new(group_order.as_ref(), a_fp, b_fp, fp_params.as_ref()).unwrap();
//     //     // let twist = WeierstrassCurve::new(group_order.as_ref(), a_fp2, b_fp2, fp2_params.as_ref()).unwrap();

//     //     // let curve = Arc::from(curve);
//     //     // let twist = Arc::from(twist);

//     //     // let p_x = BigUint::from_str_radix("1", 10).unwrap().to_bytes_be();
//     //     // let p_y = BigUint::from_str_radix("2", 10).unwrap().to_bytes_be();

//     //     // let q_x_0 = BigUint::from_str_radix("10857046999023057135944570762232829481370756359578518086990519993285655852781", 10).unwrap().to_bytes_be();
//     //     // let q_x_1 = BigUint::from_str_radix("11559732032986387107991004021392285783925812861821192530917403151452391805634", 10).unwrap().to_bytes_be();
//     //     // let q_y_0 = BigUint::from_str_radix("8495653923123431417604973247489272438418190587263600148770280649306958101930", 10).unwrap().to_bytes_be();
//     //     // let q_y_1 = BigUint::from_str_radix("4082367875863433681332203403145435568316851327593401208105741076214120093531", 10).unwrap().to_bytes_be();

//     //     // let p_x = Fp::from_be_bytes(base_field.as_ref(), &p_x, true).unwrap();
//     //     // let p_y = Fp::from_be_bytes(base_field.as_ref(), &p_y, true).unwrap();

//     //     // let q_x_0 = Fp::from_be_bytes(base_field.as_ref(), &q_x_0, true).unwrap();
//     //     // let q_x_1 = Fp::from_be_bytes(base_field.as_ref(), &q_x_1, true).unwrap();
//     //     // let q_y_0 = Fp::from_be_bytes(base_field.as_ref(), &q_y_0, true).unwrap();
//     //     // let q_y_1 = Fp::from_be_bytes(base_field.as_ref(), &q_y_1, true).unwrap();

//     //     // let mut q_x = Fp2::zero(extension_2.as_ref());
//     //     // q_x.c0 = q_x_0;
//     //     // q_x.c1 = q_x_1;

//     //     // let mut q_y = Fp2::zero(extension_2.as_ref());
//     //     // q_y.c0 = q_y_0;
//     //     // q_y.c1 = q_y_1;

//     //     // let p = CurvePoint::point_from_xy(curve.as_ref(), p_x, p_y);
//     //     // let q = CurvePoint::point_from_xy(twist.as_ref(), q_x, q_y);

//     //     // assert!(p.is_on_curve());
//     //     // assert!(q.is_on_curve());

//     //     // let mut minus_one_over_2 = Fp::one(base_field.as_ref());
//     //     // minus_one_over_2.negate();
//     //     // let mut two = Fp::one(base_field.as_ref());
//     //     // two.double();
//     //     // let two_inv = two.inverse().unwrap();
//     //     // minus_one_over_2.mul_assign(&two_inv);

//     //     // let non_residue_in_p_minus_one_over_2 = fp2_non_residue.pow(&minus_one_over_2.into_repr());

//     //     // let u = U256Repr::from(4965661367192848881);
//     //     // let mut six_u_plus_2 = u;
//     //     // six_u_plus_2.mul2();
//     //     // let two_u = six_u_plus_2;
//     //     // six_u_plus_2.mul2();
//     //     // six_u_plus_2.add_nocarry(&two_u);
//     //     // let mut two = U256Repr::from(1);
//     //     // two.mul2();
//     //     // six_u_plus_2.add_nocarry(&two);

//     //     // let mut six_u_plus_two_array = [0u64; 2];
//     //     // six_u_plus_two_array.copy_from_slice(&six_u_plus_2.0[..2]);
//     //     // let six_u_plus_2 = Arc::from(six_u_plus_two_array);

//     //     // use crate::pairings::bn::*;

//     //     // let aux = (
//     //     //     Arc::clone(&group_order), 
//     //     //     Arc::clone(&six_u_plus_2),
//     //     //     Arc::clone(&base_field),
//     //     //     Arc::clone(&extension_2),
//     //     //     Arc::clone(&extension_6),
//     //     //     Arc::clone(&extension_12)
//     //     // );

//     //     // let e = Engine::<_, _, _, _, _, _, _> {
//     //     //     // base_field: Arc::clone(&base_field),
//     //     //     // ext_field: Arc::clone(&extension_2),
//     //     //     curve: Arc::clone(&curve),
//     //     //     twist: Arc::clone(&twist),
//     //     //     g1_generator: p,
//     //     //     g2_generator: q,
//     //     //     pairing_engine: None,

//     //     //     base_params: Arc::clone(&fp_params),
//     //     //     twist_params: Arc::clone(&fp2_params),
//     //     //     aux: aux,
//     //     // };

//     //     // let engine = BnInstanceParams {
//     //     //     u: &[4965661367192848881],
//     //     //     u_is_negative: false,
//     //     //     six_u_plus_2: e.aux.1.as_ref(),
//     //     //     twist_type: super::TwistType::D,
//     //     //     base_field: e.aux.2.as_ref(),
//     //     //     curve: e.curve.as_ref(),
//     //     //     curve_twist: e.twist.as_ref(),
//     //     //     fp2_extension: e.aux.3.as_ref(),
//     //     //     fp6_extension: e.aux.4.as_ref(),
//     //     //     fp12_extension: e.aux.5.as_ref(),
//     //     //     non_residue_in_p_minus_one_over_2: non_residue_in_p_minus_one_over_2,
//     //     //     force_no_naf: true
//     //     // };

//     //     // let engine = BnInstance::from_params(engine);

//     //     // drop(six_u_plus_2);

//     //     // e
//     // }

//     #[test]
//     fn t() {
//         let base_field = construct_field();
//         let _ = construct_ext(&*base_field);
//         // let _ = create_bn254_engine();
//     }
// }