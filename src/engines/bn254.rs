use crate::field::*;
use crate::traits::*;
use crate::representation::*;
use crate::fp::*;
use crate::extension_towers::fp2::*;
use crate::extension_towers::fp6_as_3_over_2::*;
use crate::extension_towers::fp12_as_2_over3_over_2::*;
use crate::weierstrass::*;
use crate::weierstrass::curve::*;
use crate::pairings::{TwistType};
use crate::integers::MaxFieldUint;
use super::generic::*;
use crate::pairings::bn::*;

type Bn254Engine<'a> = Engine<'a,
    Fp<'a, U256Repr, PrimeField<U256Repr> >,
    Fp2<'a, U256Repr, PrimeField<U256Repr> >,
    Fp12<'a, U256Repr, PrimeField<U256Repr> >,
    CurveOverFpParameters<'a, U256Repr, PrimeField<U256Repr> >,
    CurveOverFp2Parameters<'a, U256Repr, PrimeField<U256Repr> >,
    BnInstance<'a,
        U256Repr, 
        PrimeField<U256Repr>, 
        CurveOverFpParameters<'a, U256Repr, PrimeField<U256Repr> >,
        CurveOverFp2Parameters<'a, U256Repr, PrimeField<U256Repr> >
        >,
    ()
>;
use once_cell::sync::Lazy;

pub static BN254_MODULUS: Lazy<MaxFieldUint> = Lazy::new(|| {
    use num_bigint::BigUint;
    use num_traits::*;

    let modulus = BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    modulus
});

pub static BN254_SUBGROUP_ORDER: Lazy<[u64; 4]> = Lazy::new(|| {
    use num_bigint::BigUint;
    use num_traits::*;

    let group_order = BigUint::from_str_radix("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    let group_order_uint = MaxFieldUint::from_big_endian(&group_order.to_bytes_be());

    let mut group_order = [0u64; 4];
    group_order.copy_from_slice(&group_order_uint.as_ref()[..4]);

    group_order
});

pub static BN254_BASE_FIELD: Lazy<PrimeField<U256Repr>> = Lazy::new(|| {
    field_from_modulus(&*BN254_MODULUS).unwrap()
});

pub static BN254_EXT2_FIELD: Lazy<Extension2<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    let mut fp_non_residue = Fp::one(&*BN254_BASE_FIELD);
    fp_non_residue.negate(); // non-residue is -1

    use num_bigint::BigUint;
    use num_traits::*;

    let modulus = BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    let mut extension_2 = Extension2::new(fp_non_residue);
    extension_2.calculate_frobenius_coeffs(&modulus).expect("must work");

    extension_2
});

pub static BN254_FP_NONRESIDUE: Lazy<Fp<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    let mut fp_non_residue = Fp::one(&*BN254_BASE_FIELD);
    fp_non_residue.negate(); // non-residue is -1

    fp_non_residue
});

pub static BN254_FP2_NONRESIDUE: Lazy<Fp2<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    let one = Fp::one(&*BN254_BASE_FIELD);
    
    // non-residue is u+9
    let mut fp2_non_residue = Fp2::zero(&*BN254_EXT2_FIELD);
    let fp_9_repr = U256Repr::from(9u64);
    let fp_9 = Fp::from_repr(&*BN254_BASE_FIELD, fp_9_repr).unwrap(); 
    fp2_non_residue.c0 = fp_9.clone();
    fp2_non_residue.c1 = one.clone();

    fp2_non_residue
});

pub static BN254_EXT6_FIELD: Lazy<Extension3Over2<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    let one = Fp::one(&*BN254_BASE_FIELD);
    
    // non-residue is u+9
    let mut fp2_non_residue = Fp2::zero(&*BN254_EXT2_FIELD);
    let fp_9_repr = U256Repr::from(9u64);
    let fp_9 = Fp::from_repr(&*BN254_BASE_FIELD, fp_9_repr).unwrap(); 
    fp2_non_residue.c0 = fp_9.clone();
    fp2_non_residue.c1 = one.clone();

    let mut extension_6 = Extension3Over2::new(fp2_non_residue.clone());
    extension_6.calculate_frobenius_coeffs_optimized(&*BN254_MODULUS).expect("must work");

    extension_6
});

pub static BN254_EXT12_FIELD: Lazy<Extension2Over3Over2<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(&*BN254_EXT6_FIELD));
    extension_12.calculate_frobenius_coeffs_optimized(&*BN254_MODULUS).expect("must work");

    extension_12
});


pub static BN254_G1_A_COEFF: Lazy<Fp<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    let a_fp = Fp::zero(&*BN254_BASE_FIELD);

    a_fp
});

pub static BN254_G1_B_COEFF: Lazy<Fp<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    let b_fp = Fp::from_repr(&*BN254_BASE_FIELD, U256Repr::from(3)).unwrap();

    b_fp
});

pub static BN254_G2_A_COEFF: Lazy<Fp2<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    let a_fp2 = Fp2::zero(&*BN254_EXT2_FIELD);

    a_fp2
});

pub static BN254_G2_B_COEFF: Lazy<Fp2<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    // here it's b/(u+9)
    let mut b_fp2 = BN254_FP2_NONRESIDUE.inverse().unwrap();
    b_fp2.mul_by_fp(&*BN254_G1_B_COEFF);

    b_fp2
});

pub static BN254_G1_PARAMS: Lazy<CurveOverFpParameters<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    let fp_params = CurveOverFpParameters::new(&*BN254_BASE_FIELD);

    fp_params
});

pub static BN254_G2_PARAMS: Lazy<CurveOverFp2Parameters<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    let fp2_params = CurveOverFp2Parameters::new(&*BN254_EXT2_FIELD);

    fp2_params
});

pub static BN254_G1_CURVE: Lazy<WeierstrassCurve<'static, CurveOverFpParameters<'static, U256Repr, PrimeField<U256Repr>>>> = Lazy::new(|| {
    let curve = WeierstrassCurve::new(&*BN254_SUBGROUP_ORDER, BN254_G1_A_COEFF.clone(), BN254_G1_B_COEFF.clone(), &*BN254_G1_PARAMS).unwrap();
    
    curve
});

pub static BN254_G2_CURVE: Lazy<WeierstrassCurve<'static, CurveOverFp2Parameters<'static, U256Repr, PrimeField<U256Repr>>>> = Lazy::new(|| {
    let curve = WeierstrassCurve::new(&*BN254_SUBGROUP_ORDER, BN254_G2_A_COEFF.clone(), BN254_G2_B_COEFF.clone(), &*BN254_G2_PARAMS).unwrap();
    
    curve
});

pub static BN254_G1_GENERATOR: Lazy<CurvePoint<'static, CurveOverFpParameters<'static, U256Repr, PrimeField<U256Repr>>>> = Lazy::new(|| {
    use num_bigint::BigUint;
    use num_traits::*;
    
    let p_x = BigUint::from_str_radix("1", 10).unwrap().to_bytes_be();
    let p_y = BigUint::from_str_radix("2", 10).unwrap().to_bytes_be();

    let p_x = Fp::from_be_bytes(&*BN254_BASE_FIELD, &p_x, true).unwrap();
    let p_y = Fp::from_be_bytes(&*BN254_BASE_FIELD, &p_y, true).unwrap();

    let p = CurvePoint::point_from_xy(&*BN254_G1_CURVE, p_x, p_y);

    p
});

pub static BN254_G2_GENERATOR: Lazy<CurvePoint<'static, CurveOverFp2Parameters<'static, U256Repr, PrimeField<U256Repr>>>> = Lazy::new(|| {
    use num_bigint::BigUint;
    use num_traits::*;

    let q_x_0 = BigUint::from_str_radix("10857046999023057135944570762232829481370756359578518086990519993285655852781", 10).unwrap().to_bytes_be();
    let q_x_1 = BigUint::from_str_radix("11559732032986387107991004021392285783925812861821192530917403151452391805634", 10).unwrap().to_bytes_be();
    let q_y_0 = BigUint::from_str_radix("8495653923123431417604973247489272438418190587263600148770280649306958101930", 10).unwrap().to_bytes_be();
    let q_y_1 = BigUint::from_str_radix("4082367875863433681332203403145435568316851327593401208105741076214120093531", 10).unwrap().to_bytes_be();

    let q_x_0 = Fp::from_be_bytes(&*BN254_BASE_FIELD, &q_x_0, true).unwrap();
    let q_x_1 = Fp::from_be_bytes(&*BN254_BASE_FIELD, &q_x_1, true).unwrap();
    let q_y_0 = Fp::from_be_bytes(&*BN254_BASE_FIELD, &q_y_0, true).unwrap();
    let q_y_1 = Fp::from_be_bytes(&*BN254_BASE_FIELD, &q_y_1, true).unwrap();

    let mut q_x = Fp2::zero(&*BN254_EXT2_FIELD);
    q_x.c0 = q_x_0;
    q_x.c1 = q_x_1;

    let mut q_y = Fp2::zero(&*BN254_EXT2_FIELD);
    q_y.c0 = q_y_0;
    q_y.c1 = q_y_1;

    let q = CurvePoint::point_from_xy(&*BN254_G2_CURVE, q_x, q_y);

    q
});

static BN254_FP2_NONRESIDUE_IN_P_MINUS_ONE_OVER_TWO: Lazy<Fp2<'static, U256Repr, PrimeField<U256Repr>>> = Lazy::new(|| {
    let mut minus_one_over_2 = Fp::one(&*BN254_BASE_FIELD);
    minus_one_over_2.negate();
    let mut two = Fp::one(&*BN254_BASE_FIELD);
    two.double();
    let two_inv = two.inverse().unwrap();
    minus_one_over_2.mul_assign(&two_inv);

    let non_residue_in_p_minus_one_over_2 = BN254_FP2_NONRESIDUE.pow(&minus_one_over_2.into_repr());

    non_residue_in_p_minus_one_over_2
});

static BN254_SIX_U_PLUS_TWO: Lazy<[u64; 2]> = Lazy::new(|| {
    let u = U256Repr::from(4965661367192848881);
    let mut six_u_plus_2 = u;
    six_u_plus_2.mul2();
    let two_u = six_u_plus_2;
    six_u_plus_2.mul2();
    six_u_plus_2.add_nocarry(&two_u);
    let mut two = U256Repr::from(1);
    two.mul2();
    six_u_plus_2.add_nocarry(&two);

    let mut six_u_plus_two_array = [0u64; 2];
    six_u_plus_two_array.copy_from_slice(&six_u_plus_2.0[..2]);

    six_u_plus_two_array
});

pub const BN254_U: u64 = 4965661367192848881;

pub static BN254_PAIRING_ENGINE: Lazy<
BnInstance<
    'static,
    U256Repr, 
    PrimeField<U256Repr>, 
    CurveOverFpParameters<'static, U256Repr, PrimeField<U256Repr> >,
    CurveOverFp2Parameters<'static, U256Repr, PrimeField<U256Repr> >
>
> = Lazy::new(|| {
    let engine = BnInstanceParams {
        u: &[BN254_U],
        u_is_negative: false,
        six_u_plus_2: &*BN254_SIX_U_PLUS_TWO,
        twist_type: TwistType::D,
        base_field: &*BN254_BASE_FIELD,
        curve: &*BN254_G1_CURVE,
        curve_twist: &*BN254_G2_CURVE,
        fp2_extension: &*BN254_EXT2_FIELD,
        fp6_extension: &*BN254_EXT6_FIELD,
        fp12_extension: &*BN254_EXT12_FIELD,
        non_residue_in_p_minus_one_over_2: (&*BN254_FP2_NONRESIDUE_IN_P_MINUS_ONE_OVER_TWO).clone(),
        force_no_naf: true
    };

    let engine = BnInstance::from_params(engine);

    engine
});

pub static BN254_ENGINE: Lazy<Bn254Engine<'static>> = Lazy::new(|| {
    let e = Engine::<_, _, _, _, _, _, _> {
        curve: &*BN254_G1_CURVE,
        twist: &*BN254_G2_CURVE,
        g1_generator: BN254_G1_GENERATOR.clone(),
        g2_generator: BN254_G2_GENERATOR.clone(),
        pairing_engine: &*BN254_PAIRING_ENGINE,

        base_params: &*BN254_G1_PARAMS,
        twist_params: &*BN254_G2_PARAMS,
        aux: (),
    };

    e
});