// expected structure:

/// Every call has common parameters (may be redundant):
/// - Curve type
/// - Lengths of modulus (in bytes)
/// - Field modulus
/// - Curve A
/// - Curve B
/// - Lengths of group size (in bytes)
/// - Group size
/// - Type specific params
///
/// Assumptions:
/// - one byte for length encoding
/// 
/// 

use crate::weierstrass::curve::WeierstrassCurve;
use crate::weierstrass::{Group, CurveOverFpParameters, CurveOverFp2Parameters, CurveOverFp3Parameters};
use crate::pairings::*;
use crate::pairings::bls12::{Bls12Instance};
use crate::pairings::bn::{BnInstance};
use crate::pairings::mnt4::{MNT4Instance};
use crate::pairings::mnt6::{MNT6Instance};
use crate::representation::{ElementRepr};
use crate::traits::{FieldElement, ZeroAndOne};
use crate::field::biguint_to_u64_vec;
use crate::sliding_window_exp::WindowExpBase;
use crate::extension_towers::*;
use crate::fp::Fp;

use super::decode_g1::*;
use super::decode_utils::*;
use super::decode_fp::*;
use super::decode_g2::*;
use super::constants::*;

#[macro_use]
use super::api_specialization_macro::*;

use crate::errors::ApiError;

pub struct PublicPairingApi;

impl PairingApi for PublicPairingApi {
    fn pair(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        use crate::field::*;
        let (_curve_type, rest) = split(bytes, CURVE_TYPE_LENGTH, "Input should be longer than curve type encoding")?;
        let (_, modulus, _) = parse_modulus_and_length(&rest)?;
        // let (modulus, _, _, _, _order, _, _) = parse_encodings(&rest)?;
        let modulus_limbs = (modulus.bits() / 64) + 1;
        // let order_limbs = (order.bits() / 64) + 1;

        let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, PairingApiImplementation, bytes, pair); 

        result
    }
}

pub trait PairingApi {
    fn pair(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
}

struct PairingApiImplementation<FE: ElementRepr> {
    _marker_fe: std::marker::PhantomData<FE>,
}

impl<FE: ElementRepr> PairingApi for PairingApiImplementation<FE> {
    fn pair(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (curve_type, rest) = split(bytes, CURVE_TYPE_LENGTH, "Input should be longer than curve type encoding")?;

        match curve_type[0] {
            BLS12 => {
                PairingApiImplementation::<FE>::pair_bls12(&rest)
            },
            BN => {
                PairingApiImplementation::<FE>::pair_bn(&rest)
            },
            MNT4 => {
                PairingApiImplementation::<FE>::pair_mnt4(&rest)
            },
            MNT6 => {
                PairingApiImplementation::<FE>::pair_mnt6(&rest)
            },
            _ => {
                return Err(ApiError::InputError("Unknown curve type".to_owned()));
            }
        }
    }
}

impl<FE: ElementRepr>PairingApiImplementation<FE> {
    fn pair_bls12(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        use crate::extension_towers::fp2::{Fp2, Extension2};
        use crate::extension_towers::fp6_as_3_over_2::{Fp6, Extension3Over2};
        use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};

        let (base_field, modulus_len, modulus, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (a_fp, b_fp, rest) = parse_ab_in_base_field_from_encoding(&rest, modulus_len, &base_field)?;
        if !a_fp.is_zero() {
            return Err(ApiError::UnknownParameter("A parameter must be zero for BLS12 curve".to_owned()));
        }
        let (order_repr, _order_len, _order, rest) = parse_group_order_from_encoding(rest)?;
        let fp_params = CurveOverFpParameters::new(&base_field);
        let g1_curve = WeierstrassCurve::new(order_repr.clone(), a_fp, b_fp.clone(), &fp_params);


        // Now we need to expect:
        // - non-residue for Fp2
        // - non-residue for Fp6
        // - twist type M/D
        // - parameter X
        // - sign of X
        // - number of pairs
        // - list of encoded pairs

        let (fp_non_residue, rest) = decode_fp(&rest, modulus_len, &base_field)?;

        {
            let is_not_a_square = is_non_nth_root(&fp_non_residue, modulus.clone(), 2u64);
            if !is_not_a_square {
                return Err(ApiError::InputError(format!("Non-residue for Fp2 is actually a residue file {}, line {}", file!(), line!())));
            }
        }

        // build an extension field
        let mut extension_2 = Extension2::new(fp_non_residue);
        extension_2.calculate_frobenius_coeffs(modulus.clone()).map_err(|_| {
            ApiError::InputError("Failed to calculate Frobenius coeffs for Fp2".to_owned())
        })?;

        let (fp2_non_residue, rest) = decode_fp2(&rest, modulus_len, &extension_2)?;

        {
            let is_not_a_6th_root = is_non_nth_root_fp2(&fp2_non_residue, modulus.clone(), 6u64);
            if !is_not_a_6th_root {
                return Err(ApiError::InputError(format!("Non-residue for Fp6(12) is actually a residue, file {}, line {}", file!(), line!())));
            }
        }

        let (twist_type_encoding, rest) = split(rest, TWIST_TYPE_LENGTH, "Input is not long enough to get twist type")?;

        let twist_type = match twist_type_encoding[0] {
            TWIST_TYPE_D => TwistType::D,
            TWIST_TYPE_M => TwistType::M, 
            _ => {
                return Err(ApiError::UnknownParameter("Unknown twist type supplied".to_owned()));
            },
        };

        let mut extension_6 = Extension3Over2::new(fp2_non_residue.clone());

        let exp_base = WindowExpBase::new(&fp2_non_residue, Fp2::one(&extension_2), 8, 7);

        {
            extension_6.calculate_frobenius_coeffs(modulus.clone(), &exp_base).map_err(|_| {
                ApiError::UnknownParameter("Can not calculate Frobenius coefficients for Fp6".to_owned())
            })?;
        }

        let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(&extension_6));

        {
            extension_12.calculate_frobenius_coeffs(modulus.clone(), &exp_base).map_err(|_| {
                ApiError::InputError("Can not calculate Frobenius coefficients for Fp12".to_owned())
            })?;
        }

        let fp2_non_residue_inv = fp2_non_residue.inverse().ok_or(ApiError::UnexpectedZero("Fp2 non-residue must be invertible".to_owned()))?;

        let b_fp2 = match twist_type {
            TwistType::D => {
                let mut b_fp2 = fp2_non_residue_inv.clone();
                b_fp2.mul_by_fp(&b_fp);

                b_fp2
            },
            TwistType::M => {
                let mut b_fp2 = fp2_non_residue.clone();
                b_fp2.mul_by_fp(&b_fp);

                b_fp2
            },
        };

        let a_fp2 = Fp2::zero(&extension_2);

        let fp2_params = CurveOverFp2Parameters::new(&extension_2);
        let g2_curve = WeierstrassCurve::new(order_repr, a_fp2, b_fp2, &fp2_params);

        let (x, rest) = decode_biguint_with_length(&rest)?;
        let (x_sign, rest) = split(rest, SIGN_ENCODING_LENGTH, "Input is not long enough to get X sign encoding")?;
        let x_is_negative = match x_sign[0] {
            SIGN_PLUS => false,
            SIGN_MINUS => true,
            _ => {
                return Err(ApiError::InputError("X sign is not encoded properly".to_owned()));
            },
        };

        let (num_pairs_encoding, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
        let num_pairs = num_pairs_encoding[0] as usize;

        let mut global_rest = rest;

        let mut g1_points = vec![];
        let mut g2_points = vec![];

        for _ in 0..num_pairs {
            let (g1, rest) = decode_g1_point_from_xy(&global_rest, modulus_len, &g1_curve)?;
            let (g2, rest) = decode_g2_point_from_xy_in_fp2(&rest, modulus_len, &g2_curve)?;

            global_rest = rest;
            if !g1.check_on_curve() || !g2.check_on_curve() {
                return Err(ApiError::InputError("G1 or G2 point is not on curve".to_owned()));
            }

            if !g1.check_correct_subgroup() || !g2.check_correct_subgroup() {
                return Err(ApiError::InputError("G1 or G2 point is not in the expected subgroup".to_owned()));
            }

            g1_points.push(g1);
            g2_points.push(g2);
        }

        let engine = Bls12Instance {
            x: biguint_to_u64_vec(x),
            x_is_negative: x_is_negative,
            twist_type: twist_type,
            base_field: &base_field,
            curve: &g1_curve,
            curve_twist: &g2_curve,
            fp2_extension: &extension_2,
            fp6_extension: &extension_6,
            fp12_extension: &extension_12,
        };

        let pairing_result = engine.pair(&g1_points, &g2_points);

        if pairing_result.is_none() {
            return Err(ApiError::UnknownParameter("Pairing engine returned no value".to_owned()));
        }

        let one_fp12 = Fp12::one(&extension_12);
        let pairing_result = pairing_result.unwrap();
        let result = if pairing_result == one_fp12 {
            vec![1u8]
        } else {
            vec![0u8]
        };

        Ok(result)
    }

    fn pair_bn(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        use num_bigint::BigUint;
        use crate::extension_towers::fp2::{Fp2, Extension2};
        use crate::extension_towers::fp6_as_3_over_2::{Fp6, Extension3Over2};
        use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};

        let (base_field, modulus_len, modulus, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (a_fp, b_fp, rest) = parse_ab_in_base_field_from_encoding(&rest, modulus_len, &base_field)?;
        if !a_fp.is_zero() {
            return Err(ApiError::UnknownParameter("A parameter must be zero for BN curve".to_owned()));
        }
        let (order_repr, _order_len, _order, rest) = parse_group_order_from_encoding(rest)?;
        let fp_params = CurveOverFpParameters::new(&base_field);
        let g1_curve = WeierstrassCurve::new(order_repr.clone(), a_fp, b_fp.clone(), &fp_params);


        // Now we need to expect:
        // - non-residue for Fp2
        // - non-residue for Fp6
        // - twist type M/D
        // - parameter U
        // - sign of U
        // - number of pairs
        // - list of encoded pairs
        // U is used instead of x for convention of go-ethereum people :)

        let (fp_non_residue, rest) = decode_fp(&rest, modulus_len, &base_field)?;

        {
            let is_not_a_square = is_non_nth_root(&fp_non_residue, modulus.clone(), 2u64);
            if !is_not_a_square {
                return Err(ApiError::InputError(format!("Non-residue for Fp2 is actually a residue file {}, line {}", file!(), line!())));
            }
        }

        // build an extension field
        let mut extension_2 = Extension2::new(fp_non_residue);
        extension_2.calculate_frobenius_coeffs(modulus.clone()).map_err(|_| {
            ApiError::InputError("Failed to calculate Frobenius coeffs for Fp2".to_owned())
        })?;

        let (fp2_non_residue, rest) = decode_fp2(&rest, modulus_len, &extension_2)?;

        {
            let is_not_a_6th_root = is_non_nth_root_fp2(&fp2_non_residue, modulus.clone(), 6u64);
            if !is_not_a_6th_root {
                return Err(ApiError::InputError(format!("Non-residue for Fp6(12) is actually a residue, file {}, line {}", file!(), line!())));
            }
        }

        let (twist_type_encoding, rest) = split(rest, TWIST_TYPE_LENGTH, "Input is not long enough to get twist type")?;

        let twist_type = match twist_type_encoding[0] {
            TWIST_TYPE_D => TwistType::D,
            TWIST_TYPE_M => TwistType::M, 
            _ => {
                return Err(ApiError::UnknownParameter("Unknown twist type supplied".to_owned()));
            },
        };

        let mut extension_6 = Extension3Over2::new(fp2_non_residue.clone());

        let exp_base = WindowExpBase::new(&fp2_non_residue, Fp2::one(&extension_2), 8, 7);

        {
            extension_6.calculate_frobenius_coeffs(modulus.clone(), &exp_base).map_err(|_| {
                ApiError::UnknownParameter("Can not calculate Frobenius coefficients for Fp6".to_owned())
            })?;
        }

        let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(&extension_6));

        {
            extension_12.calculate_frobenius_coeffs(modulus.clone(), &exp_base).map_err(|_| {
                ApiError::InputError("Can not calculate Frobenius coefficients for Fp12".to_owned())
            })?;
        }

        let fp2_non_residue_inv = fp2_non_residue.inverse().ok_or(ApiError::UnexpectedZero("Fp2 non-residue must be invertible".to_owned()))?;

        let b_fp2 = match twist_type {
            TwistType::D => {
                let mut b_fp2 = fp2_non_residue_inv.clone();
                b_fp2.mul_by_fp(&b_fp);

                b_fp2
            },
            TwistType::M => {
                let mut b_fp2 = fp2_non_residue.clone();
                b_fp2.mul_by_fp(&b_fp);

                b_fp2
            },
        };

        let a_fp2 = Fp2::zero(&extension_2);

        let fp2_params = CurveOverFp2Parameters::new(&extension_2);
        let g2_curve = WeierstrassCurve::new(order_repr, a_fp2, b_fp2, &fp2_params);

        let (u, rest) = decode_biguint_with_length(&rest)?;
        let (u_sign, rest) = split(rest, SIGN_ENCODING_LENGTH, "Input is not long enough to get X sign encoding")?;
        let u_is_negative = match u_sign[0] {
            SIGN_PLUS => false,
            SIGN_MINUS => true,
            _ => {
                return Err(ApiError::InputError("X sign is not encoded properly".to_owned()));
            },
        };

        let one = BigUint::from(1u64);
        let six = BigUint::from(6u64);
        let two = BigUint::from(2u64);

        let six_u_plus_two = (six * &u) + &two;
        let p_minus_one_over_2 = (modulus.clone() - &one) >> 1;

        let fp2_non_residue_in_p_minus_one_over_2 = fp2_non_residue.pow(&biguint_to_u64_vec(p_minus_one_over_2));

        let (num_pairs_encoding, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
        let num_pairs = num_pairs_encoding[0] as usize;

        let mut global_rest = rest;

        let mut g1_points = vec![];
        let mut g2_points = vec![];

        for _ in 0..num_pairs {
            let (g1, rest) = decode_g1_point_from_xy(&global_rest, modulus_len, &g1_curve)?;
            let (g2, rest) = decode_g2_point_from_xy_in_fp2(&rest, modulus_len, &g2_curve)?;

            global_rest = rest;
            if !g1.check_on_curve() || !g2.check_on_curve() {
                return Err(ApiError::InputError("G1 or G2 point is not on curve".to_owned()));
            }

            if !g1.check_correct_subgroup() || !g2.check_correct_subgroup() {
                return Err(ApiError::InputError("G1 or G2 point is not in the expected subgroup".to_owned()));
            }

            g1_points.push(g1);
            g2_points.push(g2);
        }

        let engine = BnInstance {
            u: biguint_to_u64_vec(u),
            six_u_plus_2: biguint_to_u64_vec(six_u_plus_two),
            u_is_negative: u_is_negative,
            twist_type: twist_type,
            base_field: &base_field,
            curve: &g1_curve,
            curve_twist: &g2_curve,
            fp2_extension: &extension_2,
            fp6_extension: &extension_6,
            fp12_extension: &extension_12,
            non_residue_in_p_minus_one_over_2: fp2_non_residue_in_p_minus_one_over_2
        };

        let pairing_result = engine.pair(&g1_points, &g2_points);

        if pairing_result.is_none() {
            return Err(ApiError::UnknownParameter("Pairing engine returned no value".to_owned()));
        }

        let one_fp12 = Fp12::one(&extension_12);
        let pairing_result = pairing_result.unwrap();
        let result = if pairing_result == one_fp12 {
            vec![1u8]
        } else {
            vec![0u8]
        };

        Ok(result)
    }

    fn pair_mnt6(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        use crate::extension_towers::fp3::{Fp3, Extension3};
        use crate::extension_towers::fp6_as_2_over_3::{Fp6, Extension2Over3};

        let (base_field, modulus_len, modulus, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (a_fp, b_fp, rest) = parse_ab_in_base_field_from_encoding(&rest, modulus_len, &base_field)?;
        if !a_fp.is_zero() {
            return Err(ApiError::UnknownParameter("A parameter must be zero for BLS12 curve".to_owned()));
        }
        let (order_repr, _order_len, _order, rest) = parse_group_order_from_encoding(rest)?;
        let fp_params = CurveOverFpParameters::new(&base_field);
        let g1_curve = WeierstrassCurve::new(order_repr.clone(), a_fp.clone(), b_fp.clone(), &fp_params);

        // Now we need to expect:
        // - non-residue for Fp3
        // now separate Miller loop params
        // - parameter X
        // - sign of X 
        // Final exp params
        // - exp_w0
        // - exp_w1
        // - exp_w0_is_negative
        // - number of pairs
        // - list of encoded pairs

        let (fp_non_residue, rest) = decode_fp(&rest, modulus_len, &base_field)?;

        {
            let is_not_a_root = is_non_nth_root(&fp_non_residue, modulus.clone(), 6u64);
            if !is_not_a_root {
                return Err(ApiError::InputError(format!("Non-residue for Fp3 is actually a residue, file {}, line {}", file!(), line!())));
            }
        }

        // build an extension field
        let mut extension_3 = Extension3::new(fp_non_residue);
        extension_3.calculate_frobenius_coeffs(modulus.clone()).map_err(|_| {
            ApiError::InputError("Failed to calculate Frobenius coeffs for Fp3".to_owned())
        })?;

        let mut extension_6 = Extension2Over3::new(Fp3::zero(&extension_3));

        {
            extension_6.calculate_frobenius_coeffs(modulus.clone()).map_err(|_| {
                ApiError::UnknownParameter("Can not calculate Frobenius coefficients for Fp6".to_owned())
            })?;
        }

        let one = Fp::one(&base_field);

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

        let fp3_params = CurveOverFp3Parameters::new(&extension_3);
        let g2_curve = WeierstrassCurve::new(order_repr, a_fp3, b_fp3, &fp3_params);

        let (x, rest) = decode_biguint_with_length(&rest)?;
        let (x_sign, rest) = split(rest, SIGN_ENCODING_LENGTH, "Input is not long enough to get X sign encoding")?;
        let x_is_negative = match x_sign[0] {
            SIGN_PLUS => false,
            SIGN_MINUS => true,
            _ => {
                return Err(ApiError::InputError("X sign is not encoded properly".to_owned()));
            },
        };

        let (exp_w0, rest) = decode_biguint_with_length(&rest)?;
        let (exp_w1, rest) = decode_biguint_with_length(&rest)?;

        let (exp_w0_sign, rest) = split(rest, SIGN_ENCODING_LENGTH, "Input is not long enough to get exp_w0 sign encoding")?;
        let exp_w0_is_negative = match exp_w0_sign[0] {
            SIGN_PLUS => false,
            SIGN_MINUS => true,
            _ => {
                return Err(ApiError::InputError("Exp_w0 sign is not encoded properly".to_owned()));
            },
        };

        let (num_pairs_encoding, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
        let num_pairs = num_pairs_encoding[0] as usize;

        let mut global_rest = rest;

        let mut g1_points = vec![];
        let mut g2_points = vec![];

        for _ in 0..num_pairs {
            let (g1, rest) = decode_g1_point_from_xy(&global_rest, modulus_len, &g1_curve)?;
            let (g2, rest) = decode_g2_point_from_xy_in_fp3(&rest, modulus_len, &g2_curve)?;

            global_rest = rest;
            if !g1.check_on_curve() || !g2.check_on_curve() {
                return Err(ApiError::InputError("G1 or G2 point is not on curve".to_owned()));
            }

            if !g1.check_correct_subgroup() || !g2.check_correct_subgroup() {
                return Err(ApiError::InputError("G1 or G2 point is not in the expected subgroup".to_owned()));
            }

            g1_points.push(g1);
            g2_points.push(g2);
        }

        let engine = MNT6Instance {
            x: biguint_to_u64_vec(x),
            x_is_negative: x_is_negative,
            exp_w0: biguint_to_u64_vec(exp_w0),
            exp_w1: biguint_to_u64_vec(exp_w1),
            exp_w0_is_negative: exp_w0_is_negative,
            base_field: &base_field,
            curve: &g1_curve,
            curve_twist: &g2_curve,
            twist: twist,
            fp3_extension: &extension_3,
            fp6_extension: &extension_6,
        };

        let pairing_result = engine.pair(&g1_points, &g2_points);

        if pairing_result.is_none() {
            return Err(ApiError::UnknownParameter("Pairing engine returned no value".to_owned()));
        }

        let one_fp6 = Fp6::one(&extension_6);
        let pairing_result = pairing_result.unwrap();
        let result = if pairing_result == one_fp6 {
            vec![1u8]
        } else {
            vec![0u8]
        };

        Ok(result)
    }

    fn pair_mnt4(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        use crate::extension_towers::fp2::{Fp2, Extension2};
        use crate::extension_towers::fp4_as_2_over_2::{Fp4, Extension2Over2};

        let (base_field, modulus_len, modulus, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (a_fp, b_fp, rest) = parse_ab_in_base_field_from_encoding(&rest, modulus_len, &base_field)?;
        if !a_fp.is_zero() {
            return Err(ApiError::UnknownParameter("A parameter must be zero for BLS12 curve".to_owned()));
        }
        let (order_repr, _order_len, _order, rest) = parse_group_order_from_encoding(rest)?;
        let fp_params = CurveOverFpParameters::new(&base_field);
        let g1_curve = WeierstrassCurve::new(order_repr.clone(), a_fp.clone(), b_fp.clone(), &fp_params);

        // Now we need to expect:
        // - non-residue for Fp2
        // now separate Miller loop params
        // - parameter X
        // - sign of X 
        // Final exp params
        // - exp_w0
        // - exp_w1
        // - exp_w0_is_negative
        // - number of pairs
        // - list of encoded pairs

        let (fp_non_residue, rest) = decode_fp(&rest, modulus_len, &base_field)?;

        {
            let is_not_a_root = is_non_nth_root(&fp_non_residue, modulus.clone(), 4u64);
            if !is_not_a_root {
                return Err(ApiError::InputError(format!("Non-residue for Fp2 is actually a residue, file {}, line {}", file!(), line!())));
            }
        }

        // build an extension field
        let mut extension_2 = Extension2::new(fp_non_residue);
        extension_2.calculate_frobenius_coeffs(modulus.clone()).map_err(|_| {
            ApiError::InputError("Failed to calculate Frobenius coeffs for Fp2".to_owned())
        })?;

        let mut extension_4 = Extension2Over2::new(Fp2::zero(&extension_2));

        {
            extension_4.calculate_frobenius_coeffs(modulus.clone()).map_err(|_| {
                ApiError::UnknownParameter("Can not calculate Frobenius coefficients for Fp4".to_owned())
            })?;
        }

        let one = Fp::one(&base_field);

        let mut twist = Fp2::zero(&extension_2);
        twist.c1 = one.clone();

        let mut twist_squared = twist.clone();
        twist_squared.square();

        let mut twist_cubed = twist_squared.clone();
        twist_cubed.mul_assign(&twist);

        let mut a_fp2 = twist_squared.clone();
        a_fp2.mul_by_fp(&a_fp);

        let mut b_fp2 = twist_cubed.clone();
        b_fp2.mul_by_fp(&b_fp);

        let fp2_params = CurveOverFp2Parameters::new(&extension_2);
        let g2_curve = WeierstrassCurve::new(order_repr, a_fp2, b_fp2, &fp2_params);

        let (x, rest) = decode_biguint_with_length(&rest)?;
        let (x_sign, rest) = split(rest, SIGN_ENCODING_LENGTH, "Input is not long enough to get X sign encoding")?;
        let x_is_negative = match x_sign[0] {
            SIGN_PLUS => false,
            SIGN_MINUS => true,
            _ => {
                return Err(ApiError::InputError("X sign is not encoded properly".to_owned()));
            },
        };

        let (exp_w0, rest) = decode_biguint_with_length(&rest)?;
        let (exp_w1, rest) = decode_biguint_with_length(&rest)?;

        let (exp_w0_sign, rest) = split(rest, SIGN_ENCODING_LENGTH, "Input is not long enough to get exp_w0 sign encoding")?;
        let exp_w0_is_negative = match exp_w0_sign[0] {
            SIGN_PLUS => false,
            SIGN_MINUS => true,
            _ => {
                return Err(ApiError::InputError("Exp_w0 sign is not encoded properly".to_owned()));
            },
        };

        let (num_pairs_encoding, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
        let num_pairs = num_pairs_encoding[0] as usize;

        let mut global_rest = rest;

        let mut g1_points = vec![];
        let mut g2_points = vec![];

        for _ in 0..num_pairs {
            let (g1, rest) = decode_g1_point_from_xy(&global_rest, modulus_len, &g1_curve)?;
            let (g2, rest) = decode_g2_point_from_xy_in_fp2(&rest, modulus_len, &g2_curve)?;

            global_rest = rest;
            if !g1.check_on_curve() || !g2.check_on_curve() {
                return Err(ApiError::InputError("G1 or G2 point is not on curve".to_owned()));
            }

            if !g1.check_correct_subgroup() || !g2.check_correct_subgroup() {
                return Err(ApiError::InputError("G1 or G2 point is not in the expected subgroup".to_owned()));
            }

            g1_points.push(g1);
            g2_points.push(g2);
        }

        let engine = MNT4Instance {
            x: biguint_to_u64_vec(x),
            x_is_negative: x_is_negative,
            exp_w0: biguint_to_u64_vec(exp_w0),
            exp_w1: biguint_to_u64_vec(exp_w1),
            exp_w0_is_negative: exp_w0_is_negative,
            base_field: &base_field,
            curve: &g1_curve,
            curve_twist: &g2_curve,
            twist: twist,
            fp2_extension: &extension_2,
            fp4_extension: &extension_4,
        };

        let pairing_result = engine.pair(&g1_points, &g2_points);

        if pairing_result.is_none() {
            return Err(ApiError::UnknownParameter("Pairing engine returned no value".to_owned()));
        }

        let one_fp4 = Fp4::one(&extension_4);
        let pairing_result = pairing_result.unwrap();
        let result = if pairing_result == one_fp4 {
            vec![1u8]
        } else {
            vec![0u8]
        };

        Ok(result)
    }
}