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

use crate::weierstrass::curve;
use crate::weierstrass::twist;
use crate::weierstrass::Group;
use crate::pairings::*;
use crate::pairings::bls12::{Bls12Instance};
use crate::extension_towers::fp2::{Fp2, Extension2};
use crate::extension_towers::fp6_as_3_over_2::{Fp6, Extension3Over2};
use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};
use crate::representation::{ElementRepr};
use crate::traits::FieldElement;
use crate::field::biguint_to_u64_vec;
use crate::sliding_window_exp::WindowExpBase;
use crate::extension_towers::*;

use num_bigint::BigUint;
use num_traits::FromPrimitive;

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
        let (modulus, _, _, _, _order, _, _) = parse_encodings(&rest)?;
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
            _ => {
                unimplemented!("Not implemented for not BLS12 curves");
            }
        }
    }
}

impl<FE: ElementRepr>PairingApiImplementation<FE> {
    fn pair_bls12(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (base_field, modulus_len, modulus, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (a_fp, b_fp, rest) = parse_ab_in_base_field_from_encoding(&rest, modulus_len, &base_field)?;
        if !a_fp.is_zero() {
            return Err(ApiError::UnknownParameter("A parameter must be zero for BLS12 curve".to_owned()));
        }
        let (order_repr, _order_len, _order, rest) = parse_group_order_from_encoding(rest)?;
        let g1_curve = curve::WeierstrassCurve::new(order_repr.clone(), a_fp, b_fp.clone());

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

        // let coeffs = frobenius_calculator_fp2(&extension_2).map_err(|_| {
        //     ApiError::InputError("Failed to calculate Frobenius coeffs for Fp2".to_owned())
        // })?;
        // extension_2.frobenius_coeffs_c1 = coeffs;

        let (fp2_non_residue, rest) = decode_fp2(&rest, modulus_len, &extension_2)?;

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

        // let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6_as_3_over_2_using_sliding_window(modulus.clone(), &exp_base, &extension_6).map_err(|_| {
        //     ApiError::UnknownParameter("Can not calculate Frobenius coefficients for Fp6".to_owned())
        // })?;

        // extension_6.frobenius_coeffs_c1 = coeffs_c1;
        // extension_6.frobenius_coeffs_c2 = coeffs_c2;

        {
            extension_6.calculate_frobenius_coeffs(modulus.clone(), &exp_base).map_err(|_| {
                ApiError::UnknownParameter("Can not calculate Frobenius coefficients for Fp6".to_owned())
            })?;
        }

        let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(&extension_6));

        // let coeffs = frobenius_calculator_fp12_using_sliding_window(modulus, &exp_base, &extension_12).map_err(|_| {
        //     ApiError::InputError("Can not calculate Frobenius coefficients for Fp12".to_owned())
        // })?;
        // extension_12.frobenius_coeffs_c1 = coeffs;

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
        let g2_curve = twist::WeierstrassCurveTwist::new(order_repr, &extension_2, a_fp2, b_fp2);

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
}