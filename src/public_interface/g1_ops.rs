/// This api should consist of 
/// - Point decompression
/// - Addition
/// - Multiplication
/// - Multiexponentiations
/// 
/// Every call has common parameters (may be redundant):
/// - Lengths of modulus (in bytes)
/// - Field modulus
/// - Curve A
/// - Curve B
/// - Length of a scalar field (curve order) (in bytes)
/// - Curve order
///
/// Assumptions:
/// - one byte for length encoding

use crate::weierstrass::{Group, CurveOverFpParameters};
use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
use crate::representation::ElementRepr;
use crate::multiexp::peppinger;
use crate::field::*;
use super::constants::*;

use super::decode_g1::*;
use super::decode_utils::*;
use super::decode_fp::*;

use crate::errors::ApiError;

pub trait G1Api {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
}

pub struct G1ApiImplementation<FE: ElementRepr> {
    _marker_fe: std::marker::PhantomData<FE>,
}

impl<FE: ElementRepr> G1Api for G1ApiImplementation<FE> {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, _, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (a, b, rest) = parse_ab_in_base_field_from_encoding(&rest, modulus_len, &field)?;
        let (_order_len, order, rest) = parse_group_order_from_encoding(rest)?;

        let fp_params = CurveOverFpParameters::new(&field);

        let curve = WeierstrassCurve::new(&order.as_ref(), a, b, &fp_params).map_err(|_| {
            ApiError::InputError("Curve shape is not supported".to_owned())
        })?;

        let (mut p_0, rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;
        let (p_1, rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;

        if rest.len() != 0 {
            return Err(ApiError::InputError("Input contains garbage at the end".to_owned()));
        }

        if !p_0.is_on_curve() {
            if !crate::features::in_fuzzing_or_gas_metering() {
                return Err(ApiError::InputError(format!("Point 0 is not on curve, file {}, line {}", file!(), line!())));
            }
        }
        if !p_1.is_on_curve() {
            if !crate::features::in_fuzzing_or_gas_metering() {
                return Err(ApiError::InputError(format!("Point 1 is not on curve, file {}, line {}", file!(), line!())));
            }
        }

        p_0.add_assign(&p_1);

        serialize_g1_point(modulus_len, &p_0)   
    }

    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, _, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (a, b, rest) = parse_ab_in_base_field_from_encoding(&rest, modulus_len, &field)?;
        let (order_len, order, rest) = parse_group_order_from_encoding(rest)?;

        let fp_params = CurveOverFpParameters::new(&field);

        let curve = WeierstrassCurve::new(&order.as_ref(), a, b, &fp_params).map_err(|_| {
            ApiError::InputError("Curve shape is not supported".to_owned())
        })?;

        let (p_0, rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;
        let (scalar, rest) = decode_scalar_representation(rest, order_len, &order)?;

        if rest.len() != 0 {
            return Err(ApiError::InputError("Input contains garbage at the end".to_owned()));
        }

        if !p_0.is_on_curve() {
            if !crate::features::in_fuzzing_or_gas_metering() {
                return Err(ApiError::InputError(format!("Point is not on curve, file {}, line {}", file!(), line!())));
            }
        }

        let p = p_0.mul(&scalar);

        serialize_g1_point(modulus_len, &p)   
    }

    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, _, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (a, b, rest) = parse_ab_in_base_field_from_encoding(&rest, modulus_len, &field)?;
        let (order_len, order, rest) = parse_group_order_from_encoding(rest)?;

        let fp_params = CurveOverFpParameters::new(&field);

        let curve = WeierstrassCurve::new(&order.as_ref(), a, b, &fp_params).map_err(|_| {
            ApiError::InputError("Curve shape is not supported".to_owned())
        })?;

        let (num_pairs_encoding, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
        let num_pairs = num_pairs_encoding[0] as usize;

        if num_pairs == 0 {
            return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
        }

        let expected_pair_len = 2*modulus_len + order_len;
        if rest.len() != expected_pair_len * num_pairs {
            return Err(ApiError::InputError("Input length is invalid for number of pairs".to_owned()));
        }

        let mut global_rest = rest;
        let mut bases = Vec::with_capacity(num_pairs);
        let mut scalars = Vec::with_capacity(num_pairs);

        for _ in 0..num_pairs {
            let (p, local_rest) = decode_g1_point_from_xy(global_rest, modulus_len, &curve)?;
            let (scalar, local_rest) = decode_scalar_representation(local_rest, order_len, &order)?;
            if !p.is_on_curve() {
                if !crate::features::in_fuzzing_or_gas_metering() {
                    return Err(ApiError::InputError(format!("Point is not on curve, file {}, line {}", file!(), line!())));
                }
            }
            bases.push(p);
            scalars.push(scalar);
            global_rest = local_rest;
        }

        if global_rest.len() != 0 {
            return Err(ApiError::InputError("Input contains garbage at the end".to_owned()));
        }

        if bases.len() != scalars.len() || bases.len() == 0 {
            if !crate::features::in_gas_metering() {
                return Err(ApiError::InputError(format!("Multiexp with empty input pairs, file {}, line {}", file!(), line!())));
            } else {
                let result = CurvePoint::zero(&curve);
                return serialize_g1_point(modulus_len, &result);
            }
        } 

        let result = peppinger(&bases, scalars);

        serialize_g1_point(modulus_len, &result)   
    }
}

pub struct PublicG1Api;

impl G1Api for PublicG1Api {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (_, modulus, _) = parse_modulus_and_length(&bytes)?;
        let modulus_limbs = num_limbs_for_modulus(&modulus)?;

        let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, G1ApiImplementation, bytes, add_points); 

        result
    }

    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (_, modulus, _) = parse_modulus_and_length(&bytes)?;
        let modulus_limbs = num_limbs_for_modulus(&modulus)?;
        
        let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, G1ApiImplementation, bytes, mul_point); 

        result
    }

    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (_, modulus, _) = parse_modulus_and_length(&bytes)?;
        let modulus_limbs = num_limbs_for_modulus(&modulus)?;

        let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, G1ApiImplementation, bytes, multiexp); 

        result
    }
}