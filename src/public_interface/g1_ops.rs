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
use crate::weierstrass::curve::{WeierstrassCurve};
use crate::representation::ElementRepr;
use crate::multiexp::peppinger;
use crate::field::*;
use super::constants::*;

// #[macro_use]
// use super::api_specialization_macro::*;

use super::decode_g1::*;
use super::decode_utils::*;

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
        let (order_repr, _order_len, _, rest) = parse_group_order_from_encoding(rest)?;

        let fp_params = CurveOverFpParameters::new(&field);

        let curve = WeierstrassCurve::new(order_repr, a, b, &fp_params);

        let (mut p_0, rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;
        let (p_1, _rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;

        p_0.add_assign(&p_1);

        serialize_g1_point(modulus_len, &p_0)   
    }

    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, _, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (a, b, rest) = parse_ab_in_base_field_from_encoding(&rest, modulus_len, &field)?;
        let (order_repr, order_len, order, rest) = parse_group_order_from_encoding(rest)?;

        let fp_params = CurveOverFpParameters::new(&field);

        let curve = WeierstrassCurve::new(order_repr.clone(), a, b, &fp_params);

        let (p_0, rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;
        let (scalar, _rest) = decode_scalar_representation(rest, order_len, &order, &order_repr)?;

        let p = p_0.mul(&scalar);

        serialize_g1_point(modulus_len, &p)   
    }

    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, _, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (a, b, rest) = parse_ab_in_base_field_from_encoding(&rest, modulus_len, &field)?;
        let (order_repr, order_len, order, rest) = parse_group_order_from_encoding(rest)?;

        let fp_params = CurveOverFpParameters::new(&field);

        let curve = WeierstrassCurve::new(order_repr.clone(), a, b, &fp_params);

        let (num_pairs_encoding, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
        let num_pairs = num_pairs_encoding[0] as usize;

        let expected_pair_len = 2*modulus_len + order_len;
        if rest.len() != expected_pair_len * num_pairs {
            return Err(ApiError::InputError("Input length is invalid for number of pairs".to_owned()));
        }

        // let mut acc = CurvePoint::zero(&curve);

        let mut global_rest = rest;
        let mut pairs = vec![];

        for _ in 0..num_pairs {
            let (p, local_rest) = decode_g1_point_from_xy(global_rest, modulus_len, &curve)?;
            let (scalar, local_rest) = decode_scalar_representation(local_rest, order_len, &order, &order_repr)?;
            pairs.push((p, scalar));
            // acc.add_assign(&p.mul(&scalar));
            global_rest = local_rest;
        }

        let result = peppinger(pairs);

        serialize_g1_point(modulus_len, &result)   
    }
}

pub struct PublicG1Api;

impl G1Api for PublicG1Api {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (_, modulus, _) = parse_modulus_and_length(&bytes)?;
        let modulus_limbs = (modulus.bits() / 64) + 1;

        let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, G1ApiImplementation, bytes, add_points); 

        result
    }

    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (_, modulus, _) = parse_modulus_and_length(&bytes)?;
        let modulus_limbs = (modulus.bits() / 64) + 1;
        
        let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, G1ApiImplementation, bytes, mul_point); 

        result
    }

    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (_, modulus, _) = parse_modulus_and_length(&bytes)?;
        let modulus_limbs = (modulus.bits() / 64) + 1;

        let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, G1ApiImplementation, bytes, multiexp); 

        result
    }
}