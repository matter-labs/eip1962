use crate::weierstrass::{Group, CurveOverFp2Parameters};
use crate::weierstrass::curve::WeierstrassCurve;
use crate::representation::ElementRepr;
use crate::multiexp::peppinger;

use crate::field::*;

use super::decode_utils::*;
use super::decode_g2::*;
use super::decode_g1::*;
use super::constants::*;

use crate::errors::ApiError;

/// Every call has common parameters (may be redundant):
/// - Lengths of modulus (in bytes)
/// - Field modulus
/// - Extension degree (2/3)
/// - Non-residue
/// - Curve A in Fpk
/// - Curve B in Fpk
/// - Length of a scalar field (curve order) (in bytes)
/// - Curve order

pub trait G2Api {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
}

pub struct G2ApiImplementationFp2<FE: ElementRepr> {
    _marker_fe: std::marker::PhantomData<FE>,
}

impl<FE: ElementRepr> G2Api for G2ApiImplementationFp2<FE> {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, modulus, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (extension_2, rest) = create_fp2_extension(rest, modulus, modulus_len, &field)?;
        let (a, b, rest) = parse_ab_in_fp2_from_encoding(&rest, modulus_len, &extension_2)?;
        let (order_repr, _order_len, _, rest) = parse_group_order_from_encoding(rest)?;

        let fp2_params = CurveOverFp2Parameters::new(&extension_2);

        let curve = WeierstrassCurve::new(order_repr, a, b, &fp2_params);

        let (mut p_0, rest) = decode_g2_point_from_xy_in_fp2(rest, modulus_len, &curve)?;
        let (p_1, _rest) = decode_g2_point_from_xy_in_fp2(rest, modulus_len, &curve)?;

        p_0.add_assign(&p_1);

        serialize_g2_point_in_fp2(modulus_len, &p_0)   
    }

    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, modulus, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (extension_2, rest) = create_fp2_extension(rest, modulus, modulus_len, &field)?;
        let (a, b, rest) = parse_ab_in_fp2_from_encoding(&rest, modulus_len, &extension_2)?;
        let (order_repr, order_len, order, rest) = parse_group_order_from_encoding(rest)?;

        let fp2_params = CurveOverFp2Parameters::new(&extension_2);

        let curve = WeierstrassCurve::new(order_repr.clone(), a, b, &fp2_params);

        let (p_0, rest) = decode_g2_point_from_xy_in_fp2(rest, modulus_len, &curve)?;
        let (scalar, _rest) = decode_scalar_representation(rest, order_len, &order, &order_repr)?;

        let p = p_0.mul(&scalar);

        serialize_g2_point_in_fp2(modulus_len, &p)   
    }

    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, modulus, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (extension_2, rest) = create_fp2_extension(&rest, modulus, modulus_len, &field)?;
        let (a, b, rest) = parse_ab_in_fp2_from_encoding(&rest, modulus_len, &extension_2)?;
        let (order_repr, order_len, order, rest) = parse_group_order_from_encoding(rest)?;

        let fp2_params = CurveOverFp2Parameters::new(&extension_2);

        let curve = WeierstrassCurve::new(order_repr.clone(), a, b, &fp2_params);

        let (num_pairs_encoding, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
        let num_pairs = num_pairs_encoding[0] as usize;

        if num_pairs == 0 {
            return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
        }

        let expected_pair_len = 4*modulus_len + order_len;
        if rest.len() != expected_pair_len * num_pairs {
            return Err(ApiError::InputError("Input length is invalid for number of pairs".to_owned()));
        }

        let mut global_rest = rest;
        let mut pairs = vec![];

        for _ in 0..num_pairs {
            let (p, local_rest) = decode_g2_point_from_xy_in_fp2(global_rest, modulus_len, &curve)?;
            let (scalar, local_rest) = decode_scalar_representation(local_rest, order_len, &order, &order_repr)?;
            pairs.push((p, scalar));
            global_rest = local_rest;
        }

        let result = peppinger(pairs);

        serialize_g2_point_in_fp2(modulus_len, &result)   
    }
}

pub struct PublicG2Api;

impl G2Api for PublicG2Api {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (modulus, _, extension_degree, _, _) = parse_modulus_and_extension_degree(&bytes)?;
        let modulus_limbs = num_limbs_for_modulus(&modulus)?;

        let result: Result<Vec<u8>, ApiError> = match extension_degree {
            EXTENSION_DEGREE_2 => {
                let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, G2ApiImplementationFp2, bytes, add_points); 

                result
            },
            EXTENSION_DEGREE_3 => {
                return Err(ApiError::InputError("Extension degree 3 is not yet implemented".to_owned()));
            },
            _ => {
                return Err(ApiError::InputError("Invalid extension degree".to_owned()));
            }
        };

        result
    }

    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (modulus, _, extension_degree, _, _) = parse_modulus_and_extension_degree(&bytes)?;
        let modulus_limbs = num_limbs_for_modulus(&modulus)?;

        let result: Result<Vec<u8>, ApiError> = match extension_degree {
            EXTENSION_DEGREE_2 => {
                let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, G2ApiImplementationFp2, bytes, mul_point); 

                result
            },
            EXTENSION_DEGREE_3 => {
                return Err(ApiError::InputError("Extension degree 3 is not yet implemented".to_owned()));
            },
            _ => {
                return Err(ApiError::InputError("Invalid extension degree".to_owned()));
            }
        };

        result
    }

    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (modulus, _, extension_degree, _, _) = parse_modulus_and_extension_degree(&bytes)?;
        let modulus_limbs = num_limbs_for_modulus(&modulus)?;

        let result: Result<Vec<u8>, ApiError> = match extension_degree {
            EXTENSION_DEGREE_2 => {
                let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, G2ApiImplementationFp2, bytes, multiexp); 

                result
            },
            EXTENSION_DEGREE_3 => {
                return Err(ApiError::InputError("Extension degree 3 is not yet implemented".to_owned()));
            },
            _ => {
                return Err(ApiError::InputError("Invalid extension degree".to_owned()));
            }
        };

        result
    }
}