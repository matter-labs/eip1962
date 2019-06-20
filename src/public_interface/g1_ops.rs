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

use crate::weierstrass::Group;
use crate::weierstrass::curve::{WeierstrassCurve};
use crate::field::{field_from_modulus};
use crate::fp::Fp;
use crate::field::{
    U256Repr, 
    U320Repr,
    U384Repr,
    U448Repr,
    U512Repr,
    U576Repr,
    U640Repr,
    U704Repr,
    U768Repr,
    U832Repr,
    U896Repr,
    // U960Repr,
    // U1024Repr
};
use crate::representation::ElementRepr;
use crate::multiexp::peppinger;

use num_bigint::BigUint;

use super::decode_utils::parse_encodings;

#[macro_use]
use super::decode_g1::*;

use crate::errors::ApiError;

pub trait G1Api {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
}

pub struct G1ApiImplementation<FE: ElementRepr, GE: ElementRepr> {
    _marker_fe: std::marker::PhantomData<FE>,
    _marker_ge: std::marker::PhantomData<GE>
}

impl<FE: ElementRepr, GE: ElementRepr> G1Api for G1ApiImplementation<FE, GE> {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, rest) = create_base_field!(bytes, FE);
        let (a, b, rest) = get_ab_in_base_field!(rest, field, modulus_len);
        let (group, order_len, rest) = create_group!(rest, GE);

        let curve = WeierstrassCurve::new(&group, a, b);

        let (mut p_0, rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;
        let (p_1, _rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;

        p_0.add_assign(&p_1);

        serialize_g1_point(modulus_len, &p_0)   
    }

    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, rest) = create_base_field!(bytes, FE);
        let (a, b, rest) = get_ab_in_base_field!(rest, field, modulus_len);
        let (group, order_len, rest) = create_group!(rest, GE);

        let curve = WeierstrassCurve::new(&group, a, b);

        let (p_0, rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;
        let (scalar, _rest) = decode_scalar_representation(rest, order_len, &group)?;

        let p = p_0.mul(&scalar);

        serialize_g1_point(modulus_len, &p)   
    }

    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, rest) = create_base_field!(bytes, FE);
        let (a, b, rest) = get_ab_in_base_field!(rest, field, modulus_len);
        let (group, order_len, rest) = create_group!(rest, GE);

        let curve = WeierstrassCurve::new(&group, a, b);

        let expected_pair_len = 2*modulus_len + order_len;
        if rest.len() % expected_pair_len != 0 {
            return Err(ApiError::InputError("Input length is invalid for number of pairs".to_owned()));
        }

        let expected_pairs = rest.len() / expected_pair_len;
        if expected_pairs == 0 {
            return Err(ApiError::InputError("Number of pairs must be > 0".to_owned()));
        }

        // let mut acc = CurvePoint::zero(&curve);

        let mut global_rest = rest;
        let mut pairs = vec![];

        for _ in 0..expected_pairs {
            let (p, local_rest) = decode_g1_point_from_xy(global_rest, modulus_len, &curve)?;
            let (scalar, local_rest) = decode_scalar_representation(local_rest, order_len, &group)?;
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
        let (modulus, _, _, _, order, _, _) = parse_encodings(&bytes)?;
        let modulus_limbs = (modulus.bits() / 64) + 1;
        let order_limbs = (order.bits() / 64) + 1;

        let result: Result<Vec<u8>, ApiError> = match (modulus_limbs, order_limbs) {
            (4, 4) => {
                G1ApiImplementation::<U256Repr, U256Repr>::add_points(&bytes)
            },
            (5, 4) => {
                G1ApiImplementation::<U320Repr, U256Repr>::add_points(&bytes)
            },
            (5, 5) => {
                G1ApiImplementation::<U320Repr, U320Repr>::add_points(&bytes)
            },
            _ => {
                unimplemented!();
            }
        };

        result
    }

    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (modulus, _, _, _, order, _, _) = parse_encodings(&bytes)?;
        let modulus_limbs = (modulus.bits() / 64) + 1;
        let order_limbs = (order.bits() / 64) + 1;
        
        let result: Result<Vec<u8>, ApiError> = match (modulus_limbs, order_limbs) {
            (4, 4) => {
                G1ApiImplementation::<U256Repr, U256Repr>::mul_point(&bytes)
            },
            (5, 4) => {
                G1ApiImplementation::<U320Repr, U256Repr>::mul_point(&bytes)
            },
            (5, 5) => {
                G1ApiImplementation::<U320Repr, U320Repr>::mul_point(&bytes)
            },
            (10, 7) => {
                G1ApiImplementation::<U640Repr, U448Repr>::mul_point(&bytes)
            },
            (field_limbs, group_limbs) => {
                unimplemented!("unimplemented for {} modulus and {} group limbs", field_limbs, group_limbs);
            }
            // _ => {
            //     unimplemented!();
            // }
        };

        result
    }

    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (modulus, _, _, _, order, _, _) = parse_encodings(&bytes)?;
        let modulus_limbs = (modulus.bits() / 64) + 1;
        let order_limbs = (order.bits() / 64) + 1;

        let result: Result<Vec<u8>, ApiError> = match (modulus_limbs, order_limbs) {
            (4, 4) => {
                G1ApiImplementation::<U256Repr, U256Repr>::multiexp(&bytes)
            },
            (5, 4) => {
                G1ApiImplementation::<U320Repr, U256Repr>::multiexp(&bytes)
            },
            (6, 4) => {
                G1ApiImplementation::<U320Repr, U256Repr>::multiexp(&bytes)
            },
            (7, 4) => {
                G1ApiImplementation::<U320Repr, U256Repr>::multiexp(&bytes)
            },
            (8, 4) => {
                G1ApiImplementation::<U320Repr, U256Repr>::multiexp(&bytes)
            },
            (9, 4) => {
                G1ApiImplementation::<U320Repr, U256Repr>::multiexp(&bytes)
            },
            (10, 4) => {
                G1ApiImplementation::<U320Repr, U256Repr>::multiexp(&bytes)
            },
            (5, 5) => {
                G1ApiImplementation::<U320Repr, U320Repr>::multiexp(&bytes)
            },
            (field_limbs, group_limbs) => {
                unimplemented!("unimplemented for {} modulus and {} group limbs", field_limbs, group_limbs);
            }
        };

        result
    }
}