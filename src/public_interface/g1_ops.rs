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
use crate::weierstrass::curve::{WeierstrassCurve, CurvePoint};
use crate::field::{SizedPrimeField, field_from_modulus};
use crate::fp::Fp;
use crate::field::{U256Repr, U320Repr};
use crate::representation::ElementRepr;
use crate::multiexp::peppinger;

use num_bigint::BigUint;
use num_traits::{Zero};

#[macro_use]
use super::decode_g1::*;

pub trait G1Api {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ()>;
    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ()>;
    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ()>;
}

pub struct G1ApiImplementation<FE: ElementRepr, GE: ElementRepr> {
    _marker_fe: std::marker::PhantomData<FE>,
    _marker_ge: std::marker::PhantomData<GE>
}

impl<FE: ElementRepr, GE: ElementRepr> G1Api for G1ApiImplementation<FE, GE> {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ()> {
        let (field, modulus_len, rest) = create_base_field!(bytes, FE);
        let (a, b, rest) = get_ab_in_base_field!(rest, field, modulus_len);
        let (group, order_len, rest) = create_group!(rest, GE);

        let curve = WeierstrassCurve::new(&group, a, b);

        let (mut p_0, rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;
        let (p_1, _rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;

        p_0.add_assign(&p_1);

        serialize_g1_point(modulus_len, &p_0)   
    }

    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ()> {
        let (field, modulus_len, rest) = create_base_field!(bytes, FE);
        let (a, b, rest) = get_ab_in_base_field!(rest, field, modulus_len);
        let (group, order_len, rest) = create_group!(rest, GE);

        let curve = WeierstrassCurve::new(&group, a, b);

        let (p_0, rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;
        let (scalar, _rest) = decode_scalar_representation(rest, order_len, &group)?;

        let p = p_0.mul(&scalar);

        serialize_g1_point(modulus_len, &p)   
    }

    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ()> {
        let (field, modulus_len, rest) = create_base_field!(bytes, FE);
        let (a, b, rest) = get_ab_in_base_field!(rest, field, modulus_len);
        let (group, order_len, rest) = create_group!(rest, GE);

        let curve = WeierstrassCurve::new(&group, a, b);

        let expected_pair_len = 2*modulus_len + order_len;
        if rest.len() % expected_pair_len != 0 {
            return Err(());
        }

        let expected_pairs = rest.len() / expected_pair_len;
        if expected_pairs == 0 {
            return Err(());
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