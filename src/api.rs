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
use crate::weierstrass::{WeierstrassCurve, CurvePoint, CurveType, Group};
use crate::field::{SizedPrimeField, field_from_modulus};
use crate::field::PrimeFieldElement;
use crate::representation::ElementRepr;

use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};

const BYTES_FOR_LENGTH_ENCODING: usize = 1;

pub trait PrecompileAPI {
    // fn decompress_point(bytes: &[u8]) -> Vec<u8>;
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ()>;
    // fn mul_point(bytes: &[u8]) -> Vec<u8>;
    // fn multiexp(bytes: &[u8]) -> Vec<u8>;
}

pub struct ApiImplementation;

impl PrecompileAPI for ApiImplementation {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ()> {
        let ((modulus, modulus_len), rest) = get_field_params(bytes)?;
        let field = field_from_modulus(modulus).ok_or(())?;
        if rest.len() < modulus_len {
            return Err(());
        }
        let (a_encoding, rest) = rest.split_at(modulus_len);
        let a = PrimeFieldElement::from_be_bytes(&field, a_encoding).map_err(|_| ())?;
        if rest.len() < modulus_len {
            return Err(());
        }
        let (b_encoding, rest) = rest.split_at(modulus_len);
        let b = PrimeFieldElement::from_be_bytes(&field, b_encoding).map_err(|_| ())?;
        let ((order, order_len), rest) = get_curve_params(rest)?;
        let order = BigUint::from_bytes_be(&order);
        let group_field = field_from_modulus(order).ok_or(())?;

        let curve = WeierstrassCurve::new(&group_field, a, b, CurveType::Generic);

        let mut p_0 = decode_point_from_xy(rest, modulus_len, &curve)?;
        let p_1 = decode_point_from_xy(rest, modulus_len, &curve)?;

        p_0.add_assign(&p_1);

        let (x, y) = p_0.into_xy();
        let x = x.into_repr();
        let y = y.into_repr();

        let mut x_bytes: Vec<u8> = vec![];
        x.write_be(&mut x_bytes).map_err(|_| ())?;
        if x_bytes.len() < modulus_len {
            x_bytes.reverse();
            x_bytes.resize(modulus_len, 0u8);
            x_bytes.reverse();
        }

        let mut y_bytes: Vec<u8> = vec![];
        y.write_be(&mut y_bytes).map_err(|_| ())?;
        if y_bytes.len() < modulus_len {
            y_bytes.reverse();
            y_bytes.resize(modulus_len, 0u8);
            y_bytes.reverse();
        }

        let mut result = vec![];
        result.append(&mut x_bytes);
        result.append(&mut y_bytes);

        Ok(result)
    }
}

fn get_field_params(bytes: &[u8]) -> Result<((BigUint, usize), &[u8]), ()> {
    if bytes.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(());
    }

    let (modulus_len, rest) = bytes.split_at(BYTES_FOR_LENGTH_ENCODING);
    let modulus_len = modulus_len[0] as usize;
    if rest.len() < modulus_len {
        return Err(());
    }
    let (modulus_encoding, rest) = rest.split_at(BYTES_FOR_LENGTH_ENCODING);
    let modulus = BigUint::from_bytes_be(&modulus_encoding);
    if modulus.is_zero() {
        return Err(());
    }
    Ok(((modulus, modulus_len), rest))
}

fn get_curve_params(bytes: &[u8]) -> Result<((&[u8], usize), &[u8]), ()> {
    if bytes.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(());
    }

    let (order_len, rest) = bytes.split_at(BYTES_FOR_LENGTH_ENCODING);
    let order_len = order_len[0] as usize;
    if rest.len() < order_len {
        return Err(());
    }
    let (order_encoding, rest) = rest.split_at(BYTES_FOR_LENGTH_ENCODING);

    Ok(((order_encoding, order_len), rest)
}

fn decode_point_from_xy<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>,
    GE: ElementRepr,
    G: SizedPrimeField<Repr = GE>
    >
    (
        bytes: &[u8], 
        field_byte_len: usize,
        curve: &'a WeierstrassCurve<'a, FE, F, GE, G>
    ) -> Result<CurvePoint<'a, FE, F, GE, G>, ()>
{
    if bytes.len() < field_byte_len {
        return Err(());
    }
    let (x_encoding, rest) = bytes.split_at(field_byte_len);
    let x = PrimeFieldElement::from_be_bytes(curve.field, x_encoding).map_err(|_| ())?;
    if rest.len() < field_byte_len {
        return Err(());
    }
    let (y_encoding, rest) = bytes.split_at(field_byte_len);
    let y = PrimeFieldElement::from_be_bytes(curve.field, y_encoding).map_err(|_| ())?;
    
    let p: CurvePoint<'a, FE, F, GE, G> = CurvePoint::point_from_xy(&curve, x, y);
    
    Ok(p)
}