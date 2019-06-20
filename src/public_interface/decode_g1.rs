use crate::weierstrass::Group;
use crate::weierstrass::curve::{WeierstrassCurve, CurvePoint};
use crate::field::{SizedPrimeField, field_from_modulus};
use crate::fp::Fp;
use crate::representation::ElementRepr;

use super::constants::*;
use super::decode_fp::*;

use num_bigint::BigUint;
use num_traits::{Zero};

use crate::errors::ApiError;

macro_rules! create_base_field {
    ($bytes:expr, $repr:tt) => {
        {
            let ((modulus, modulus_len), rest) = get_base_field_params($bytes)?;
            let field = field_from_modulus::<$repr>(modulus).map_err(|_| {
                ApiError::InputError("Failed to create prime field from modulus".to_owned())
            })?;
            if rest.len() < modulus_len {
                return Err(ApiError::InputError("Input is not long enough".to_owned()));
            }

            (field, modulus_len, rest)
        }
    }
}

macro_rules! create_base_field_with_modulus {
    ($bytes:expr, $repr:tt) => {
        {
            let ((modulus, modulus_len), rest) = get_base_field_params($bytes)?;
            let field = field_from_modulus::<$repr>(modulus.clone()).map_err(|_| {
                ApiError::InputError("Failed to create prime field from modulus".to_owned())
            })?;
            if rest.len() < modulus_len {
                return Err(ApiError::InputError("Input is not long enough".to_owned()));
            }

            (field, modulus_len, modulus, rest)
        }
    }
}

macro_rules! get_ab_in_base_field {
    ($rest:expr, $field:expr, $modulus_len: expr) => {
        {
            if $rest.len() < $modulus_len {
                return Err(ApiError::InputError("Input is not long enough to A parameter".to_owned()));
            }
            let (a_encoding, rest) = $rest.split_at($modulus_len);
            let a = Fp::from_be_bytes(&$field, a_encoding, true).map_err(|_| {
                ApiError::InputError("Failed to parse A in Fp".to_owned())
            })?;
            if rest.len() < $modulus_len {
                return Err(ApiError::InputError("Input is not long enough to B paramter".to_owned()));
            }
            let (b_encoding, rest) = rest.split_at($modulus_len);
            let b = Fp::from_be_bytes(&$field, b_encoding, true).map_err(|_| {
                ApiError::InputError("Failed to parse B in Fp".to_owned())
            })?;
            
            (a, b, rest)
        }
    }
}

macro_rules! create_group {
    ($bytes:expr, $repr:tt) => {
        {
            let ((order, order_len), rest) = get_g1_curve_params($bytes)?;
            let order = BigUint::from_bytes_be(&order);
            let group = field_from_modulus::<$repr>(order).map_err(|_| {
                ApiError::InputError("Failed to create scalar field from group order".to_owned())
            })?;

            (group, order_len, rest)
        }
    }
}

pub(crate) fn serialize_g1_point<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>,
    GE: ElementRepr,
    G: SizedPrimeField<Repr = GE>
    >
    (
        modulus_len: usize,
        point: &CurvePoint<'a, FE, F, GE, G>
    ) -> Result<Vec<u8>, ApiError>
{
    let (x, y) = point.into_xy();
    let mut result = serialize_fp_fixed_len(modulus_len, &x)?;
    result.extend(serialize_fp_fixed_len(modulus_len, &y)?);

    Ok(result)
}

pub(crate) fn get_base_field_params(bytes: &[u8]) -> Result<((BigUint, usize), &[u8]), ApiError> {
    if bytes.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(ApiError::InputError("Input is not long enough to get modulus length".to_owned()));
    }
    let (modulus_len, rest) = bytes.split_at(BYTES_FOR_LENGTH_ENCODING);
    let modulus_len = modulus_len[0] as usize;

    if rest.len() < modulus_len {
        return Err(ApiError::InputError("Input is not long enough to get modulus".to_owned()));
    }
    let (modulus_encoding, rest) = rest.split_at(modulus_len);
    let modulus = BigUint::from_bytes_be(&modulus_encoding);
    if modulus.is_zero() {
        return Err(ApiError::UnexpectedZero("Modulus can not be zero".to_owned()));
    }

    Ok(((modulus, modulus_len), rest))
}

pub(crate) fn get_g1_curve_params(bytes: &[u8]) -> Result<((&[u8], usize), &[u8]), ApiError> {
    if bytes.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(ApiError::InputError("Input is not long enough to get group size length".to_owned()));
    }

    let (order_len, rest) = bytes.split_at(BYTES_FOR_LENGTH_ENCODING);
    let order_len = order_len[0] as usize;
    if rest.len() < order_len {
        return Err(ApiError::InputError("Input is not long enough to get main group order size".to_owned()));
    }
    let (order_encoding, rest) = rest.split_at(order_len);

    Ok(((order_encoding, order_len), rest))
}

pub(crate) fn decode_g1_point_from_xy<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>,
    GE: ElementRepr,
    G: SizedPrimeField<Repr = GE>
    >
    (
        bytes: &'a [u8], 
        field_byte_len: usize,
        curve: &'a WeierstrassCurve<'a, FE, F, GE, G>
    ) -> Result<(CurvePoint<'a, FE, F, GE, G>, &'a [u8]), ApiError>
{
    if bytes.len() < field_byte_len {
        return Err(ApiError::InputError("Input is not long enough to get X".to_owned()));
    }
    let (x_encoding, rest) = bytes.split_at(field_byte_len);
    let x = Fp::from_be_bytes(curve.base_field, x_encoding, true).map_err(|_| {
        ApiError::InputError("Failed to parse X".to_owned())
    })?;
    if rest.len() < field_byte_len {
        return Err(ApiError::InputError("Input is not long enough to get Y".to_owned()));
    }
    let (y_encoding, rest) = rest.split_at(field_byte_len);
    let y = Fp::from_be_bytes(curve.base_field, y_encoding, true).map_err(|_| {
        ApiError::InputError("Failed to parse Y".to_owned())
    })?;
    
    let p: CurvePoint<'a, FE, F, GE, G> = CurvePoint::point_from_xy(&curve, x, y);
    
    Ok((p, rest))
}

pub(crate) fn decode_scalar_representation<
    'a,
    GE: ElementRepr,
    G: SizedPrimeField<Repr = GE>
    >
    (
        bytes: &'a [u8], 
        order_byte_len: usize,
        group: &G,
    ) -> Result<(GE, &'a [u8]), ApiError>
{
    if bytes.len() < order_byte_len {
        return Err(ApiError::InputError("Input is not long enough to get scalar".to_owned()));
    }
    let (encoding, rest) = bytes.split_at(order_byte_len);
    let mut repr = GE::default();
    if encoding.len() >= repr.as_ref().len() * 8 {
        repr.read_be(encoding).map_err(|_| {
            ApiError::InputError("Failed to parse scalar".to_owned())
        })?;
    } else {
        let mut padded = vec![0u8; repr.as_ref().len() * 8 - encoding.len()];
        padded.extend_from_slice(encoding);
        repr.read_be(&padded[..]).map_err(|_| {
            ApiError::InputError("Failed to parse scalar".to_owned())
        })?;
    }

    Ok((repr, rest))
}

