use crate::weierstrass::Group;
use crate::weierstrass::twist;
use crate::weierstrass::cubic_twist;
use crate::field::{SizedPrimeField, field_from_modulus};
use crate::fp::Fp;
use crate::extension_towers::fp2;
use crate::extension_towers::fp3;
use crate::representation::ElementRepr;

use num_bigint::BigUint;
use num_traits::{Zero};

use super::constants::*;

use crate::errors::ApiError;

pub(crate) fn decode_biguint_with_length<
    'a
    >
    (
        bytes: &'a [u8], 
    ) -> Result<(BigUint, &'a [u8]), ApiError>
{
    if bytes.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(ApiError::InputError("Input is not long enough to get modulus length".to_owned()));
    }
    let (length_encoding, rest) = bytes.split_at(BYTES_FOR_LENGTH_ENCODING);
    let length = length_encoding[0] as usize;
    if rest.len() < length {
        return Err(ApiError::InputError("Input is not long enough to get modulus".to_owned()));
    }
    let (be_encoding, rest) = rest.split_at(length);
    let x = BigUint::from_bytes_be(&be_encoding);

    Ok((x, rest))
}

/// return:
/// - modulus, 
/// - modulus_len, 
/// - a_bytes, 
/// - b_bytes, 
/// - scalar field modulus 
/// - scalar field length
/// - rest
pub(crate) fn parse_encodings<'a>(bytes: &'a [u8]) -> Result<(
        BigUint, 
        usize,
        &'a [u8],
        &'a [u8],
        BigUint,
        usize,
        &'a [u8]), ApiError> {
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
    if rest.len() < modulus_len {
        return Err(ApiError::InputError("Input is not long enough to get A parameter".to_owned()));
    }
    let (a_encoding, rest) = rest.split_at(modulus_len);

    if rest.len() < modulus_len {
        return Err(ApiError::InputError("Input is not long enough to get B parameter".to_owned()));
    }
    let (b_encoding, rest) = rest.split_at(modulus_len);

    if rest.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(ApiError::InputError("Input is not long enough to get order length".to_owned()));
    }
    let (order_len, rest) = rest.split_at(BYTES_FOR_LENGTH_ENCODING);
    let order_len = order_len[0] as usize;
    if rest.len() < order_len {
        return Err(ApiError::InputError("Input is not long enough to get main group order".to_owned()));
    }
    let (order_encoding, rest) = rest.split_at(order_len);
    let order = BigUint::from_bytes_be(&order_encoding);
    if order.is_zero() {
        return Err(ApiError::UnexpectedZero("Main group order can not be zero".to_owned()));
    }
    if rest.len() == 0 {
        return Err(ApiError::InputError("Input is not long enough".to_owned()));
    }

    Ok(
        (
            modulus,
            modulus_len,
            a_encoding,
            b_encoding,
            order,
            order_len,
            rest
        )
    )
}

/// return:
/// - modulus, 
/// - modulus_len, 
/// - extension degree
/// - a_bytes, 
/// - b_bytes, 
/// - scalar field modulus 
/// - scalar field length
/// - rest
pub(crate) fn parse_encodings_in_extension<'a>(bytes: &'a [u8]) -> Result<(
        BigUint, 
        usize,
        u8,
        &'a [u8],
        &'a [u8],
        BigUint,
        usize,
        &'a [u8]), ApiError> {
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
        return Err(ApiError::InputError("Modulus can not be zero".to_owned()));
    }
    if rest.len() < EXTENSION_DEGREE_ENCODING_LENGTH {
        return Err(ApiError::InputError("Input is not long enough to get extension degree".to_owned()));
    }
    let (extension_degree, rest) = rest.split_at(EXTENSION_DEGREE_ENCODING_LENGTH);
    let extension_degree = extension_degree[0];
    if !(extension_degree == EXTENSION_DEGREE_2 || extension_degree == EXTENSION_DEGREE_3) {
        return Err(ApiError::InputError("Extension degree must be 2 or 3".to_owned()));
    }

    let extension_element_len = (extension_degree as usize) * modulus_len;

    if rest.len() < extension_element_len {
        return Err(ApiError::InputError("Input is not long enough to get A in extension".to_owned()));
    }
    let (a_encoding, rest) = rest.split_at(extension_element_len);

    if rest.len() < extension_element_len {
        return Err(ApiError::InputError("Input is not long enough to get B in extension".to_owned()));
    }
    let (b_encoding, rest) = rest.split_at(extension_element_len);

    if rest.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(ApiError::InputError("Input is not long enough to get group order length".to_owned()));
    }
    let (order_len, rest) = rest.split_at(BYTES_FOR_LENGTH_ENCODING);
    let order_len = order_len[0] as usize;
    if rest.len() < order_len {
        return Err(ApiError::InputError("Input is not long enough to size of the main group".to_owned()));
    }
    let (order_encoding, rest) = rest.split_at(order_len);
    let order = BigUint::from_bytes_be(&order_encoding);
    if order.is_zero() {
        return Err(ApiError::InputError("Main group size can not be zero".to_owned()));
    }
    if rest.len() == 0 {
        return Err(ApiError::InputError("Input is not long enough".to_owned()));
    }

    Ok(
        (
            modulus,
            modulus_len,
            extension_degree,
            a_encoding,
            b_encoding,
            order,
            order_len,
            rest
        )
    )
}

/// return:
/// - modulus
pub(crate) fn parse_curve_type_and_modulus<'a>(bytes: &'a [u8]) -> Result<(u8, BigUint), ApiError> {
    if bytes.len() < CURVE_TYPE_LENGTH {
        return Err(ApiError::InputError("Input is not long enough to get curve type".to_owned()));
    }
    let (curve_type, rest) = bytes.split_at(CURVE_TYPE_LENGTH);
    let curve_type = curve_type[0];
    if bytes.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(ApiError::InputError("Input is not long enough to get modulus length".to_owned()));
    }
    let (modulus_len, rest) = rest.split_at(BYTES_FOR_LENGTH_ENCODING);
    let modulus_len = modulus_len[0] as usize;
    if rest.len() < modulus_len {
        return Err(ApiError::InputError("Input is not long enough to get modulus".to_owned()));
    }
    let (modulus_encoding, rest) = rest.split_at(modulus_len);
    let modulus = BigUint::from_bytes_be(&modulus_encoding);
    if modulus.is_zero() {
        return Err(ApiError::InputError("Modulus can not be zero".to_owned()));
    }
    if rest.len() < modulus_len {
        return Err(ApiError::InputError("Input is not long enough".to_owned()));
    }

    Ok((curve_type, modulus))
}