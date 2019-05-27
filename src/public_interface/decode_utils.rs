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

pub(crate) fn decode_biguint_with_length<
    'a
    >
    (
        bytes: &'a [u8], 
    ) -> Result<(BigUint, &'a [u8]), ()>
{
    if bytes.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(());
    }
    let (length_encoding, rest) = bytes.split_at(BYTES_FOR_LENGTH_ENCODING);
    let length = length_encoding[0] as usize;
    if rest.len() < length {
        return Err(());
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
        &'a [u8]), ()> {
    if bytes.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(());
    }
    let (modulus_len, rest) = bytes.split_at(BYTES_FOR_LENGTH_ENCODING);
    let modulus_len = modulus_len[0] as usize;
    if rest.len() < modulus_len {
        return Err(());
    }
    let (modulus_encoding, rest) = rest.split_at(modulus_len);
    let modulus = BigUint::from_bytes_be(&modulus_encoding);
    if modulus.is_zero() {
        return Err(());
    }
    if rest.len() < modulus_len {
        return Err(());
    }
    let (a_encoding, rest) = rest.split_at(modulus_len);

    if rest.len() < modulus_len {
        return Err(());
    }
    let (b_encoding, rest) = rest.split_at(modulus_len);

    if rest.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(());
    }
    let (order_len, rest) = rest.split_at(BYTES_FOR_LENGTH_ENCODING);
    let order_len = order_len[0] as usize;
    if rest.len() < order_len {
        return Err(());
    }
    let (order_encoding, rest) = rest.split_at(order_len);
    let order = BigUint::from_bytes_be(&order_encoding);
    if order.is_zero() {
        return Err(());
    }
    if rest.len() == 0 {
        return Err(());
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
/// - modulus
pub(crate) fn parse_curve_type_and_modulus<'a>(bytes: &'a [u8]) -> Result<(u8, BigUint), ()> {
    if bytes.len() < CURVE_TYPE_LENGTH {
        return Err(());
    }
    let (curve_type, rest) = bytes.split_at(CURVE_TYPE_LENGTH);
    let curve_type = curve_type[0];
    if bytes.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(());
    }
    let (modulus_len, rest) = rest.split_at(BYTES_FOR_LENGTH_ENCODING);
    let modulus_len = modulus_len[0] as usize;
    if rest.len() < modulus_len {
        return Err(());
    }
    let (modulus_encoding, rest) = rest.split_at(modulus_len);
    let modulus = BigUint::from_bytes_be(&modulus_encoding);
    if modulus.is_zero() {
        return Err(());
    }
    if rest.len() < modulus_len {
        return Err(());
    }

    Ok((curve_type, modulus))
}