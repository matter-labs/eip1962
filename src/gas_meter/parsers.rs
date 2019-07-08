use num_bigint::BigUint;
use num_traits::{Zero};

use crate::public_interface::decode_utils::*;
use crate::public_interface::constants::*;
use crate::errors::ApiError;


/// return:
/// - modulus, 
/// - scalar field modulus
/// - rest
/// eats up to the operation-specific parameters
pub(crate) fn parse_g1_curve_parameters<'a>(bytes: &'a [u8]) -> Result<(
        BigUint, 
        BigUint,
        &'a [u8]), ApiError> {
    let ((modulus, modulus_len), rest) = get_base_field_params(&bytes)?;
    let (_, rest) = split(rest, modulus_len, "Input is not long enough to get A parameter")?;
    let (_, rest) = split(rest, modulus_len, "Input is not long enough to get B parameter")?;

    let (order_len, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get order length")?;
    let order_len = order_len[0] as usize;
    if order_len == 0 {
        return Err(ApiError::UnexpectedZero("Order encoding length is 0".to_owned()));
    }
    let (order_encoding, rest) = split(rest, order_len, "Input is not long enough to get main group order")?;
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
            order,
            rest
        )
    )
}