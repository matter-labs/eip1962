use num_bigint::BigUint;
use num_traits::{Zero};
use num_integer::Integer;

use super::constants::*;

use crate::errors::ApiError;

pub(crate) fn split<'a>(bytes: &'a [u8], at: usize, err: &'static str) 
    -> Result<(&'a [u8], &'a [u8]), ApiError> 
{
    if bytes.len() < at {
        Err(ApiError::InputError(err.to_owned()))
    } else {
        Ok(bytes.split_at(at))
    }
}

pub(crate) fn decode_biguint_with_length<
    'a
    >
    (
        bytes: &'a [u8], 
    ) -> Result<(BigUint, &'a [u8]), ApiError>
{
    let (length_encoding, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
    let length = length_encoding[0] as usize;
    let (be_encoding, rest) = split(rest, length, "Input is not long enough to get modulus")?;
    let x = BigUint::from_bytes_be(&be_encoding);

    Ok((x, rest))
}

pub(crate) fn parse_modulus_and_length<
    'a
    >
    (
        bytes: &'a [u8], 
    ) -> Result<(usize, BigUint, &'a [u8]), ApiError>
{
    let (length_encoding, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
    let length = length_encoding[0] as usize;
    let (be_encoding, rest) = split(rest, length, "Input is not long enough to get modulus")?;
    let x = BigUint::from_bytes_be(&be_encoding);

    Ok((length, x, rest))
}

/// return:
/// - modulus, 
/// - modulus_len, 
/// - extension degree
/// - non-residue encoding
/// - rest
pub(crate) fn parse_modulus_and_extension_degree<'a>(bytes: &'a [u8]) -> Result<(
        BigUint, 
        usize,
        u8,
        &'a [u8],
        &'a [u8]), ApiError> {
    let ((modulus, modulus_len), rest) = get_base_field_params(&bytes)?;
    let (extension_degree, rest) = split(rest, EXTENSION_DEGREE_ENCODING_LENGTH, "Input is not long enough to get extension degree")?;
    let extension_degree = extension_degree[0];
    if !(extension_degree == EXTENSION_DEGREE_2 || extension_degree == EXTENSION_DEGREE_3) {
        return Err(ApiError::InputError("Extension degree must be 2 or 3".to_owned()));
    }

    let (nonresidue_encoding, rest) = split(rest, modulus_len, "Input is not long enough to Fp non-residue")?;
    if rest.len() == 0 {
        return Err(ApiError::InputError("Input is not long enough".to_owned()));
    }

    Ok(
        (
            modulus,
            modulus_len,
            extension_degree,
            nonresidue_encoding,
            rest
        )
    )
}

pub(crate) fn get_base_field_params(bytes: &[u8]) -> Result<((BigUint, usize), &[u8]), ApiError> {
    use crate::constants::THREE_BIGUINT;
    let (modulus_len, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
    let modulus_len = modulus_len[0] as usize;
    if modulus_len == 0 {
        return Err(ApiError::InputError(format!("Modulus is length is zero, file {}, line {}", file!(), line!())));
    }
    let (modulus_encoding, rest) = split(rest, modulus_len, "Input is not long enough to get modulus")?;
    if modulus_encoding[0] == 0u8 {
        return Err(ApiError::InputError(format!("In modulus encoding highest byte is zero, file {}, line {}", file!(), line!())));
    }
    let modulus = BigUint::from_bytes_be(&modulus_encoding);
    if modulus.is_zero() {
        return Err(ApiError::UnexpectedZero("Modulus can not be zero".to_owned()));
    }
    if modulus.is_even() {
        return Err(ApiError::InputError(format!("Modulus is even, file {}, line {}", file!(), line!())));
    }
    if modulus < *THREE_BIGUINT {
        return Err(ApiError::InputError(format!("Modulus is less than 3, file {}, line {}", file!(), line!())));
    }

    Ok(((modulus, modulus_len), rest))
}

pub(crate) fn num_limbs_for_modulus(modulus: &BigUint) -> Result<usize, ApiError> {
    use crate::field::calculate_num_limbs;

    let modulus_limbs = calculate_num_limbs(&modulus).map_err(|_| ApiError::InputError(format!("Modulus is too large, file {}, line {}", file!(), line!())) )?;

    Ok(modulus_limbs)
}

pub(crate) fn num_units_for_group_order(order: &BigUint) -> Result<usize, ApiError> {
    let limbs = (order.bits() + 63) / 64;
    if limbs > 16 {
        return Err(ApiError::InputError(format!("Group order is too large, file {}, line {}", file!(), line!())));
    }

    Ok(limbs)
}

pub(crate) fn decode_loop_parameter_scalar_with_bit_limit<
    'a
    >
    (
        bytes: &'a [u8], 
        bit_limit: usize,
    ) -> Result<(BigUint, &'a [u8]), ApiError>
{
    let (length_encoding, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
    let max_length_for_bits = (bit_limit + 7) / 8;
    let length = length_encoding[0] as usize;
    if length == 0 {
        return Err(ApiError::InputError(format!("Loop parameter scalar has zero length, file {}, line {}", file!(), line!())));
    }
    if length > max_length_for_bits {
        return Err(ApiError::InputError(format!("Scalar is too large for bit length, file {}, line {}", file!(), line!())));
    }
    let (be_encoding, rest) = split(rest, length, "Input is not long enough to get modulus")?;
    let x = BigUint::from_bytes_be(&be_encoding);
    let num_bits = x.bits();
    if num_bits > bit_limit {
        return Err(ApiError::InputError(format!("Number of bits for scalar is too large, file {}, line {}", file!(), line!())));
    }

    Ok((x, rest))
}

pub(crate) fn calculate_hamming_weight(representation: &[u64]) -> u32 {
    let mut weight = 0;
    for el in representation.iter() {
        weight += el.count_ones();
    }

    weight
}