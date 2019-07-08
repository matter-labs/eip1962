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
/// - a_bytes, 
/// - b_bytes, 
/// - scalar field modulus 
/// - scalar field length
/// - rest
// pub(crate) fn parse_encodings<'a>(bytes: &'a [u8]) -> Result<(
//         BigUint, 
//         usize,
//         &'a [u8],
//         &'a [u8],
//         BigUint,
//         usize,
//         &'a [u8]), ApiError> {
//     let ((modulus, modulus_len), rest) = get_base_field_params(&bytes)?;
//     // let (modulus_len, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
//     // let modulus_len = modulus_len[0] as usize;
//     // let (modulus_encoding, rest) = split(rest, modulus_len, "Input is not long enough to get modulus")?;
//     // let modulus = BigUint::from_bytes_be(&modulus_encoding);
//     // if modulus.is_zero() {
//     //     return Err(ApiError::UnexpectedZero("Modulus can not be zero".to_owned()));
//     // }
//     let (a_encoding, rest) = split(rest, modulus_len, "Input is not long enough to get A parameter")?;
//     let (b_encoding, rest) = split(rest, modulus_len, "Input is not long enough to get B parameter")?;

//     let (order_len, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get order length")?;
//     let order_len = order_len[0] as usize;
//     let (order_encoding, rest) = split(rest, order_len, "Input is not long enough to get main group order")?;
//     let order = BigUint::from_bytes_be(&order_encoding);
//     if order.is_zero() {
//         return Err(ApiError::UnexpectedZero("Main group order can not be zero".to_owned()));
//     }
//     if rest.len() == 0 {
//         return Err(ApiError::InputError("Input is not long enough".to_owned()));
//     }

//     Ok(
//         (
//             modulus,
//             modulus_len,
//             a_encoding,
//             b_encoding,
//             order,
//             order_len,
//             rest
//         )
//     )
// }

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
    // let (modulus_len, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
    // let modulus_len = modulus_len[0] as usize;
    // let (modulus_encoding, rest) = split(rest, modulus_len, "Input is not long enough to get modulus")?;
    // let modulus = BigUint::from_bytes_be(&modulus_encoding);
    // if modulus.is_zero() {
    //     return Err(ApiError::InputError("Modulus can not be zero".to_owned()));
    // }
    // if modulus.is_even() {
    //     return Err(ApiError::InputError(format!("Modulus is even, file {}, line {}", file!(), line!())));
    // }
    // if modulus < *THREE_BIGUINT {
    //     return Err(ApiError::InputError(format!("Modulus is less than 3, file {}, line {}", file!(), line!())));
    // }
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

/// return:
/// - modulus, 
/// - modulus_len, 
/// - extension degree
/// - non-residue
/// - a_bytes, 
/// - b_bytes, 
/// - scalar field modulus 
/// - scalar field length
/// - rest
// pub(crate) fn parse_encodings_in_extension<'a>(bytes: &'a [u8]) -> Result<(
//         BigUint, 
//         usize,
//         u8,
//         &'a [u8],
//         &'a [u8],
//         &'a [u8],
//         BigUint,
//         usize,
//         &'a [u8]), ApiError> {
//     let ((modulus, modulus_len), rest) = get_base_field_params(&bytes)?;
//     // let (modulus_len, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
//     // let modulus_len = modulus_len[0] as usize;
//     // let (modulus_encoding, rest) = split(rest, modulus_len, "Input is not long enough to get modulus")?;
//     // let modulus = BigUint::from_bytes_be(&modulus_encoding);
//     // if modulus.is_zero() {
//     //     return Err(ApiError::InputError("Modulus can not be zero".to_owned()));
//     // }
//     let (extension_degree, rest) = split(rest, EXTENSION_DEGREE_ENCODING_LENGTH, "Input is not long enough to get extension degree")?;
//     let extension_degree = extension_degree[0];
//     if !(extension_degree == EXTENSION_DEGREE_2 || extension_degree == EXTENSION_DEGREE_3) {
//         return Err(ApiError::InputError("Extension degree must be 2 or 3".to_owned()));
//     }

//     let (nonresidue_encoding, rest) = split(rest, modulus_len, "Input is not long enough to Fp non-residue")?;

//     let extension_element_len = (extension_degree as usize) * modulus_len;
//     let (a_encoding, rest) = split(rest, extension_element_len, "Input is not long enough to get A in extension")?;
//     let (b_encoding, rest) = split(rest, extension_element_len, "Input is not long enough to get B in extension")?;

//     let (order_len, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get group order length")?;
//     let order_len = order_len[0] as usize;
//     let (order_encoding, rest) = split(rest, order_len, "Input is not long enough to size of the main group")?;
//     let order = BigUint::from_bytes_be(&order_encoding);
//     if order.is_zero() {
//         return Err(ApiError::InputError(format!("Main group size can not be zero, file {}, line {}", file!(), line!())));
//     }
//     if rest.len() == 0 {
//         return Err(ApiError::InputError("Input is not long enough".to_owned()));
//     }

//     Ok(
//         (
//             modulus,
//             modulus_len,
//             extension_degree,
//             nonresidue_encoding,
//             a_encoding,
//             b_encoding,
//             order,
//             order_len,
//             rest
//         )
//     )
// }

/// return:
/// - modulus
// pub(crate) fn parse_curve_type_and_modulus<'a>(bytes: &'a [u8]) -> Result<(u8, BigUint), ApiError> {
//     let (curve_type, rest) = split(bytes, CURVE_TYPE_LENGTH, "Input is not long enough to get curve type")?;
//     let curve_type = curve_type[0];
//     let (modulus_len, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
//     let modulus_len = modulus_len[0] as usize;
//     let (modulus_encoding, rest) = split(rest, modulus_len, "Input is not long enough to get modulus")?;
//     let modulus = BigUint::from_bytes_be(&modulus_encoding);
//     if modulus.is_zero() {
//         return Err(ApiError::InputError("Modulus can not be zero".to_owned()));
//     }
//     if rest.len() < modulus_len {
//         return Err(ApiError::InputError("Input is not long enough".to_owned()));
//     }

//     Ok((curve_type, modulus))
// }

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

pub(crate) fn num_libs_for_modulus(modulus: &BigUint) -> Result<usize, ApiError> {
    let mut modulus_limbs = (modulus.bits() / 64) + 1;
    if modulus_limbs > 16 {
        return Err(ApiError::InputError(format!("Modulus is too large, file {}, line {}", file!(), line!())));
    }
    if modulus_limbs < 4 {
        modulus_limbs = 4;
    }

    Ok(modulus_limbs)
}