use crate::integers::{MaxFieldUint, MaxGroupSizeUint, MaxLoopParametersUint};

use crate::public_interface::constants::*;

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

fn is_even(x: &MaxFieldUint) -> bool {
    x.low_u64() & 1 == 0
}

pub(crate) fn decode_group_order_with_length<
    'a
    >
    (
        bytes: &'a [u8], 
    ) -> Result<((usize, MaxGroupSizeUint), &'a [u8]), ApiError>
{
    use crate::public_interface::constants::*;

    let (length_encoding, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
    let length = length_encoding[0] as usize;
    if length == 0 {
        return Err(ApiError::InputError(format!("Encoded group length is zero, file {}, line {}", file!(), line!())));
    }
    if length > MAX_GROUP_BYTE_LEN {
        return Err(ApiError::InputError(format!("Encoded group length is too large, file {}, line {}", file!(), line!())));
    }
    let (be_encoding, rest) = split(rest, length, "Input is not long enough to get modulus")?;
    // let first_byte = be_encoding[0];
    // if first_byte == 0 {
    //     return Err(ApiError::InputError(format!("Encoded group length has zero top byte, file {}, line {}", file!(), line!())));
    // }
    let x = MaxGroupSizeUint::from_big_endian(&be_encoding);

    Ok( ((length, x), rest) )
}

pub(crate) fn decode_sign_is_negative<
    'a
    >
    (
        bytes: &'a [u8], 
    ) -> Result<(bool, &'a [u8]), ApiError>
{
    let (x_sign, rest) = split(bytes, SIGN_ENCODING_LENGTH, "Input is not long enough to get sign encoding")?;
    let x_is_negative = match x_sign[0] {
        SIGN_PLUS => false,
        SIGN_MINUS => true,
        _ => {
            return Err(ApiError::InputError("sign is not encoded properly".to_owned()));
        },
    };

    Ok((x_is_negative, rest))
}

use crate::pairings::TwistType;

pub(crate) fn decode_twist_type<
    'a
    >
    (
        bytes: &'a [u8], 
    ) -> Result<(TwistType, &'a [u8]), ApiError>
{
    let (twist_type_encoding, rest) = split(bytes, TWIST_TYPE_LENGTH, "Input is not long enough to get twist type")?;

    let twist_type = match twist_type_encoding[0] {
        TWIST_TYPE_D => TwistType::D,
        TWIST_TYPE_M => TwistType::M, 
        _ => {
            return Err(ApiError::UnknownParameter("Unknown twist type supplied".to_owned()));
        },
    };

    Ok((twist_type, rest))
}

pub(crate) fn decode_boolean<
    'a
    >
    (
        bytes: &'a [u8], 
    ) -> Result<(bool, &'a [u8]), ApiError>
{
    let (boolean_encoding, rest) = split(bytes, BOOLEAN_ENCODING_LENGTH, "Input is not long enough to get boolean")?;
    let boolean = match boolean_encoding[0] {
        BOOLEAN_FALSE => false,
        BOOLEAN_TRUE => true,
        _ => {
            return Err(ApiError::InputError("boolean is not encoded properly".to_owned()));
        },
    };

    Ok((boolean, rest))
}

pub(crate) fn parse_modulus_and_length<
    'a
    >
    (
        bytes: &'a [u8], 
    ) -> Result<(usize, MaxFieldUint, &'a [u8]), ApiError>
{
    use crate::public_interface::constants::*;

    let (length_encoding, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
    let length = length_encoding[0] as usize;
    if length > MAX_MODULUS_BYTE_LEN {
        return Err(ApiError::InputError(format!("Encoded modulus length is too large, file {}, line {}", file!(), line!())));
    }
    let (be_encoding, rest) = split(rest, length, "Input is not long enough to get modulus")?;
    let x = MaxFieldUint::from_big_endian(&be_encoding);

    Ok((length, x, rest))
}

/// return:
/// - modulus, 
/// - modulus_len, 
/// - extension degree
/// - non-residue encoding
/// - rest
pub(crate) fn parse_modulus_and_extension_degree<'a>(bytes: &'a [u8]) -> Result<(
        MaxFieldUint, 
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

pub(crate) fn get_base_field_params(bytes: &[u8]) -> Result<((MaxFieldUint, usize), &[u8]), ApiError> {
    use crate::integers::*;

    let (modulus_len, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
    let modulus_len = modulus_len[0] as usize;
    if modulus_len == 0 {
        return Err(ApiError::InputError(format!("Modulus is length is zero, file {}, line {}", file!(), line!())));
    }
    if modulus_len > MAX_MODULUS_BYTE_LEN {
        return Err(ApiError::InputError(format!("Encoded modulus length is too large, file {}, line {}", file!(), line!())));
    }
    let (modulus_encoding, rest) = split(rest, modulus_len, "Input is not long enough to get modulus")?;
    if modulus_encoding[0] == 0u8 {
        return Err(ApiError::InputError(format!("In modulus encoding highest byte is zero, file {}, line {}", file!(), line!())));
    }
    let modulus = MaxFieldUint::from_big_endian(&modulus_encoding);
    if modulus.is_zero() {
        return Err(ApiError::UnexpectedZero("Modulus can not be zero".to_owned()));
    }
    if is_even(&modulus) {
        return Err(ApiError::InputError(format!("Modulus is even, file {}, line {}", file!(), line!())));
    }
    if modulus < MaxFieldUint::from(3u64) {
        return Err(ApiError::InputError(format!("Modulus is less than 3, file {}, line {}", file!(), line!())));
    }

    Ok(((modulus, modulus_len), rest))
}

pub(crate) fn num_limbs_for_modulus(modulus: &MaxFieldUint) -> Result<usize, ApiError> {
    use crate::field::calculate_num_limbs;

    let modulus_limbs = calculate_num_limbs(modulus.bits())
        .map_err(|_| ApiError::InputError(format!("Modulus is too large, file {}, line {}", file!(), line!())) )?;

    Ok(modulus_limbs)
}

// pub(crate) fn num_units_for_group_order(order: &MaxGroupSizeUint) -> Result<usize, ApiError> {
//     use crate::public_interface::constants::*;

//     let limbs = (order.bits() + 63) / 64;

//     if limbs < NUM_GROUP_LIMBS_MIN {
//         return Err(ApiError::InputError(format!("Group has zero limbs, file {}, line {}", file!(), line!())));
//     }

//     if limbs > NUM_GROUP_LIMBS_MAX {
//         return Err(ApiError::InputError(format!("Group order has too many limbs, file {}, line {}", file!(), line!())));
//     }

//     Ok(limbs)
// }

pub(crate) fn num_units_for_group_order_length(order_len: usize) -> Result<usize, ApiError> {
    use crate::public_interface::constants::*;

    let limbs = (order_len + 7) / 8;

    if limbs < NUM_GROUP_LIMBS_MIN {
        return Err(ApiError::InputError(format!("Group has zero limbs, file {}, line {}", file!(), line!())));
    }

    if limbs > NUM_GROUP_LIMBS_MAX {
        return Err(ApiError::InputError(format!("Group order has too many limbs, file {}, line {}", file!(), line!())));
    }

    Ok(limbs)
}

pub(crate) fn decode_loop_parameter_scalar_with_bit_limit<
    'a
    >
    (
        bytes: &'a [u8], 
        bit_limit: usize,
    ) -> Result<(MaxLoopParametersUint, &'a [u8]), ApiError>
{
    let (length_encoding, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
    let max_length_for_bits = (bit_limit + 7) / 8;
    let length = length_encoding[0] as usize;
    if length == 0 {
        return Err(ApiError::InputError(format!("Loop parameter scalar has zero length, file {}, line {}", file!(), line!())));
    }
    if length > max_length_for_bits {
        return Err(ApiError::InputError(format!("Loop parameter is too large for bit length, max {} bits, got {} bytes, file {}, line {}", bit_limit, length, file!(), line!())));
    }
    let (be_encoding, rest) = split(rest, length, "Input is not long enough to get modulus")?;
    let first_byte = be_encoding[0];
    if first_byte == 0 {
        return Err(ApiError::InputError(format!("Encoded loop parameter has zero top byte, file {}, line {}", file!(), line!())));
    }
    let x = MaxLoopParametersUint::from_big_endian(&be_encoding);
    let num_bits = x.bits();
    if num_bits > bit_limit {
        return Err(ApiError::InputError(format!("Number of bits for scalar is too large, file {}, line {}", file!(), line!())));
    }

    Ok((x, rest))
}

