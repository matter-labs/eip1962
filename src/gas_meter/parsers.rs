use crate::public_interface::decode_utils::*;
use crate::public_interface::decode_g1::*;
use crate::public_interface::constants::*;
use crate::errors::ApiError;
use crate::constants::*;


/// return:
/// - modulus, 
/// - scalar field modulus
/// - rest
/// eats up to the operation-specific parameters
pub(crate) fn parse_g1_curve_parameters<'a>(bytes: &'a [u8]) -> Result<(
    MaxFieldUint, 
    MaxGroupSizeUint,
    &'a [u8]), ApiError> 
{
    let ((modulus, modulus_len), rest) = get_base_field_params(&bytes)?;
    let (_, rest) = split(rest, modulus_len, "Input is not long enough to get A parameter")?;
    let (_, rest) = split(rest, modulus_len, "Input is not long enough to get B parameter")?;

    let (_, order, rest) = parse_group_order_from_encoding(rest)?;

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

/// return:
/// - modulus, 
/// - scalar field modulus
/// - extension degree
/// - rest
/// eats up to the operation-specific parameters
pub(crate) fn parse_g2_curve_parameters<'a>(bytes: &'a [u8]) -> Result<(
    MaxFieldUint, 
    MaxGroupSizeUint,
    u8,
    &'a [u8]), ApiError> 
{
    let ((modulus, modulus_len), rest) = get_base_field_params(&bytes)?;
    let (ext_degree_encoding, rest) = split(&rest, EXTENSION_DEGREE_ENCODING_LENGTH, "Input is not long enough to get extension degree")?;
    let extension_degree = ext_degree_encoding[0];
    if !(extension_degree == EXTENSION_DEGREE_2 || extension_degree == EXTENSION_DEGREE_3) {
        return Err(ApiError::InputError("Invalid extension degree".to_owned()));
    }
    let (_, rest) = split(rest, modulus_len, "Input is not long enough to get non-residue")?;
    let extension_field_element_len = modulus_len * (extension_degree as usize);
    let (_, rest) = split(rest, extension_field_element_len, "Input is not long enough to get A parameter")?;
    let (_, rest) = split(rest, extension_field_element_len, "Input is not long enough to get B parameter")?;

    let (_, order, rest) = parse_group_order_from_encoding(rest)?;
    if rest.len() == 0 {
        return Err(ApiError::InputError("Input is not long enough".to_owned()));
    }

    Ok(
        (
            modulus,
            order,
            extension_degree,
            rest
        )
    )
}

pub(crate) fn parse_mnt_pairing_parameters<'a>(bytes: &'a [u8]) -> Result<(
    MaxFieldUint, 
    MaxGroupSizeUint,
    usize,
    (u64, u64),
    (u64, u64),
    (u64, u64),
    &'a [u8]), ApiError> 
{
    use crate::public_interface::sane_limits::*;

    let ((modulus, modulus_len), rest) = get_base_field_params(&bytes)?;
    let (_, rest) = split(rest, modulus_len, "Input is not long enough to get A parameter")?;
    let (_, rest) = split(rest, modulus_len, "Input is not long enough to get B parameter")?;
    let (_, rest) = split(rest, modulus_len, "Input is not long enough to get non-residue")?;

    let (_, order, rest) = parse_group_order_from_encoding(rest)?;

    let (x, rest) = decode_loop_parameter_scalar_with_bit_limit(&rest, MAX_ATE_PAIRING_ATE_LOOP_COUNT)?;
    if x.is_zero() {
        return Err(ApiError::InputError("Ate pairing loop count parameters can not be zero".to_owned()));
    }

    let ate_loop_bits = x.bits();
    let ate_loop_hamming = calculate_hamming_weight(&x.as_ref());

    if ate_loop_hamming > MAX_ATE_PAIRING_ATE_LOOP_COUNT_HAMMING {
        return Err(ApiError::InputError("X has too large hamming weight".to_owned()));
    }

    let (x_sign, rest) = split(rest, SIGN_ENCODING_LENGTH, "Input is not long enough to get X sign encoding")?;
    let _ = match x_sign[0] {
        SIGN_PLUS => false,
        SIGN_MINUS => true,
        _ => {
            return Err(ApiError::InputError("X sign is not encoded properly".to_owned()));
        },
    };

    let (exp_w0, rest) = decode_loop_parameter_scalar_with_bit_limit(&rest, MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH)?;
    if exp_w0.is_zero() {
        return Err(ApiError::InputError("Final exp w0 loop count parameters can not be zero".to_owned()));
    }
    let exp_w0_bits = exp_w0.bits();
    let exp_w0_hamming = calculate_hamming_weight(&exp_w0.as_ref());

    let (exp_w1, rest) = decode_loop_parameter_scalar_with_bit_limit(&rest, MAX_ATE_PAIRING_FINAL_EXP_W1_BIT_LENGTH)?;
    if exp_w1.is_zero() {
        return Err(ApiError::InputError("Final exp w1 loop count parameters can not be zero".to_owned()));
    }
    let exp_w1_bits = exp_w1.bits();
    let exp_w1_hamming = calculate_hamming_weight(&exp_w1.as_ref());

    let (exp_w0_sign, rest) = split(rest, SIGN_ENCODING_LENGTH, "Input is not long enough to get exp_w0 sign encoding")?;
    let _ = match exp_w0_sign[0] {
        SIGN_PLUS => false,
        SIGN_MINUS => true,
        _ => {
            return Err(ApiError::InputError("Exp_w0 sign is not encoded properly".to_owned()));
        },
    };

    let (num_pairs_encoding, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
    let num_pairs = num_pairs_encoding[0] as usize;

    if num_pairs == 0 {
        return Err(ApiError::InputError("Zero pairs encoded".to_owned()));
    }
    
    if rest.len() == 0 {
        return Err(ApiError::InputError("Input is not long enough".to_owned()));
    }

    Ok(
        (
            modulus,
            order,
            num_pairs,
            (ate_loop_bits as u64, ate_loop_hamming as u64),
            (exp_w0_bits as u64, exp_w0_hamming as u64),
            (exp_w1_bits as u64, exp_w1_hamming as u64),
            rest
        )
    )
}