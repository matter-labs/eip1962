mod constants;
mod parsers;
mod utils;
mod g1_constants;

use crate::errors::ApiError;

use crate::public_interface::decode_utils::*;
use crate::public_interface::constants::*;
use self::parsers::*;
use self::g1_constants::*;
use self::utils::*;

pub struct GasMeter;

impl GasMeter {
    pub fn meter(bytes: &[u8]) -> Result<u64, ApiError> {
        let (op_type, rest) = split(bytes, OPERATION_ENCODING_LENGTH , "Input should be longer than operation type encoding")?;
        let result = match op_type[0] {
            OPERATION_G1_ADD => {
                Self::meter_g1_addition(&rest)
            },
            OPERATION_G1_MUL => {
                Self::meter_g1_multiplication(&rest)
            },
            OPERATION_G1_MULTIEXP => {
                Self::meter_g1_multiexp(&rest)
            },
            // OPERATION_G2_ADD => {
            //     PublicG2Api::multiexp(&rest)
            // },
            // OPERATION_G2_MUL => {
            //     PublicG2Api::multiexp(&rest)
            // },
            // OPERATION_G2_MULTIEXP => {
            //     PublicG2Api::multiexp(&rest)
            // },
            // OPERATION_PAIRING => {
            //     PublicPairingApi::pair(&rest)
            // },
            _ => {
                Err(ApiError::InputError("Unknown operation type".to_owned()))
            }
        };

        result
    }

    fn meter_g1_addition(bytes: &[u8]) -> Result<u64, ApiError> {
        use std::u64;
        // apriori formula:
        // CONSTANT + A*LIBMS^1 + B*LIMBS^2 + C*LIMBS^3
        // linear term account for additions and should be small
        // quadratic term accounts for multiplications and should give the largest contribution
        // cubic term reflects inversions

        let (modulus, _, _) = parse_g1_curve_parameters(&bytes)?;
        let num_limbs = num_limbs_for_modulus(&modulus)? as u64;
        let price = evaluate_poly(&vec![
            G1_ADDITION_CONSTANT,
            FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER,
            FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER,
            FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER
        ], num_limbs)?;
        return Ok(price)
    }

    fn meter_g1_multiplication(bytes: &[u8]) -> Result<u64, ApiError> {
        use std::u64;
        // apriori formula:
        // (CONSTANT + A*LIBMS^1 + B*LIMBS^2)*ORDER_WORDS + C*LIMBS^3
        // linear term account for additions and should be small
        // quadratic term accounts for multiplications and should give the largest contribution
        // cubic term reflects inversions

        let (modulus, order, _) = parse_g1_curve_parameters(&bytes)?;
        let num_limbs = num_limbs_for_modulus(&modulus)? as u64;
        let num_order_words = num_units_for_group_order(&order)? as u64;

        let inner_term = evaluate_poly(&vec![
            G1_MULTIPLICATION_CONSTANT,
            FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER,
            FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER,
        ], num_limbs)?;
        let order_words_factor = num_order_words.checked_mul(ORDER_UNITS_LINEAR_TERM).ok_or(ApiError::Overflow)?;
        let inner_term = inner_term.checked_mul(order_words_factor).ok_or(ApiError::Overflow)?;

        let cubic_term = evaluate_poly(&vec![
            0,
            0,
            0,
            FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER
        ], num_limbs)?;

        let price = inner_term.checked_add(cubic_term).ok_or(ApiError::Overflow)?;

        return Ok(price)
    }

    fn meter_g1_multiexp(bytes: &[u8]) -> Result<u64, ApiError> {
        use std::u64;
        // apriori formula:
        // (CONSTANT + A*LIBMS^1 + B*LIMBS^2)*ORDER_WORDS*NUM_PAIR/DISCOUNT + C*LIMBS^3
        // linear term account for additions and should be small
        // quadratic term accounts for multiplications and should give the largest contribution
        // cubic term reflects inversions

        let (modulus, order, rest) = parse_g1_curve_parameters(&bytes)?;
        let (num_pairs_encoding, _) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
        let num_pairs = num_pairs_encoding[0] as u64;
        if num_pairs == 0 {
            return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
        }

        let num_limbs = num_limbs_for_modulus(&modulus)? as u64;
        let num_order_words = num_units_for_group_order(&order)? as u64;

        let inner_term = evaluate_poly(&vec![
            G1_MULTIEXP_CONSTANT,
            FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER,
            FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER,
        ], num_limbs)?;
        let order_words_factor = num_order_words.checked_mul(ORDER_UNITS_LINEAR_TERM).ok_or(ApiError::Overflow)?;
        let num_pair_factor = order_words_factor.checked_mul(num_pairs).ok_or(ApiError::Overflow)?;
        let num_pair_factor = num_pairs / MULTIEXP_DISCOUNT;
        let inner_term = inner_term.checked_mul(num_pair_factor).ok_or(ApiError::Overflow)?;

        let cubic_term = evaluate_poly(&vec![
            0,
            0,
            0,
            FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER
        ], num_limbs)?;

        let price = inner_term.checked_add(cubic_term).ok_or(ApiError::Overflow)?;

        return Ok(price)
    }
}