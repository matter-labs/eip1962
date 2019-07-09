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
        use std::f64;
        // apriori formula:
        // CONSTANT + A*LIBMS^1 + B*LIMBS^2 + C*LIMBS^3
        // linear term account for additions and should be small
        // quadratic term accounts for multiplications and should give the largest contribution
        // cubic term reflects inversions

        let (modulus, _, _) = parse_g1_curve_parameters(&bytes)?;
        let num_limbs = num_limbs_for_modulus(&modulus)? as f64;
        let price = evaluate_poly(&vec![
            G1_ADDITION_CONSTANT,
            FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER,
            FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER,
            FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER
        ], num_limbs)?;
        let price = price as u64;
        return Ok(price)
    }

    fn meter_g1_multiplication(bytes: &[u8]) -> Result<u64, ApiError> {
        use std::f64;
        // apriori formula:
        // CONSTANT + (A*LIBMS^1 + B*LIMBS^2)*ORDER_WORDS + C*LIMBS^3
        // linear term account for additions and should be small
        // quadratic term accounts for multiplications and should give the largest contribution
        // cubic term reflects inversions

        let (modulus, order, _) = parse_g1_curve_parameters(&bytes)?;
        let num_limbs = num_limbs_for_modulus(&modulus)? as f64;
        let num_order_words = num_units_for_group_order(&order)? as f64;

        let inner_term = evaluate_poly(&vec![
            0f64,
            FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER,
            FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER,
        ], num_limbs)?;
        let order_words_factor = num_order_words* ORDER_UNITS_LINEAR_TERM;
        let inner_term = inner_term * order_words_factor;

        let cubic_term = evaluate_poly(&vec![
            0f64,
            0f64,
            0f64,
            FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER
        ], num_limbs)?;

        let price = G1_MULTIPLICATION_CONSTANT;
        let price = price + inner_term;
        let price = price + cubic_term;
        let price = price as u64;

        return Ok(price)
    }

    fn meter_g1_multiexp(bytes: &[u8]) -> Result<u64, ApiError> {
        use std::f64;
        // apriori formula:
        // CONSTANT + (A*LIBMS^1 + B*LIMBS^2)*ORDER_WORDS*NUM_PAIR*DISCOUNT + C*LIMBS^3
        // linear term account for additions and should be small
        // quadratic term accounts for multiplications and should give the largest contribution
        // cubic term reflects inversions

        let (modulus, order, rest) = parse_g1_curve_parameters(&bytes)?;
        let (num_pairs_encoding, _) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
        let num_pairs = num_pairs_encoding[0] as f64;
        if num_pairs == 0f64 {
            return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
        }

        let num_limbs = num_limbs_for_modulus(&modulus)? as f64;
        let num_order_words = num_units_for_group_order(&order)? as f64;

        let inner_term = evaluate_poly(&vec![
            0f64,
            FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER,
            FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER,
        ], num_limbs)?;
        let order_words_factor = num_order_words * ORDER_UNITS_LINEAR_TERM;
        let num_pair_factor = order_words_factor * num_pairs;
        let num_pair_factor = num_pair_factor* G1_MULTIEXP_DISCOUNT;
        let inner_term = inner_term * num_pair_factor;

        let cubic_term = evaluate_poly(&vec![
            0f64,
            0f64,
            0f64,
            FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER
        ], num_limbs)?;

        let price = G1_MULTIEXP_CONSTANT;
        let price = price + inner_term;
        let price = price + cubic_term;

        let price = price as u64;
        return Ok(price)
    }
}