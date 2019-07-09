mod parsers;
mod utils;

use num_bigint::BigUint;

use crate::errors::ApiError;
use crate::public_interface::decode_utils::*;
use crate::public_interface::constants::*;
use self::parsers::*;
use self::utils::*;

pub struct GasMeter;

trait GroupGasMeteringParams {
    const ADD_FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER: f64;
    const ADD_FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER: f64;
    const ADD_FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER: f64;
    const MUL_FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER: f64;
    const MUL_FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER: f64;
    const MUL_FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER: f64;
    const ORDER_UNITS_LINEAR_TERM: f64;
    const ADDITION_CONSTANT: f64;
    const MULTIPLICATION_CONSTANT: f64;
    const MULTIEXP_CONSTANT: f64;
    const MULTIEXP_DISCOUNT: f64;
}

struct G1GasMeteringParams;

impl GroupGasMeteringParams for G1GasMeteringParams {
    const ADD_FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER: f64 = 2.5f64;
    const ADD_FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER: f64 = 2f64;
    const ADD_FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER: f64 = 1.1f64;
    const MUL_FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER: f64 = 45f64;
    const MUL_FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER: f64 = 14f64;
    const MUL_FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER: f64 = 1.2f64;
    const ORDER_UNITS_LINEAR_TERM: f64 = 1f64;
    const ADDITION_CONSTANT: f64 = 300f64;
    const MULTIPLICATION_CONSTANT: f64 = 100f64;
    const MULTIEXP_CONSTANT: f64 = 100f64;
    const MULTIEXP_DISCOUNT: f64 = 0.25f64;
}

struct G2GasMeteringParamsFp2;

impl GroupGasMeteringParams for G2GasMeteringParamsFp2 {
    const ADD_FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER: f64 = 15f64;
    const ADD_FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER: f64 = 12f64;
    const ADD_FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER: f64 = 7f64;
    const MUL_FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER: f64 = 100f64;
    const MUL_FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER: f64 = 50f64;
    const MUL_FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER: f64 = 50f64;
    const ORDER_UNITS_LINEAR_TERM: f64 = 1f64;
    const ADDITION_CONSTANT: f64 = 500f64;
    const MULTIPLICATION_CONSTANT: f64 = 300f64;
    const MULTIEXP_CONSTANT: f64 = 300f64;
    const MULTIEXP_DISCOUNT: f64 = 0.25f64;
}

struct G2GasMeteringParamsFp3;

impl GroupGasMeteringParams for G2GasMeteringParamsFp3 {
    const ADD_FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER: f64 = 30f64;
    const ADD_FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER: f64 = 24f64;
    const ADD_FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER: f64 = 14f64;
    const MUL_FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER: f64 = 30f64;
    const MUL_FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER: f64 = 24f64;
    const MUL_FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER: f64 = 14f64;
    const ORDER_UNITS_LINEAR_TERM: f64 = 12f64;
    const ADDITION_CONSTANT: f64 = 500f64;
    const MULTIPLICATION_CONSTANT: f64 = 300f64;
    const MULTIEXP_CONSTANT: f64 = 300f64;
    const MULTIEXP_DISCOUNT: f64 = 0.25f64;
}

impl GasMeter {
    pub fn meter(bytes: &[u8]) -> Result<u64, ApiError> {
        let (op_type, rest) = split(bytes, OPERATION_ENCODING_LENGTH , "Input should be longer than operation type encoding")?;
        let operation = op_type[0];
        let result = match operation {
            OPERATION_G1_ADD | OPERATION_G1_MUL | OPERATION_G1_MULTIEXP => {
                Self::meter_in_base_field(operation, &rest)
            },
            OPERATION_G2_ADD | OPERATION_G2_MUL | OPERATION_G2_MULTIEXP => {
                Self::meter_in_extension(operation, &rest)
            },
            // OPERATION_PAIRING => {
            //     PublicPairingApi::pair(&rest)
            // },
            _ => {
                Err(ApiError::InputError("Unknown operation type".to_owned()))
            }
        };

        result
    }

    fn meter_addition<P: GroupGasMeteringParams>(
        modulus: &BigUint
    ) -> Result<u64, ApiError> {
        use std::f64;
        // apriori formula:
        // CONSTANT + A*LIBMS^1 + B*LIMBS^2 + C*LIMBS^3
        // linear term account for additions and should be small
        // quadratic term accounts for multiplications and should give the largest contribution
        // cubic term reflects inversions

        // let (modulus, _, _) = parse_g1_curve_parameters(&bytes)?;
        let num_limbs = num_limbs_for_modulus(&modulus)? as f64;
        let price = evaluate_poly(&vec![
            P::ADDITION_CONSTANT,
            P::ADD_FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER,
            P::ADD_FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER,
            P::ADD_FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER
        ], num_limbs)?;
        let price = price as u64;
        return Ok(price)
    }

    fn meter_multiplication<P: GroupGasMeteringParams>(
        modulus: &BigUint,
        order: &BigUint
    ) -> Result<u64, ApiError> {
        use std::f64;
        // apriori formula:
        // CONSTANT + (A*LIBMS^1 + B*LIMBS^2)*ORDER_WORDS + C*LIMBS^3
        // linear term account for additions and should be small
        // quadratic term accounts for multiplications and should give the largest contribution
        // cubic term reflects inversions

        // let (modulus, order, _) = parse_g1_curve_parameters(&bytes)?;
        let num_limbs = num_limbs_for_modulus(&modulus)? as f64;
        let num_order_words = num_units_for_group_order(&order)? as f64;

        let inner_term = evaluate_poly(&vec![
            0f64,
            P::MUL_FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER,
            P::MUL_FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER,
        ], num_limbs)?;
        let order_words_factor = num_order_words * P::ORDER_UNITS_LINEAR_TERM;
        let inner_term = inner_term * order_words_factor;

        let cubic_term = evaluate_poly(&vec![
            0f64,
            0f64,
            0f64,
            P::MUL_FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER
        ], num_limbs)?;

        let price = P::MULTIPLICATION_CONSTANT;
        let price = price + inner_term;
        let price = price + cubic_term;
        let price = price as u64;

        return Ok(price)
    }

    fn meter_multiexp<P: GroupGasMeteringParams>(
        modulus: &BigUint,
        order: &BigUint,
        rest: &[u8]
    ) -> Result<u64, ApiError> {
        use std::f64;
        // apriori formula:
        // CONSTANT + (A*LIBMS^1 + B*LIMBS^2)*ORDER_WORDS*NUM_PAIR*DISCOUNT + C*LIMBS^3
        // linear term account for additions and should be small
        // quadratic term accounts for multiplications and should give the largest contribution
        // cubic term reflects inversions

        // let (modulus, order, rest) = parse_g1_curve_parameters(&bytes)?;
        let (num_pairs_encoding, _) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
        let num_pairs = num_pairs_encoding[0] as f64;
        if num_pairs == 0f64 {
            return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
        }

        let num_limbs = num_limbs_for_modulus(&modulus)? as f64;
        let num_order_words = num_units_for_group_order(&order)? as f64;

        let inner_term = evaluate_poly(&vec![
            0f64,
            P::MUL_FIELD_MODULUS_LIMBS_LINEAR_TERM_MULTIPLIER,
            P::MUL_FIELD_MODULUS_LIMBS_QUADRATIC_TERM_MULTIPLIER,
        ], num_limbs)?;
        let order_words_factor = num_order_words * P::ORDER_UNITS_LINEAR_TERM;
        let num_pair_factor = order_words_factor * num_pairs;
        let num_pair_factor = num_pair_factor* P::MULTIEXP_DISCOUNT;
        let inner_term = inner_term * num_pair_factor;

        let cubic_term = evaluate_poly(&vec![
            0f64,
            0f64,
            0f64,
            P::MUL_FIELD_MODULUS_LIMBS_CUBIC_TERM_MULTIPLIER
        ], num_limbs)?;

        let price = P::MULTIEXP_CONSTANT;
        let price = price + inner_term;
        let price = price + cubic_term;

        let price = price as u64;
        return Ok(price)
    }

    fn meter_in_base_field(operation: u8, bytes: &[u8]) -> Result<u64, ApiError> {
        let (modulus, order, rest) = parse_g1_curve_parameters(&bytes)?;
        let result = match operation {
            OPERATION_G1_ADD => {
                Self::meter_addition::<G1GasMeteringParams>(&modulus)
            },
            OPERATION_G1_MUL => {
                Self::meter_multiplication::<G1GasMeteringParams>(&modulus, &order)
            },
            OPERATION_G1_MULTIEXP => {
                Self::meter_multiexp::<G1GasMeteringParams>(&modulus, &order, &rest)
            },
            _ => {
                Err(ApiError::InputError("Unknown operation type".to_owned()))
            }
        };

        result
    }

    fn meter_in_extension(operation: u8, bytes: &[u8]) -> Result<u64, ApiError> {
        let (modulus, order, extension_degree, rest) = parse_g2_curve_parameters(&bytes)?;
        let result = match extension_degree {
            EXTENSION_DEGREE_2 => {
                match operation {
                    OPERATION_G2_ADD => {
                        Self::meter_addition::<G2GasMeteringParamsFp2>(&modulus)
                    },
                    OPERATION_G2_MUL => {
                        Self::meter_multiplication::<G2GasMeteringParamsFp2>(&modulus, &order)
                    },
                    OPERATION_G2_MULTIEXP => {
                        Self::meter_multiexp::<G2GasMeteringParamsFp2>(&modulus, &order, &rest)
                    },
                    _ => {
                        Err(ApiError::InputError("Unknown operation type".to_owned()))
                    }
                }
            },
            EXTENSION_DEGREE_3 => {
                match operation {
                    OPERATION_G2_ADD => {
                        Self::meter_addition::<G2GasMeteringParamsFp3>(&modulus)
                    },
                    OPERATION_G2_MUL => {
                        Self::meter_multiplication::<G2GasMeteringParamsFp3>(&modulus, &order)
                    },
                    OPERATION_G2_MULTIEXP => {
                        Self::meter_multiexp::<G2GasMeteringParamsFp3>(&modulus, &order, &rest)
                    },
                    _ => {
                        Err(ApiError::InputError("Unknown operation type".to_owned()))
                    }
                }
            },
            _ => {
                Err(ApiError::InputError("Unknown extension degree".to_owned()))
            }
        };

        result
    }
}