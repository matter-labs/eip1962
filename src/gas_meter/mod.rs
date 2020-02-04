mod parsers;
mod utils;
mod meter_arith;
mod meter_pairing;

extern crate serde;
extern crate serde_json;
extern crate once_cell;

use crate::errors::ApiError;
use crate::public_interface::decode_utils::*;
use crate::public_interface::constants::*;
use self::parsers::*;
use crate::public_interface::OperationType;

pub struct GasMeter;

// This is pure rust API
pub fn meter_operation(operation: OperationType, input: &[u8]) -> Result<u64, ApiError> {
    match operation {
        OperationType::G1ADD => {
            meter_addition_g1(&input)
        },
        OperationType::G1MUL => {
            meter_multiplication_g1(&input)
        },
        OperationType::G1MULTIEXP => {
            meter_multiexp_g1(&input)
        },
        OperationType::G2ADD => {
            meter_addition_g2(&input)
        },
        OperationType::G2MUL => {
            meter_multiplication_g2(&input)
        },
        OperationType::G2MULTIEXP => {
            meter_multiexp_g2(&input)
        },
        OperationType::MNT4PAIR => {
            meter_mnt4(&input)
        },
        OperationType::MNT6PAIR => {
            meter_mnt6(&input)
        },
        OperationType::BLS12PAIR => {
            meter_bls12(&input)
        },
        OperationType::BNPAIR => {
            meter_bn(&input)
        } 
    }
}

fn meter_addition_g1(input: &[u8]) -> Result<u64, ApiError> {

    let (modulus, modulus_len, _, rest) = parse_g1_curve_parameters(&input)?;
    if rest.len() != modulus_len * 4 {
        return Err(ApiError::InputError("Input is either too short or contains garbage for g1 addition metering".to_owned()));
    }
    let modulus_limbs = num_limbs_for_modulus(&modulus)?;

    let params = &*meter_arith::G1_ADDITION_PARAMS_INSTANCE;

    meter_arith::meter_addition(modulus_limbs, params)
}

fn meter_addition_g2(input: &[u8]) -> Result<u64, ApiError> {

    let (modulus, modulus_len, _, ext_degree, rest) = parse_g2_curve_parameters(&input)?;
    if rest.len() != modulus_len * 4 * (ext_degree as usize) {
        return Err(ApiError::InputError("Input is either too short or contains garbage for g2 addition metering".to_owned()));
    }
    let modulus_limbs = num_limbs_for_modulus(&modulus)?;

    let params = if ext_degree == EXTENSION_DEGREE_2 {
        &*meter_arith::G2_EXT_2_ADDITION_PARAMS_INSTANCE
    } else if ext_degree == EXTENSION_DEGREE_3 {
        &*meter_arith::G2_EXT_3_ADDITION_PARAMS_INSTANCE
    } else {
        unreachable!();
    };

    meter_arith::meter_addition(modulus_limbs, params)
}


fn meter_multiplication_g1(input: &[u8]) -> Result<u64, ApiError> {
    let (modulus, modulus_len, order_len, rest) = parse_g1_curve_parameters(&input)?;
    if rest.len() != modulus_len * 2 + order_len {
        return Err(ApiError::InputError("Input is either too short or contains garbage for g1 multiplication metering".to_owned()));
    }
    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    // let order_limbs = num_units_for_group_order(&order)?;
    let order_limbs = num_units_for_group_order_length(order_len)?;

    let params = &*meter_arith::G1_MULTIPLICATION_PARAMS_INSTANCE;

    meter_arith::meter_multiplication(modulus_limbs, order_limbs, params, true)
}

fn meter_multiplication_g2(input: &[u8]) -> Result<u64, ApiError> {
    let (modulus, modulus_len, order_len, ext_degree, rest) = parse_g2_curve_parameters(&input)?;
    if rest.len() != modulus_len * 2 * (ext_degree as usize) + order_len {
        return Err(ApiError::InputError("Input is either too short or contains garbage for g2 multiplication metering".to_owned()));
    }

    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    // let order_limbs = num_units_for_group_order(&order)?;
    let order_limbs = num_units_for_group_order_length(order_len)?;

    let params = if ext_degree == EXTENSION_DEGREE_2 {
        &*meter_arith::G2_EXT_2_MULTIPLICATION_PARAMS_INSTANCE
    } else if ext_degree == EXTENSION_DEGREE_3 {
        &*meter_arith::G2_EXT_3_MULTIPLICATION_PARAMS_INSTANCE
    } else {
        unreachable!();
    };

    meter_arith::meter_multiplication(modulus_limbs, order_limbs, params, true)
}

fn meter_multiexp_g1(input: &[u8]) -> Result<u64, ApiError> {
    let (modulus, modulus_len, order_len, rest) = parse_g1_curve_parameters(&input)?;
    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    // let order_limbs = num_units_for_group_order(&order)?;
    let order_limbs = num_units_for_group_order_length(order_len)?;

    let (num_pairs_encoding, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
    let num_pairs = num_pairs_encoding[0] as usize;

    if num_pairs == 0 {
        return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
    }

    if rest.len() != num_pairs * (modulus_len * 2 +  order_len) {
        return Err(ApiError::InputError("Input is either too short or contains garbage for g1 multiexp metering".to_owned()));
    }

    let params = &*meter_arith::G1_MULTIPLICATION_PARAMS_INSTANCE;
    let discounts = &*meter_arith::MULTIEXP_PARAMS_INSTANCE;

    meter_arith::meter_multiexp(modulus_limbs, order_limbs, num_pairs, params, discounts)
}

fn meter_multiexp_g2(input: &[u8]) -> Result<u64, ApiError> {
    let (modulus, modulus_len, order_len, ext_degree, rest) = parse_g2_curve_parameters(&input)?;

    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    // let order_limbs = num_units_for_group_order(&order)?;
    let order_limbs = num_units_for_group_order_length(order_len)?;

    let params = if ext_degree == EXTENSION_DEGREE_2 {
        &*meter_arith::G2_EXT_2_MULTIPLICATION_PARAMS_INSTANCE
    } else if ext_degree == EXTENSION_DEGREE_3 {
        &*meter_arith::G2_EXT_3_MULTIPLICATION_PARAMS_INSTANCE
    } else {
        unreachable!();
    };

    let (num_pairs_encoding, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
    let num_pairs = num_pairs_encoding[0] as usize;

    if num_pairs == 0 {
        return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
    }

    if rest.len() != num_pairs * (modulus_len * 2 * (ext_degree as usize) +  order_len) {
        return Err(ApiError::InputError("Input is either too short or contains garbage for g2 multiexp metering".to_owned()));
    }

    let discounts = &*meter_arith::MULTIEXP_PARAMS_INSTANCE;

    meter_arith::meter_multiexp(modulus_limbs, order_limbs, num_pairs, params, discounts)
}

fn meter_bls12(input: &[u8]) -> Result<u64, ApiError> {
    self::meter_pairing::meter_bls12_pairing(input, &*self::meter_pairing::BLS12_PARAMS_INSTANCE, self::meter_pairing::BLS12_MAX_MODULUS_POWER)
}

fn meter_bn(input: &[u8]) -> Result<u64, ApiError> {
    self::meter_pairing::meter_bn_pairing(input, &*self::meter_pairing::BN_PARAMS_INSTANCE, self::meter_pairing::BN_MAX_MODULUS_POWER)
}

fn meter_mnt4(input: &[u8]) -> Result<u64, ApiError> {
    self::meter_pairing::meter_mnt_pairing(
        input, 
        &*self::meter_pairing::MNT4_PARAMS_INSTANCE, 
        self::meter_pairing::MNT4_MAX_MODULUS_POWER,
        2
    )
}

fn meter_mnt6(input: &[u8]) -> Result<u64, ApiError> {
    self::meter_pairing::meter_mnt_pairing(
        input, 
        &*self::meter_pairing::MNT6_PARAMS_INSTANCE, 
        self::meter_pairing::MNT6_MAX_MODULUS_POWER,
        3
    )
}

impl GasMeter {
    pub fn meter(bytes: &[u8]) -> Result<u64, ApiError> {
        let (op_type, rest) = split(bytes, OPERATION_ENCODING_LENGTH , "Input should be longer than operation type encoding")?;
        let operation = op_type[0];
        let result = match operation {
            OPERATION_G1_ADD => {
                meter_addition_g1(&rest)
            },
            OPERATION_G2_ADD => {
                meter_addition_g2(&rest)
            },
            OPERATION_G1_MUL => {
                meter_multiplication_g1(&rest)
            },
            OPERATION_G2_MUL => {
                meter_multiplication_g2(&rest)
            }
            OPERATION_G1_MULTIEXP => {
                meter_multiexp_g1(&rest)
            },
            OPERATION_G2_MULTIEXP => {
                meter_multiexp_g2(&rest)
            },
            OPERATION_PAIRING => {
                let (curve_type, rest) = split(rest, CURVE_TYPE_LENGTH, "Input should be longer than curve type encoding")?;

                match curve_type[0] {
                    BLS12 => {
                        meter_bls12(&rest)
                    },
                    BN => {
                        meter_bn(&rest)
                    },
                    MNT4 => {
                        meter_mnt4(&rest)
                    },
                    MNT6 => {
                        meter_mnt6(&rest)
                    },
                    _ => {
                        return Err(ApiError::InputError("Unknown curve type".to_owned()));
                    }
                }
            },
            _ => {
                Err(ApiError::InputError("Unknown operation type".to_owned()))
            }
        };

        result
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn test_calculate_example_prices_mnt4_753() {
        use crate::test::pairings::mnt4::assemble_mnt4_753;
        use crate::public_interface::OperationType;

        let calldata = assemble_mnt4_753(4);

        let price = super::meter_operation(OperationType::MNT4PAIR, &calldata[1..]).unwrap();

        println!("MNT4-753 for 4 pairs = {}", price);
        
    }
}