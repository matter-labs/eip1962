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
pub fn perform_operation(operation: OperationType, input: &[u8]) -> Result<u64, ApiError> {
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
        _ => {
            unimplemented!()
        }
        // OperationType::BLS12PAIR | OperationType::BNPAIR | OperationType::MNT4PAIR | OperationType::MNT6PAIR) => {
        //     use crate::field::*;
        //     use crate::public_interface::decode_utils::*;

        //     let modulus_limbs = {
        //         let (_, modulus, _) = parse_modulus_and_length(&input)?;
        //         let modulus_limbs = num_limbs_for_modulus(&modulus)?;

        //         modulus_limbs
        //     };

        //     match operation {
        //         OperationType::BLS12PAIR => {
        //             let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, PairingApiImplementation, input, pair_bls12); 

        //             result
        //         },
        //         OperationType::BNPAIR => {
        //             let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, PairingApiImplementation, input, pair_bn); 

        //             result
        //         },
        //         OperationType::MNT4PAIR => {
        //             let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, PairingApiImplementation, input, pair_mnt4); 

        //             result
        //         },
        //         OperationType::MNT6PAIR => {
        //             let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, PairingApiImplementation, input, pair_mnt6); 

        //             result
        //         },

        //         _ => {
        //             unreachable!()
        //         }
        //     }
        // }
    }
}

fn meter_addition_g1(input: &[u8]) -> Result<u64, ApiError> {

    let (modulus, _, _) = parse_g1_curve_parameters(&input)?;
    let modulus_limbs = num_limbs_for_modulus(&modulus)?;

    let params = &*meter_arith::G1_ADDITION_PARAMS_INSTANCE;

    meter_arith::meter_addition(modulus_limbs, params)
}

fn meter_addition_g2(input: &[u8]) -> Result<u64, ApiError> {

    let (modulus, _, ext_degree, _) = parse_g2_curve_parameters(&input)?;

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

    let (modulus, order, _) = parse_g1_curve_parameters(&input)?;
    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    let order_limbs = num_units_for_group_order(&order)?;

    let params = &*meter_arith::G1_MULTIPLICATION_PARAMS_INSTANCE;

    meter_arith::meter_multiplication(modulus_limbs, order_limbs, params)
}

fn meter_multiplication_g2(input: &[u8]) -> Result<u64, ApiError> {

    let (modulus, order, ext_degree, _) = parse_g2_curve_parameters(&input)?;

    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    let order_limbs = num_units_for_group_order(&order)?;

    let params = if ext_degree == EXTENSION_DEGREE_2 {
        &*meter_arith::G2_EXT_2_MULTIPLICATION_PARAMS_INSTANCE
    } else if ext_degree == EXTENSION_DEGREE_3 {
        &*meter_arith::G2_EXT_3_MULTIPLICATION_PARAMS_INSTANCE
    } else {
        unreachable!();
    };

    meter_arith::meter_multiplication(modulus_limbs, order_limbs, params)
}

fn meter_multiexp_g1(input: &[u8]) -> Result<u64, ApiError> {

    let (modulus, order, rest) = parse_g1_curve_parameters(&input)?;
    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    let order_limbs = num_units_for_group_order(&order)?;

    let (num_pairs_encoding, rest) = split(rest, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get number of pairs")?;
    let num_pairs = num_pairs_encoding[0] as usize;

    if num_pairs == 0 {
        return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
    }

    if rest.len() == 0 {
        return Err(ApiError::InputError("Input is not long enough".to_owned()));
    }

    let params = &*meter_arith::G1_MULTIPLICATION_PARAMS_INSTANCE;
    let discounts = &*meter_arith::MULTIEXP_PARAMS_INSTANCE;

    meter_arith::meter_multiexp(modulus_limbs, order_limbs, num_pairs, params, discounts)
}

fn meter_multiexp_g2(input: &[u8]) -> Result<u64, ApiError> {

    let (modulus, order, ext_degree, rest) = parse_g2_curve_parameters(&input)?;

    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    let order_limbs = num_units_for_group_order(&order)?;

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

    if rest.len() == 0 {
        return Err(ApiError::InputError("Input is not long enough".to_owned()));
    }

    let discounts = &*meter_arith::MULTIEXP_PARAMS_INSTANCE;

    meter_arith::meter_multiexp(modulus_limbs, order_limbs, num_pairs, params, discounts)
}

fn meter_mnt4(input: &[u8]) -> Result<u64, ApiError> {
    self::meter_pairing::meter_mnt_pairing(input, &*self::meter_pairing::MNT4_PARAMS_INSTANCE, self::meter_pairing::MNT4_MAX_MODULUS_POWER)
}

fn meter_mnt6(input: &[u8]) -> Result<u64, ApiError> {
    self::meter_pairing::meter_mnt_pairing(input, &*self::meter_pairing::MNT6_PARAMS_INSTANCE, self::meter_pairing::MNT6_MAX_MODULUS_POWER)
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
                        unimplemented!()
                    },
                    BN => {
                        unimplemented!()
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