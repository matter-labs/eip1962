mod parsers;
mod utils;
mod meter_arith;

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
        // OperationType::G1MULTIEXP => {
        //     PublicG1Api::multiexp(&input)
        // },
        OperationType::G2ADD => {
            meter_addition_g2(&input)
        },
        OperationType::G2MUL => {
            meter_multiplication_g2(&input)
        },
        _ => {
            unimplemented!()
        }
        // OperationType::G2MULTIEXP => {
        //     PublicG2Api::multiexp(&input)
        // },
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
            // OPERATION_G1_ADD | OPERATION_G1_MUL | OPERATION_G1_MULTIEXP => {
            //     Self::meter_in_base_field(operation, &rest)
            // },
            // OPERATION_G2_ADD | OPERATION_G2_MUL | OPERATION_G2_MULTIEXP => {
            //     Self::meter_in_extension(operation, &rest)
            // },
            // OPERATION_PAIRING => {
            //     Ok(1_000_000u64)
            // },
            _ => {
                Err(ApiError::InputError("Unknown operation type".to_owned()))
            }
        };

        result
    }
}