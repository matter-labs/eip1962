use crate::public_interface::pairing_ops::PairingApiImplementation;
use crate::public_interface::g1_ops::{G1Api, PublicG1Api};
use crate::public_interface::g2_ops::{G2Api, PublicG2Api};

use crate::errors::ApiError;

// For C style API caller has to preallocate some buffers for results 
pub const PREALLOCATE_FOR_ERROR_BYTES: usize = 256;
pub const PREALLOCATE_FOR_RESULT_BYTES: usize = 768;

use static_assertions::const_assert;
const_assert!(PREALLOCATE_FOR_RESULT_BYTES == crate::public_interface::constants::MAX_MODULUS_BYTE_LEN * 3 * 2);

#[repr(u8)]
#[derive(Copy, Clone, Debug)]
pub enum OperationType {
    G1ADD = 1,
    G1MUL = 2,
    G1MULTIEXP = 3,
    G2ADD = 4,
    G2MUL = 5,
    G2MULTIEXP = 6,
    BLS12PAIR = 7,
    BNPAIR = 8,
    MNT4PAIR = 9,
    MNT6PAIR = 10,
}

impl OperationType {
    pub fn from_u8(value: u8) -> Option<Self> {
        match value {
            G1ADD_OPERATION_RAW_VALUE => {
                Some(OperationType::G1ADD)
            },
            G1MUL_OPERATION_RAW_VALUE => {
                Some(OperationType::G1MUL)
            },
            G1MULTIEXP_OPERATION_RAW_VALUE => {
                Some(OperationType::G1MULTIEXP)
            },
            G2ADD_OPERATION_RAW_VALUE => {
                Some(OperationType::G2ADD)
            },
            G2MUL_OPERATION_RAW_VALUE => {
                Some(OperationType::G2MUL)
            },
            G2MULTIEXP_OPERATION_RAW_VALUE => {
                Some(OperationType::G2MULTIEXP)
            },
            BLS12PAIR_OPERATION_RAW_VALUE => {
                Some(OperationType::BLS12PAIR)
            },
            BNPAIR_OPERATION_RAW_VALUE => {
                Some(OperationType::BNPAIR)
            },
            MNT4PAIR_OPERATION_RAW_VALUE => {
                Some(OperationType::MNT4PAIR)
            },
            MNT6PAIR_OPERATION_RAW_VALUE => {
                Some(OperationType::MNT6PAIR)
            },
            _ => {
                None
            }
        }
    }

    pub fn as_u8(&self) -> u8 {
        *self as u8
    }
}

pub const G1ADD_OPERATION_RAW_VALUE: u8 = OperationType::G1ADD as u8;
pub const G1MUL_OPERATION_RAW_VALUE: u8 = OperationType::G1MUL as u8;
pub const G1MULTIEXP_OPERATION_RAW_VALUE: u8 = OperationType::G1MULTIEXP as u8;

pub const G2ADD_OPERATION_RAW_VALUE: u8 = OperationType::G2ADD as u8;
pub const G2MUL_OPERATION_RAW_VALUE: u8 = OperationType::G2MUL as u8;
pub const G2MULTIEXP_OPERATION_RAW_VALUE: u8 = OperationType::G2MULTIEXP as u8;

pub const BLS12PAIR_OPERATION_RAW_VALUE: u8 = OperationType::BLS12PAIR as u8;
pub const BNPAIR_OPERATION_RAW_VALUE: u8 = OperationType::BNPAIR as u8;
pub const MNT4PAIR_OPERATION_RAW_VALUE: u8 = OperationType::MNT4PAIR as u8;
pub const MNT6PAIR_OPERATION_RAW_VALUE: u8 = OperationType::MNT6PAIR as u8;

// This is pure rust API
pub fn perform_operation(operation: OperationType, input: &[u8]) -> Result<Vec<u8>, ApiError> {
    match operation {
        OperationType::G1ADD => {
            PublicG1Api::add_points(&input)
        },
        OperationType::G1MUL => {
            PublicG1Api::mul_point(&input)
        },
        OperationType::G1MULTIEXP => {
            PublicG1Api::multiexp(&input)
        },
        OperationType::G2ADD => {
            PublicG2Api::add_points(&input)
        },
        OperationType::G2MUL => {
            PublicG2Api::mul_point(&input)
        },
        OperationType::G2MULTIEXP => {
            PublicG2Api::multiexp(&input)
        },
        OperationType::BLS12PAIR | OperationType::BNPAIR | OperationType::MNT4PAIR | OperationType::MNT6PAIR => {
            use crate::field::*;
            use crate::public_interface::decode_utils::*;

            let modulus_limbs = {
                let (_, modulus, _) = parse_modulus_and_length(&input)?;
                let modulus_limbs = num_limbs_for_modulus(&modulus)?;

                modulus_limbs
            };

            match operation {
                OperationType::BLS12PAIR => {
                    let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, PairingApiImplementation, input, pair_bls12); 

                    result
                },
                OperationType::BNPAIR => {
                    let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, PairingApiImplementation, input, pair_bn); 

                    result
                },
                OperationType::MNT4PAIR => {
                    let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, PairingApiImplementation, input, pair_mnt4); 

                    result
                },
                OperationType::MNT6PAIR => {
                    let result: Result<Vec<u8>, ApiError> = expand_for_modulus_limbs!(modulus_limbs, PairingApiImplementation, input, pair_mnt6); 

                    result
                },

                _ => {
                    unreachable!()
                }
            }
        }
    }
}