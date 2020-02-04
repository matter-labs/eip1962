pub const BYTES_FOR_LENGTH_ENCODING: usize = 1;

pub const CURVE_TYPE_LENGTH: usize = 1;
pub const BLS12: u8 = 0x01;
pub const BN: u8 = 0x02;
pub const MNT4: u8 = 0x03;
pub const MNT6: u8 = 0x04;

pub const TWIST_TYPE_LENGTH: usize = 1;
pub const TWIST_TYPE_M: u8 = 0x01;
pub const TWIST_TYPE_D: u8 = 0x02;

pub const SIGN_ENCODING_LENGTH: usize = 1;
pub const SIGN_PLUS: u8 = 0x00;
pub const SIGN_MINUS: u8 = 0x01;

pub const BOOLEAN_ENCODING_LENGTH: usize = 1;
pub const BOOLEAN_FALSE: u8 = 0x00;
pub const BOOLEAN_TRUE: u8 = 0x01;

pub const EXTENSION_DEGREE_ENCODING_LENGTH: usize = 1;
pub const EXTENSION_DEGREE_2: u8 = 0x02;
pub const EXTENSION_DEGREE_3: u8 = 0x03;

pub const OPERATION_ENCODING_LENGTH: usize = 1;

pub const OPERATION_G1_ADD: u8 = 0x01;
pub const OPERATION_G1_MUL: u8 = 0x02;
pub const OPERATION_G1_MULTIEXP: u8 = 0x03;

pub const OPERATION_G2_ADD: u8 = 0x04;
pub const OPERATION_G2_MUL: u8 = 0x05;
pub const OPERATION_G2_MULTIEXP: u8 = 0x06;

pub const OPERATION_PAIRING: u8 = 0x07;

pub const NUM_LIMBS_MIN: usize = 4;
pub const NUM_LIMBS_MAX: usize = 16;
pub const NUM_GROUP_LIMBS_MIN: usize = 1;
pub const NUM_GROUP_LIMBS_MAX: usize = 16;

pub const MAX_MODULUS_BYTE_LEN: usize = 128;
pub const MAX_GROUP_BYTE_LEN: usize = 128;

use static_assertions::const_assert;
use crate::integers::*;

const_assert!(MAX_MODULUS_BYTE_LEN == NUM_LIMBS_MAX * 8);

const_assert!(MAX_GROUP_BYTE_LEN == NUM_GROUP_LIMBS_MAX * 8);

const_assert!(std::mem::size_of::<MaxFieldUint>() >= NUM_LIMBS_MAX * 8);
const_assert!(std::mem::size_of::<MaxFieldSquaredUint>() >= NUM_LIMBS_MAX * 8 * 2);

const_assert!(std::mem::size_of::<MaxGroupSizeUint>() >= NUM_GROUP_LIMBS_MAX * 8);
