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