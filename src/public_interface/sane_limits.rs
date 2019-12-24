pub(crate) const MAX_BLS12_X_BIT_LENGTH: usize = 128;
pub(crate) const MAX_BN_U_BIT_LENGTH: usize = 128;

pub(crate) const MAX_BLS12_X_HAMMING: u32 = 128u32;
pub(crate) const MAX_BN_SIX_U_PLUS_TWO_HAMMING: u32 = 128u32;

pub(crate) const MAX_ATE_PAIRING_ATE_LOOP_COUNT: usize = 2032;
pub(crate) const MAX_ATE_PAIRING_ATE_LOOP_COUNT_HAMMING: u32 = 2032u32;

pub(crate) const MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH: usize = 2032;
pub(crate) const MAX_ATE_PAIRING_FINAL_EXP_W1_BIT_LENGTH: usize = 2032;

pub(crate) const MAX_LOOP_PARAMETERS_BYTE_LEN: usize = MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH / 8;

use static_assertions::const_assert;
use crate::constants::*;

const_assert!(std::mem::size_of::<MaxLoopParametersUint>() >= MAX_LOOP_PARAMETERS_BYTE_LEN);