// use lazy_static::*;

pub const MIN_FIELD_LIMBS: usize = 4;
pub const MAX_FIELD_LIMBS: usize = 16;
pub const MIN_GROUP_LIMBS: usize = 1;
pub const MAX_GROUP_LIMBS: usize = 16;
pub const MAX_MODULUS_BYTE_LEN: usize = MAX_FIELD_LIMBS * 8;
pub const MAX_GROUP_BYTE_LEN: usize = MAX_GROUP_LIMBS * 8;

pub(crate) const MAX_FIELD_SQUARED_LIMBS: usize = MAX_FIELD_LIMBS * 2 + 1;
pub(crate) const MAX_FIELD_FROBENIUS_LIMBS_EXT_12: usize = MAX_FIELD_LIMBS * 12;
pub(crate) const MAX_FIELD_FROBENIUS_LIMBS_EXT_4: usize = MAX_FIELD_LIMBS * 4;
pub(crate) const MAX_FIELD_FROBENIUS_LIMBS_EXT_6: usize = MAX_FIELD_LIMBS * 6;

pub use crate::fixed_width_field::{MaxFieldUint, MaxFieldSquaredUint};
pub use crate::fixed_width_fp3_fp4::{MaxFrobeniusFp3, MaxFrobeniusFp4};
pub use crate::fixed_width_fp6::{MaxFrobeniusFp6};
pub use crate::fixed_width_fp12::{MaxFrobeniusFp12};
pub use crate::fixed_width_group_and_loop::{MaxGroupSizeUint, MaxLoopParametersUint};

// pub type GroupRepresentationArray = crate::arrayvec::ArrayVec<[u64; MAX_GROUP_LIMBS]>;
// pub type FieldRepresentationArray = crate::arrayvec::ArrayVec<[u64; MAX_FIELD_LIMBS]>;
// pub type LoopCountRepresentationArray = crate::arrayvec::ArrayVec<[u64; (crate::public_interface::sane_limits::MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH / 64) + 1]>;

// TODO: add constant asserts
