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

// TODO: add constant asserts

// mod field;
// mod field_construction;
// mod group;
// mod fp3;
// mod fp4;
// mod fp6;
// mod fp12;
// mod loop_param;

// pub use self::field::MaxFieldUint;
// pub use self::field_construction::MaxFieldSquaredUint;
// pub use self::group::MaxGroupSizeUint;
// pub use self::fp3::;
// pub use self::fp4::;
// pub use self::fp6::MaxFrobeniusFp6;
// pub use self::fp12::MaxFrobeniusFp12;
// pub use self::loop_param::MaxLoopParametersUint;








// lazy_static! {
//     pub(crate) static ref ZERO_F_SQUARED: MaxFieldSquaredUint = MaxFieldSquaredUint::zero();
//     pub(crate) static ref ONE_F_SQUARED: MaxFieldSquaredUint = MaxFieldSquaredUint::one();
//     pub(crate) static ref TWO_F_SQUARED: MaxFieldSquaredUint = MaxFieldSquaredUint::from(2);
//     pub(crate) static ref THREE_F_SQUARED: MaxFieldSquaredUint = MaxFieldSquaredUint::from(3);
//     pub(crate) static ref FOUR_F_SQUARED: MaxFieldSquaredUint = MaxFieldSquaredUint::from(4);
//     pub(crate) static ref SIX_F_SQUARED: MaxFieldSquaredUint = MaxFieldSquaredUint::from(6);
//     // pub(crate) static ref ONE_BIGUINT: BigUint = BigUint::from(1u64);
//     // pub(crate) static ref TWO_BIGUINT: BigUint = BigUint::from(2u64);
//     // pub(crate) static ref THREE_BIGUINT: BigUint = BigUint::from(3u64);
//     // pub(crate) static ref FOUR_BIGUINT: BigUint = BigUint::from(4u64);
//     // pub(crate) static ref SIX_BIGUINT: BigUint = BigUint::from(6u64);
//     // pub(crate) static ref BIGUINT_TWO_IN_64: BigUint = { BigUint::one() << 64 };
//     // pub(crate) static ref ZERO_BIGUINT: BigUint = BigUint::zero();
//     // pub(crate) static ref ONE_BIGUINT: BigUint = BigUint::from(1u64);
//     // pub(crate) static ref TWO_BIGUINT: BigUint = BigUint::from(2u64);
//     // pub(crate) static ref THREE_BIGUINT: BigUint = BigUint::from(3u64);
//     // pub(crate) static ref FOUR_BIGUINT: BigUint = BigUint::from(4u64);
//     // pub(crate) static ref SIX_BIGUINT: BigUint = BigUint::from(6u64);
//     // pub(crate) static ref BIGUINT_TWO_IN_64: BigUint = { BigUint::one() << 64 };
// }