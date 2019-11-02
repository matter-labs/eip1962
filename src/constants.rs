use lazy_static::*;
use num_bigint::BigUint;
use num_traits::{Zero, One};

lazy_static! {
    pub(crate) static ref ZERO_BIGUINT: BigUint = BigUint::zero();
    pub(crate) static ref ONE_BIGUINT: BigUint = BigUint::from(1u64);
    pub(crate) static ref TWO_BIGUINT: BigUint = BigUint::from(2u64);
    pub(crate) static ref THREE_BIGUINT: BigUint = BigUint::from(3u64);
    pub(crate) static ref FOUR_BIGUINT: BigUint = BigUint::from(4u64);
    pub(crate) static ref SIX_BIGUINT: BigUint = BigUint::from(6u64);
    pub(crate) static ref BIGUINT_TWO_IN_64: BigUint = { BigUint::one() << 64 };
}

pub const MIN_FIELD_LIMBS: usize = 4;
pub const MAX_FIELD_LIMBS: usize = 16;
pub const MIN_GROUP_LIMBS: usize = 1;
pub const MAX_GROUP_LIMBS: usize = 16;