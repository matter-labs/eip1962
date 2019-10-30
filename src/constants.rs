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