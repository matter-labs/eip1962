pub(crate) mod pairings;
pub(crate) mod g2_ops;
pub(crate) mod g1_ops;
pub(crate) mod parsers;

mod fields;
// mod fuzzing;
mod gas_meter;

use num_bigint::BigUint;
use num_traits::Zero;
use num_traits::cast::ToPrimitive;

use crate::errors::ApiError;

pub(crate) fn num_limbs_for_modulus(modulus: &BigUint) -> Result<usize, ApiError> {
    use crate::field::calculate_num_limbs;

    let modulus_limbs = calculate_num_limbs(modulus.bits()).map_err(|_| ApiError::InputError(format!("Modulus is too large, file {}, line {}", file!(), line!())) )?;

    Ok(modulus_limbs)
}

pub(crate) fn num_units_for_group_order(order: &BigUint) -> Result<usize, ApiError> {
    let limbs = (order.bits() + 63) / 64;
    if limbs > 16 {
        return Err(ApiError::InputError(format!("Group order is too large, file {}, line {}", file!(), line!())));
    }

    Ok(limbs)
}

pub(crate) fn calculate_num_limbs(modulus: &BigUint) -> Result<usize, ()> {
    let bitlength = modulus.bits();

    let mut num_limbs = (bitlength / 64) + 1;
    if num_limbs < 4 {
        num_limbs = 4;
    }

    if num_limbs > 16 {
        return Err(());
    }

    Ok(num_limbs)
}
    
pub(crate) fn biguint_to_u64_vec(mut v: BigUint) -> Vec<u64> {

    let m = BigUint::from(1u64) << 64;
    let mut ret = Vec::with_capacity((v.bits() / 64) + 1);

    while v > BigUint::zero() {
        ret.push((&v % &m).to_u64().expect("is guaranteed to fit"));
        v >>= 64;
    }

    ret
}