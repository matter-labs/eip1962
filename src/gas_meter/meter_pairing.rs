use serde::{Deserialize, Deserializer};
use std::collections::HashMap;
use crate::errors::ApiError;

use crate::once_cell::sync::Lazy;

use super::meter_arith::*;
use super::parsers::*;

use crate::public_interface::decode_utils::*;

pub(crate) const MNT4_MAX_MODULUS_POWER: usize = 4;
pub(crate) const MNT6_MAX_MODULUS_POWER: usize = 6;
pub(crate) const BN_MAX_MODULUS_POWER: usize = 6;
pub(crate) const BLS12_MAX_MODULUS_POWER: usize = 6;

#[derive(Clone, Deserialize, Debug)]
pub(crate) struct MntPairingParams {
    #[serde(deserialize_with = "parse_tuple_usize_u64")]
    #[serde(rename = "one_off")]
    one_off: Vec<(usize, u64)>,

    #[serde(deserialize_with = "parse_u64")]
    #[serde(rename = "multiplier")]
    discount_multiplier: u64,

    miller_features: Vec<(String, u64)>,

    miller: Vec<(u64, Vec<(usize, usize)>)>,

    final_exp_features: Vec<(String, u64)>,

    final_exp: Vec<(u64, Vec<(usize, usize)>)>
}

static MNT4_PARAMS_JSON: &'static str = include_str!("mnt4_model.json");
static MNT6_PARAMS_JSON: &'static str = include_str!("mnt6_model.json");

pub(crate) static MNT4_PARAMS_INSTANCE: Lazy<MntPairingParams> = Lazy::new(|| {
    crate::serde_json::from_str(MNT4_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) static MNT6_PARAMS_INSTANCE: Lazy<MntPairingParams> = Lazy::new(|| {
    crate::serde_json::from_str(MNT6_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) fn meter_mnt_pairing(input: &[u8], params: &MntPairingParams, max_power: usize) -> Result<u64, ApiError> {
    const GROUP_LIMBS_INDEX: usize = 0;
    const ATE_LOOP_BITS_INDEX: usize = 1;
    const ATE_LOOP_HAMMING_INDEX: usize = 2;
    const EXP_W0_LOOP_BITS_INDEX: usize = 3;
    const EXP_W0_HAMMING_INDEX: usize = 4;
    const EXP_W1_LOOP_BITS_INDEX: usize = 3;
    const EXP_W1_HAMMING_INDEX: usize = 4;

    debug_assert!(max_power == MNT4_MAX_MODULUS_POWER || max_power == MNT6_MAX_MODULUS_POWER);

    let (
        modulus, 
        order, 
        num_pairs, 
        (ate_loop_bits, ate_loop_hamming), 
        (exp_w0_bits, exp_w0_hamming),
        (exp_w1_bits, exp_w1_hamming),
        _
    ) = parse_mnt_pairing_parameters(&input)?;

    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    let order_limbs = num_units_for_group_order(&order)?;

    let mut one_off: Vec<_> = params.one_off.iter().filter(|(limbs, _)| *limbs == modulus_limbs).collect();
    if one_off.len() != 1 {
        return Err(ApiError::UnknownParameter(format!("Unknown number of limbs = {}", modulus_limbs)));
    }

    let one_off = one_off.pop().expect("result exists").1;

    let modulus_limbs_powers = make_powers(modulus_limbs as u64, max_power)?;
    let params_vector = vec![order_limbs as u64, ate_loop_bits, ate_loop_hamming, exp_w0_bits, exp_w0_hamming, exp_w1_bits, exp_w1_hamming];

    let miller_cost = {
        let miller_params = vec![
            &params_vector[GROUP_LIMBS_INDEX..(GROUP_LIMBS_INDEX+1)], 
            &params_vector[ATE_LOOP_BITS_INDEX..(ATE_LOOP_BITS_INDEX+1)], 
            &params_vector[ATE_LOOP_HAMMING_INDEX..(ATE_LOOP_HAMMING_INDEX+1)], 
            &modulus_limbs_powers[..] 
            ];
        let mut miller_cost = eval_model(&params.miller, &miller_params)?;
        miller_cost = miller_cost.checked_mul(num_pairs as u64).ok_or(ApiError::Overflow)?;

        miller_cost
    };

    let final_exp_cost = {
        let final_exp_params = vec![
            &params_vector[EXP_W0_LOOP_BITS_INDEX..(EXP_W0_LOOP_BITS_INDEX+1)], 
            &params_vector[EXP_W0_HAMMING_INDEX..(EXP_W0_HAMMING_INDEX+1)], 
            &params_vector[EXP_W1_LOOP_BITS_INDEX..(EXP_W1_LOOP_BITS_INDEX+1)], 
            &params_vector[EXP_W1_HAMMING_INDEX..(EXP_W1_HAMMING_INDEX+1)], 
            &modulus_limbs_powers[..] 
            ];
        let final_exp_cost = eval_model(&params.final_exp, &final_exp_params)?;

        final_exp_cost
    };

    let mut result = one_off;
    result = result.checked_add(miller_cost).ok_or(ApiError::Overflow)?;
    result = result.checked_add(final_exp_cost).ok_or(ApiError::Overflow)?;

    Ok(result)
}

fn eval_model(
    coeffs_variables_and_powers: &[(u64, Vec<(usize, usize)>)],
    variables: &[ &[u64] ]
) -> Result<u64, ApiError> {
    let mut final_result = 0u64;
    if coeffs_variables_and_powers.len() != variables.len() {
        return Err(ApiError::MissingValue);
    }

    for (coeff, var_and_power) in coeffs_variables_and_powers.iter() {
        let mut subpart = *coeff;
        for (variable, power) in var_and_power.iter() {
            let variable_powers = variables.get(*variable).ok_or(ApiError::MissingValue)?;
            let variable_power_value = variable_powers.get(*power - 1).ok_or(ApiError::MissingValue)?;
            subpart = subpart.checked_mul(*variable_power_value).ok_or(ApiError::Overflow)?;
        }
        final_result = final_result.checked_add(subpart).ok_or(ApiError::Overflow)?;
    }

    Ok(final_result)
}

fn make_powers(value: u64, required_power: usize) -> Result<Vec<u64>, ApiError> {
    debug_assert!(required_power > 0);
    let mut powers = Vec::with_capacity(required_power);
    let mut p = 1u64;
    for _ in 1..=required_power {
        p = p.checked_mul(value).ok_or(ApiError::Overflow)?;
        powers.push(p);
    }

    Ok(powers)
}

#[cfg(test)]
mod test {
    #[test]
    fn test_pairing_params_deserialization() {
        let t = &*super::MNT4_PARAMS_INSTANCE;
        println!("Params MNT4 = {:?}", t);

        let t = &*super::MNT6_PARAMS_INSTANCE;
        println!("Params MNT6 = {:?}", t);

        // let t = &*super::G2_EXT_2_ADDITION_PARAMS_INSTANCE;
        // println!("Params G2 add ext 2= {:?}", t);

        // let t = &*super::G2_EXT_3_ADDITION_PARAMS_INSTANCE;
        // println!("Params G2 add ext 3 = {:?}", t);

        // let t = &*super::G1_MULTIPLICATION_PARAMS_INSTANCE;
        // println!("Params G1 mul = {:?}", t);

        // let t = &*super::G2_EXT_2_MULTIPLICATION_PARAMS_INSTANCE;
        // println!("Params G2 mul ext 2 = {:?}", t);

        // let t = &*super::G2_EXT_3_MULTIPLICATION_PARAMS_INSTANCE;
        // println!("Params G2 mul ext 3 = {:?}", t);

        // let t = &*super::MULTIEXP_PARAMS_INSTANCE;
        // println!("Multiexp discounts = {:?}", t);
    }
}