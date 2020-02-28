use serde::{Deserialize};
use crate::errors::ApiError;

use once_cell::sync::Lazy;
use serde_json;
use std::collections::HashMap;

use super::parsers::*;

use crate::public_interface::decode_utils::*;
use crate::public_interface::sane_limits::*;

pub(crate) const MNT4_MAX_MODULUS_POWER: usize = 4;
pub(crate) const MNT6_MAX_MODULUS_POWER: usize = 6;
pub(crate) const BN_MAX_MODULUS_POWER: usize = 6;
pub(crate) const BLS12_MAX_MODULUS_POWER: usize = 6;

#[derive(Clone, Deserialize, Debug)]
pub(crate) struct MntPairingParams {
    #[serde(deserialize_with = "parse_hashmap_usize_u64_from_ints")]
    one_off: HashMap<usize, u64>,

    multiplier: u64,

    miller_features: Vec<(String, u64)>,

    miller: Vec<(u64, Vec<(usize, usize)>)>,

    final_exp_features: Vec<(String, u64)>,

    final_exp: Vec<(u64, Vec<(usize, usize)>)>
}

#[derive(Clone, Deserialize, Debug)]
pub(crate) struct Bls12PairingParams {
    #[serde(deserialize_with = "parse_hashmap_usize_u64_from_ints")]
    one_off: HashMap<usize, u64>,

    multiplier: u64,

    miller_features: Vec<(String, u64)>,

    miller: Vec<(u64, Vec<(usize, usize)>)>,

    final_exp_features: Vec<(String, u64)>,

    final_exp: Vec<(u64, Vec<(usize, usize)>)>
}

#[derive(Clone, Deserialize, Debug)]
pub(crate) struct BnPairingParams {
    #[serde(deserialize_with = "parse_hashmap_usize_u64_from_ints")]
    one_off: HashMap<usize, u64>,

    multiplier: u64,

    miller_features: Vec<(String, u64)>,

    miller: Vec<(u64, Vec<(usize, usize)>)>,

    final_exp_features: Vec<(String, u64)>,

    final_exp: Vec<(u64, Vec<(usize, usize)>)>
}

static MNT4_PARAMS_JSON: &'static str = include_str!("mnt4_model.json");
static MNT6_PARAMS_JSON: &'static str = include_str!("mnt6_model.json");
static BLS12_PARAMS_JSON: &'static str = include_str!("bls12_model.json");
static BN_PARAMS_JSON: &'static str = include_str!("bn_model.json");

pub(crate) static MNT4_PARAMS_INSTANCE: Lazy<MntPairingParams> = Lazy::new(|| {
    serde_json::from_str(MNT4_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) static MNT6_PARAMS_INSTANCE: Lazy<MntPairingParams> = Lazy::new(|| {
    serde_json::from_str(MNT6_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) static BLS12_PARAMS_INSTANCE: Lazy<Bls12PairingParams> = Lazy::new(|| {
    serde_json::from_str(BLS12_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) static BN_PARAMS_INSTANCE: Lazy<BnPairingParams> = Lazy::new(|| {
    serde_json::from_str(BN_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) fn meter_mnt_pairing(input: &[u8], params: &MntPairingParams, max_power: usize, ext_degree: usize) -> Result<u64, ApiError> {
    let (
        modulus, 
        order_len, 
        num_pairs, 
        (ate_loop_bits, ate_loop_hamming), 
        (exp_w0_bits, exp_w0_hamming),
        (exp_w1_bits, exp_w1_hamming),
        (num_g1_subgroup_checks, num_g2_subgroup_checks),
        _
    ) = parse_mnt_pairing_parameters(&input, ext_degree)?;

    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    // let order_limbs = num_units_for_group_order(&order)?;
    let order_limbs = num_units_for_group_order_length(order_len)?;

    let mut estimate = calculate_mnt_pairing_cost(
        modulus_limbs,
        order_limbs,
        num_pairs,
        (ate_loop_bits, ate_loop_hamming), 
        (exp_w0_bits, exp_w0_hamming),
        (exp_w1_bits, exp_w1_hamming),
        params,
        max_power
    )?;

    let g1_subgroup_check_cost_per_point = super::meter_arith::meter_multiplication(modulus_limbs, order_limbs, &*super::meter_arith::G1_MULTIPLICATION_PARAMS_INSTANCE, false)?;
    let g1_subgroup_checks_cost = g1_subgroup_check_cost_per_point.checked_mul(num_g1_subgroup_checks as u64).ok_or(ApiError::Overflow)?;

    estimate = estimate.checked_add(g1_subgroup_checks_cost).ok_or(ApiError::Overflow)?;

    let g2_subgroup_check_cost_per_point = match ext_degree {
        2 => {
            super::meter_arith::meter_multiplication(modulus_limbs, order_limbs, &*super::meter_arith::G2_EXT_2_MULTIPLICATION_PARAMS_INSTANCE, false)?
        },
        3 => {
            super::meter_arith::meter_multiplication(modulus_limbs, order_limbs, &*super::meter_arith::G2_EXT_3_MULTIPLICATION_PARAMS_INSTANCE, false)?
        },
        _ => {
            return Err(ApiError::InputError("Invalid extension degree for MNT4/6 pairing cost calculation".to_owned()));
        }
    };

    let g2_subgroup_checks_cost = g2_subgroup_check_cost_per_point.checked_mul(num_g2_subgroup_checks as u64).ok_or(ApiError::Overflow)?;

    estimate = estimate.checked_add(g2_subgroup_checks_cost).ok_or(ApiError::Overflow)?;

    Ok(estimate)
}

fn calculate_mnt_pairing_cost(
    modulus_limbs: usize,
    _order_limbs: usize,
    num_pairs: usize,
    (ate_loop_bits, ate_loop_hamming): (u64, u64), 
    (exp_w0_bits, exp_w0_hamming): (u64, u64),
    (exp_w1_bits, exp_w1_hamming): (u64, u64),
    params: &MntPairingParams, 
    max_power: usize

) -> Result<u64, ApiError> {
    const ATE_LOOP_BITS_INDEX: usize = 0;
    const ATE_LOOP_HAMMING_INDEX: usize = 1;
    const EXP_W0_LOOP_BITS_INDEX: usize = 2;
    const EXP_W0_HAMMING_INDEX: usize = 3;
    const EXP_W1_LOOP_BITS_INDEX: usize = 4;
    const EXP_W1_HAMMING_INDEX: usize = 5;

    debug_assert!(max_power == MNT4_MAX_MODULUS_POWER || max_power == MNT6_MAX_MODULUS_POWER);

    let one_off = *params.one_off.get(&modulus_limbs).ok_or(ApiError::MissingValue)?;

    let modulus_limbs_powers = make_powers(modulus_limbs as u64, max_power)?;
    let params_vector = vec![ate_loop_bits, ate_loop_hamming, exp_w0_bits, exp_w0_hamming, exp_w1_bits, exp_w1_hamming];

    let miller_cost = {
        let miller_params = vec![
            &params_vector[ATE_LOOP_BITS_INDEX..(ATE_LOOP_BITS_INDEX+1)], 
            &params_vector[ATE_LOOP_HAMMING_INDEX..(ATE_LOOP_HAMMING_INDEX+1)], 
            &modulus_limbs_powers[..] 
            ];
        let mut miller_cost = eval_model(&params.miller, &miller_params)?;
        miller_cost = miller_cost.checked_mul(num_pairs as u64).ok_or(ApiError::Overflow)?;

        miller_cost
    };

    println!("Miller cost = {}", miller_cost);

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

    println!("Final exp cost = {}", final_exp_cost);

    let mut result = one_off;
    result = result.checked_add(miller_cost).ok_or(ApiError::Overflow)?;
    result = result.checked_add(final_exp_cost).ok_or(ApiError::Overflow)?;
    result = result.checked_div(params.multiplier).ok_or(ApiError::Overflow)?;

    Ok(result)
}

pub(crate) fn meter_bls12_pairing(input: &[u8], params: &Bls12PairingParams, max_power: usize) -> Result<u64, ApiError> {
    let (
        modulus, 
        order_len, 
        num_pairs, 
        x,
        _,
        (num_g1_subgroup_checks, num_g2_subgroup_checks),
        _
    ) = parse_bls12_bn_pairing_parameters(&input, MAX_BLS12_X_BIT_LENGTH)?;

    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    // let order_limbs = num_units_for_group_order(&order)?;
    let order_limbs = num_units_for_group_order_length(order_len)?;

    let x_bits = x.bits();
    let x_hamming = calculate_hamming_weight(&x.as_ref());

    if x_hamming > MAX_BLS12_X_HAMMING {
        return Err(ApiError::InputError(format!("Hamming weight for scalar is too large, file {}, line {}", file!(), line!())));
    }

    let mut estimate = calculate_bls12_pairing_cost(
        modulus_limbs,
        order_limbs,
        num_pairs,
        (x_bits as u64, x_hamming as u64),
        params,
        max_power
    )?;

    let g1_subgroup_check_cost_per_point = super::meter_arith::meter_multiplication(modulus_limbs, order_limbs, &*super::meter_arith::G1_MULTIPLICATION_PARAMS_INSTANCE, false)?;
    let g1_subgroup_check_cost = g1_subgroup_check_cost_per_point.checked_mul(num_g1_subgroup_checks as u64).ok_or(ApiError::Overflow)?;

    estimate = estimate.checked_add(g1_subgroup_check_cost).ok_or(ApiError::Overflow)?;

    let g2_subgroup_check_cost_per_point = super::meter_arith::meter_multiplication(modulus_limbs, order_limbs, &*super::meter_arith::G2_EXT_2_MULTIPLICATION_PARAMS_INSTANCE, false)?;
    let g2_subgroup_check_cost = g2_subgroup_check_cost_per_point.checked_mul(num_g2_subgroup_checks as u64).ok_or(ApiError::Overflow)?;

    estimate = estimate.checked_add(g2_subgroup_check_cost).ok_or(ApiError::Overflow)?;

    Ok(estimate)
}


pub(crate) fn meter_bn_pairing(input: &[u8], params: &BnPairingParams, max_power: usize) -> Result<u64, ApiError> {
    let (
        modulus, 
        order_len, 
        num_pairs, 
        u,
        u_is_negative,
        (num_g1_subgroup_checks, num_g2_subgroup_checks),
        _
    ) = parse_bls12_bn_pairing_parameters(&input, MAX_BN_U_BIT_LENGTH)?;
    use crate::integers::MaxLoopParametersUint;

    let modulus_limbs = num_limbs_for_modulus(&modulus)?;
    // let order_limbs = num_units_for_group_order(&order)?;
    let order_limbs = num_units_for_group_order_length(order_len)?;

    let u_bits = u.bits();
    let u_hamming = calculate_hamming_weight(&u.as_ref());

    let two = MaxLoopParametersUint::from(2u64);
    let six = MaxLoopParametersUint::from(6u64);

    // we need only absolute value of 6u+2, so manually handle negative and positive U
    let six_u_plus_two = if u_is_negative {
        let six_u_plus_two = (six * u) - two;

        six_u_plus_two
    } else {
        let six_u_plus_two = (six * u) + two;

        six_u_plus_two
    };

    let six_u_plus_two_bits = six_u_plus_two.bits();

    let six_u_plus_two_hamming = calculate_hamming_weight(six_u_plus_two.as_ref());

    if six_u_plus_two_hamming > MAX_BN_SIX_U_PLUS_TWO_HAMMING {
        return Err(ApiError::InputError(format!("Hamming weight for scalar is too large, file {}, line {}", file!(), line!())));
    }

    let mut estimate = calculate_bn_pairing_cost(
        modulus_limbs,
        order_limbs,
        num_pairs,
        (six_u_plus_two_bits as u64, six_u_plus_two_hamming as u64),
        (u_bits as u64, u_hamming as u64),
        params,
        max_power
    )?;

    let g1_subgroup_check_cost_per_point = super::meter_arith::meter_multiplication(modulus_limbs, order_limbs, &*super::meter_arith::G1_MULTIPLICATION_PARAMS_INSTANCE, false)?;
    let g1_subgroup_check_cost = g1_subgroup_check_cost_per_point.checked_mul(num_g1_subgroup_checks as u64).ok_or(ApiError::Overflow)?;

    estimate = estimate.checked_add(g1_subgroup_check_cost).ok_or(ApiError::Overflow)?;

    let g2_subgroup_check_cost_per_point = super::meter_arith::meter_multiplication(modulus_limbs, order_limbs, &*super::meter_arith::G2_EXT_2_MULTIPLICATION_PARAMS_INSTANCE, false)?;
    let g2_subgroup_check_cost = g2_subgroup_check_cost_per_point.checked_mul(num_g2_subgroup_checks as u64).ok_or(ApiError::Overflow)?;

    estimate = estimate.checked_add(g2_subgroup_check_cost).ok_or(ApiError::Overflow)?;

    Ok(estimate)
}

fn calculate_bls12_pairing_cost(
    modulus_limbs: usize,
    _order_limbs: usize,
    num_pairs: usize,
    (x_bits, x_hamming): (u64, u64),
    params: &Bls12PairingParams, 
    max_power: usize

) -> Result<u64, ApiError> {
    const X_BITS_INDEX: usize = 0;
    const X_HAMMING_INDEX: usize = 1;

    debug_assert!(max_power == BLS12_MAX_MODULUS_POWER);

    let modulus_limbs_powers = make_powers(modulus_limbs as u64, max_power)?;
    let params_vector = vec![x_bits, x_hamming];

    let one_off = *params.one_off.get(&modulus_limbs).ok_or(ApiError::MissingValue)?;

    let miller_cost = {
        let miller_params = vec![
            &params_vector[X_BITS_INDEX..(X_BITS_INDEX+1)], 
            &params_vector[X_HAMMING_INDEX..(X_HAMMING_INDEX+1)], 
            &modulus_limbs_powers[..] 
            ];
        let mut miller_cost = eval_model(&params.miller, &miller_params)?;
        miller_cost = miller_cost.checked_mul(num_pairs as u64).ok_or(ApiError::Overflow)?;

        miller_cost
    };

    println!("Miller cost = {}", miller_cost);

    let final_exp_cost = {
        let final_exp_params = vec![
            &params_vector[X_BITS_INDEX..(X_BITS_INDEX+1)], 
            &params_vector[X_HAMMING_INDEX..(X_HAMMING_INDEX+1)], 
            &modulus_limbs_powers[..] 
            ];
        let final_exp_cost = eval_model(&params.final_exp, &final_exp_params)?;

        final_exp_cost
    };

    println!("Final exp cost = {}", final_exp_cost);

    let mut result = one_off;
    result = result.checked_add(miller_cost).ok_or(ApiError::Overflow)?;
    result = result.checked_add(final_exp_cost).ok_or(ApiError::Overflow)?;
    result = result.checked_div(params.multiplier).ok_or(ApiError::Overflow)?;

    Ok(result)
}

fn calculate_bn_pairing_cost(
    modulus_limbs: usize,
    _order_limbs: usize,
    num_pairs: usize,
    (six_u_plus_two_bits, six_u_plus_two_hamming): (u64, u64),
    (u_bits, u_hamming): (u64, u64),
    params: &BnPairingParams, 
    max_power: usize
) -> Result<u64, ApiError> {
    const U_BITS_INDEX: usize = 0;
    const U_HAMMING_INDEX: usize = 1;
    const SIX_U_PLUS_TWO_BITS_INDEX: usize = 2;
    const SIX_U_PLUS_TWO_HAMMING_INDEX: usize = 3;

    debug_assert!(max_power == BN_MAX_MODULUS_POWER);

    let modulus_limbs_powers = make_powers(modulus_limbs as u64, max_power)?;
    let params_vector = vec![u_bits, u_hamming, six_u_plus_two_bits, six_u_plus_two_hamming];

    let one_off = *params.one_off.get(&modulus_limbs).ok_or(ApiError::MissingValue)?;

    let miller_cost = {
        let miller_params = vec![
            &params_vector[SIX_U_PLUS_TWO_BITS_INDEX..(SIX_U_PLUS_TWO_BITS_INDEX+1)], 
            &params_vector[SIX_U_PLUS_TWO_HAMMING_INDEX..(SIX_U_PLUS_TWO_HAMMING_INDEX+1)], 
            &modulus_limbs_powers[..] 
            ];
        let mut miller_cost = eval_model(&params.miller, &miller_params)?;
        miller_cost = miller_cost.checked_mul(num_pairs as u64).ok_or(ApiError::Overflow)?;

        miller_cost
    };

    let final_exp_cost = {
        let final_exp_params = vec![
            &params_vector[U_BITS_INDEX..(U_BITS_INDEX+1)], 
            &params_vector[U_HAMMING_INDEX..(U_HAMMING_INDEX+1)], 
            &modulus_limbs_powers[..] 
            ];
        let final_exp_cost = eval_model(&params.final_exp, &final_exp_params)?;

        final_exp_cost
    };

    let mut result = one_off;
    result = result.checked_add(miller_cost).ok_or(ApiError::Overflow)?;
    result = result.checked_add(final_exp_cost).ok_or(ApiError::Overflow)?;
    result = result.checked_div(params.multiplier).ok_or(ApiError::Overflow)?;

    Ok(result)
}

fn eval_model(
    coeffs_variables_and_powers: &[(u64, Vec<(usize, usize)>)],
    variables: &[ &[u64] ]
) -> Result<u64, ApiError> {
    let mut final_result = 0u64;
    if coeffs_variables_and_powers.len() == 0 {
        return Err(ApiError::MissingValue);
    }
    let mut max_var_id = 0usize;
    for (_, var_and_power) in coeffs_variables_and_powers.iter() {
        for (variable, _) in var_and_power.iter() {
            if max_var_id < *variable {
                max_var_id = *variable;
            }
        }
    }

    if max_var_id + 1 != variables.len() {
        // println!("Max variable ID (zero enumerated) {} is missing: coeffs = {:?}, variables = {:?}", max_var_id, coeffs_variables_and_powers, variables);
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

        let t = &*super::BLS12_PARAMS_INSTANCE;
        println!("Params BLS12 = {:?}", t);

        let t = &*super::BN_PARAMS_INSTANCE;
        println!("Params BN = {:?}", t);
    }

    #[test]
    fn test_calculate_example_prices_ey_pendulum() {
        let ey_pendulum_gas_cost = super::calculate_mnt_pairing_cost(
            10, 
            5, 
            1, 
            (613, 292), 
            (613, 312), 
            (315, 157), 
            &*super::MNT6_PARAMS_INSTANCE, 
            6
        ).unwrap();

        println!("EY pendulum cost for 1 pair = {}", ey_pendulum_gas_cost);
    }

    #[test]
    fn test_calculate_example_prices_ey_sw6_bis() {
        use num_bigint::BigUint;
        use num_traits::Num;
        use crate::test::biguint_to_u64_vec;

        let ate = BigUint::from_str_radix("6891450384315732539396789682275657542479668912536150109513790160209623422243491736087683183289411687640864567753786354786735746151460237106615816805591765205438850184717968966109876910955248455129124731541829539482640472797610122", 10).unwrap();
        let w0 = BigUint::from_str_radix("2156695813352724974824326851054479880127610960548355747044807332080688727374737671308314095389122345740953981240668571898337613282699493372314698360451061276517306188376803619985090458895588556562724088277106828", 10).unwrap();
        let w1 = BigUint::from_str_radix("26642435879335816683987677701488073867751118270052650655942102502312977592501693353047140953112195348280268661194889", 10).unwrap();
    
        let ate = biguint_to_u64_vec(ate);
        let w0 = biguint_to_u64_vec(w0);
        let w1 = biguint_to_u64_vec(w1);

        let ate_bits = crate::public_interface::decode_utils::calculate_hamming_weight(&ate);
        let ate_hamming = crate::public_interface::decode_utils::calculate_hamming_weight(&ate);

        let w0_bits = crate::public_interface::decode_utils::calculate_hamming_weight(&w0);
        let w0_hamming = crate::public_interface::decode_utils::calculate_hamming_weight(&w0);

        let w1_bits = crate::public_interface::decode_utils::calculate_hamming_weight(&w1);
        let w1_hamming = crate::public_interface::decode_utils::calculate_hamming_weight(&w1);
        
        let gas_cost = super::calculate_mnt_pairing_cost(
            12, 
            6, 
            1, 
            (ate_bits as u64, ate_hamming as u64), 
            (w0_bits as u64, w0_hamming as u64), 
            (w1_bits as u64, w1_hamming as u64), 
            &*super::MNT6_PARAMS_INSTANCE, 
            6
        ).unwrap();

        println!("EY Sw6 bis cost for 1 pair = {}", gas_cost);
    }

    #[test]
    fn test_calculate_example_prices_bn254() {
        let u_hamming = crate::public_interface::decode_utils::calculate_hamming_weight(&[0x44e992b44a6909f1]);
        let six_u_plus_two_hamming = crate::public_interface::decode_utils::calculate_hamming_weight(&[0x9d797039be763ba8, 1]);
        let bn_4_pairs_cost = super::calculate_bn_pairing_cost(
            4, 
            4, 
            4, 
            (65, six_u_plus_two_hamming as u64), 
            (63, u_hamming as u64), 
            &*super::BN_PARAMS_INSTANCE, 
            6).unwrap();

        println!("BN254 for 4 pairs = {}", bn_4_pairs_cost);
        
    }

    #[test]
    fn test_calculate_example_prices_bls12_381() {
        let x_hamming = crate::public_interface::decode_utils::calculate_hamming_weight(&[0xd201000000010000]);
        let x_bits = 64 - 0xd201000000010000u64.leading_zeros();
        let bls12_4_pairs_cost = super::calculate_bls12_pairing_cost(
            6, 
            4, 
            4, 
            (x_bits as u64, x_hamming as u64), 
            &*super::BLS12_PARAMS_INSTANCE, 
            6).unwrap();

        println!("BN381 for 4 pairs = {}", bls12_4_pairs_cost);
        
    }

    #[test]
    fn test_print_example_prices_bls12_381() {
        let x_hamming = crate::public_interface::decode_utils::calculate_hamming_weight(&[0xd201000000010000]);
        let x_bits = 64 - 0xd201000000010000u64.leading_zeros();
        let bls12_1_pair_cost = super::calculate_bls12_pairing_cost(
            6, 
            4, 
            1, 
            (x_bits as u64, x_hamming as u64), 
            &*super::BLS12_PARAMS_INSTANCE, 
            6).unwrap();

        println!("BN381 for 1 pair = {}", bls12_1_pair_cost);
        
    }

    #[test]
    fn test_print_example_prices_bls12_377() {
        let x: u64 = 0x8508c00000000001;
        let x_hamming = crate::public_interface::decode_utils::calculate_hamming_weight(&[x]);
        let x_bits = 64 - x.leading_zeros();
        let bls12_1_pair_cost = super::calculate_bls12_pairing_cost(
            6, 
            4, 
            1, 
            (x_bits as u64, x_hamming as u64), 
            &*super::BLS12_PARAMS_INSTANCE, 
            6).unwrap();

        println!("BN377 for 1 pair = {}", bls12_1_pair_cost);
        
    }
}