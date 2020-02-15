use serde::{Deserialize};
use std::collections::HashMap;
use crate::errors::ApiError;

use once_cell::sync::Lazy;
use serde_json;

use super::parsers::*;

pub(crate) trait ArithmeticAdditionParams {
    fn params(&self) -> &HashMap<usize, u64>;
}

pub(crate) trait ArithmeticMultiplicationParams {
    fn params(&self) -> (&HashMap<usize, u64>, &HashMap<usize, u64>);
}

pub(crate) trait ArithmeticMultiexpParams {
    fn params(&self) -> (u64, (usize, u64), &HashMap<usize, u64>);
}

#[derive(Clone, Deserialize, Debug)]
pub(crate) struct G1G2AdditionParams {
    #[serde(deserialize_with = "parse_hashmap_usize_u64_from_ints")]
    #[serde(rename = "price")]
    lookup_parameters: HashMap<usize, u64>
}

impl ArithmeticAdditionParams for G1G2AdditionParams {
    fn params(&self) -> &HashMap<usize, u64> {
        &self.lookup_parameters
    }
}

#[derive(Clone, Deserialize, Debug)]
pub(crate) struct G1G2MultiplicationParams {
    #[serde(deserialize_with = "parse_hashmap_usize_u64_from_ints")]
    #[serde(rename = "base")]
    base: HashMap<usize, u64>,

    #[serde(deserialize_with = "parse_hashmap_usize_u64_from_ints")]
    #[serde(rename = "per_limb")]
    per_limb: HashMap<usize, u64>
}

impl ArithmeticMultiplicationParams for G1G2MultiplicationParams {
    fn params(&self) -> (&HashMap<usize, u64>, &HashMap<usize, u64>) {
        (&self.base, &self.per_limb)
    }
}

#[derive(Clone, Deserialize, Debug)]
pub(crate) struct G1G2MultiexpParams {
    // #[serde(deserialize_with = "parse_usize")]
    #[serde(rename = "max_pairs")]
    max_pairs: usize,

    // #[serde(deserialize_with = "parse_u64")]
    #[serde(rename = "max_discount")]
    max_discount: u64,

    // #[serde(deserialize_with = "parse_u64")]
    #[serde(rename = "discount_multiplier")]
    discount_multiplier: u64,

    // #[serde(deserialize_with = "parse_hashmap_usize_u64")]
    #[serde(deserialize_with = "parse_hashmap_usize_u64_from_ints")]
    #[serde(rename = "discounts")]
    discounts: HashMap<usize, u64>
}

impl ArithmeticMultiexpParams for G1G2MultiexpParams {
    fn params(&self) -> (u64, (usize, u64), &HashMap<usize, u64>) {
        (self.discount_multiplier, (self.max_pairs, self.max_discount), &self.discounts)
    }
}

static G1_ADDITION_PARAMS_JSON: &'static str = include_str!("g1_addition.json");
static G2_EXT_2_ADDITION_PARAMS_JSON: &'static str = include_str!("g2_addition_ext2.json");
static G2_EXT_3_ADDITION_PARAMS_JSON: &'static str = include_str!("g2_addition_ext3.json");
static G1_MULTIPLICATION_PARAMS_JSON: &'static str = include_str!("g1_multiplication.json");
static G2_EXT_2_MULTIPLICATION_PARAMS_JSON: &'static str = include_str!("g2_multiplication_ext2.json");
static G2_EXT_3_MULTIPLICATION_PARAMS_JSON: &'static str = include_str!("g2_multiplication_ext3.json");
static MULTIEXP_PARAMS_JSON: &'static str = include_str!("multiexp_discounts.json");

pub(crate) static G1_ADDITION_PARAMS_INSTANCE: Lazy<G1G2AdditionParams> = Lazy::new(|| {
    serde_json::from_str(G1_ADDITION_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) static G2_EXT_2_ADDITION_PARAMS_INSTANCE: Lazy<G1G2AdditionParams> = Lazy::new(|| {
    serde_json::from_str(G2_EXT_2_ADDITION_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) static G2_EXT_3_ADDITION_PARAMS_INSTANCE: Lazy<G1G2AdditionParams> = Lazy::new(|| {
    serde_json::from_str(G2_EXT_3_ADDITION_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) static G1_MULTIPLICATION_PARAMS_INSTANCE: Lazy<G1G2MultiplicationParams> = Lazy::new(|| {
    serde_json::from_str(G1_MULTIPLICATION_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) static G2_EXT_2_MULTIPLICATION_PARAMS_INSTANCE: Lazy<G1G2MultiplicationParams> = Lazy::new(|| {
    serde_json::from_str(G2_EXT_2_MULTIPLICATION_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) static G2_EXT_3_MULTIPLICATION_PARAMS_INSTANCE: Lazy<G1G2MultiplicationParams> = Lazy::new(|| {
    serde_json::from_str(G2_EXT_3_MULTIPLICATION_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) static MULTIEXP_PARAMS_INSTANCE: Lazy<G1G2MultiexpParams> = Lazy::new(|| {
    serde_json::from_str(MULTIEXP_PARAMS_JSON).expect("must deserialize parameters")
});

pub(crate) fn meter_addition<P: ArithmeticAdditionParams>(modulus_limbs: usize, parameters: &P) -> Result<u64, ApiError> {
    let found = *parameters.params().get(&modulus_limbs).ok_or(ApiError::MissingValue)?;

    return Ok(found)
}

pub(crate) fn meter_multiplication<P: ArithmeticMultiplicationParams>(
    modulus_limbs: usize, 
    group_limbs: usize, 
    parameters: &P,
    include_base: bool
) -> Result<u64, ApiError> {
    let (one_shot_params, per_limb_params) =  parameters.params();
    let one_shot = *one_shot_params.get(&modulus_limbs).ok_or(ApiError::MissingValue)?;
    let per_limb = *per_limb_params.get(&modulus_limbs).ok_or(ApiError::MissingValue)?;

    let mut result = per_limb.checked_mul(group_limbs as u64).ok_or(ApiError::Overflow)?;
    if include_base {
        result = result.checked_add(one_shot).ok_or(ApiError::Overflow)?;
    }

    return Ok(result)
}

pub(crate) fn meter_multiexp<P: ArithmeticMultiplicationParams, M: ArithmeticMultiexpParams>(
    modulus_limbs: usize, 
    group_limbs: usize, 
    num_pairs: usize, 
    parameters: &P, 
    multiexp_discounts: &M
) -> Result<u64, ApiError> {
    let per_pair = meter_multiplication(modulus_limbs, group_limbs, parameters, true)?;

    let (discount_multiplier, (max_pairs, max_discount), discount_lookup) = multiexp_discounts.params();

    let discount = if num_pairs > max_pairs {
        max_discount
    } else {
        *discount_lookup.get(&num_pairs).ok_or(ApiError::MissingValue)?
    };

    let mut result = per_pair.checked_mul(num_pairs as u64).ok_or(ApiError::Overflow)?;
    result = result.checked_mul(discount).ok_or(ApiError::Overflow)?;
    result = result.checked_div(discount_multiplier).ok_or(ApiError::Overflow)?;

    Ok(result)
}

#[cfg(test)]
mod test {
    #[test]
    fn test_deserialization() {
        let t = &*super::G1_ADDITION_PARAMS_INSTANCE;
        println!("Params G1 add = {:?}", t);

        let t = &*super::G2_EXT_2_ADDITION_PARAMS_INSTANCE;
        println!("Params G2 add ext 2= {:?}", t);

        let t = &*super::G2_EXT_3_ADDITION_PARAMS_INSTANCE;
        println!("Params G2 add ext 3 = {:?}", t);

        let t = &*super::G1_MULTIPLICATION_PARAMS_INSTANCE;
        println!("Params G1 mul = {:?}", t);

        let t = &*super::G2_EXT_2_MULTIPLICATION_PARAMS_INSTANCE;
        println!("Params G2 mul ext 2 = {:?}", t);

        let t = &*super::G2_EXT_3_MULTIPLICATION_PARAMS_INSTANCE;
        println!("Params G2 mul ext 3 = {:?}", t);

        let t = &*super::MULTIEXP_PARAMS_INSTANCE;
        println!("Multiexp discounts = {:?}", t);
    }

    #[test]
    fn test_calculate_example_arithmetic_prices_bn254() {
        let addition_price = super::meter_addition(4, &*super::G1_ADDITION_PARAMS_INSTANCE).unwrap();
        let mul_price = super::meter_multiplication(4, 4, &*super::G1_MULTIPLICATION_PARAMS_INSTANCE, true).unwrap();

        println!("BN254 addition price = {}", addition_price);
        println!("BN254 multiplication price = {}", mul_price);
        
    }

    #[test]
    fn test_calculate_example_arithmetic_prices_mnt4_753() {
        let mul_price = super::meter_multiplication(12, 12, &*super::G2_EXT_2_MULTIPLICATION_PARAMS_INSTANCE, true).unwrap();

        println!("MNT4-753 G2 multiplication price = {}", mul_price); 
    }
}