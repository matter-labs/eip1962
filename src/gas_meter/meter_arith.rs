use serde::{Deserialize, Deserializer};
use std::collections::HashMap;
use crate::errors::ApiError;

use once_cell::sync::Lazy;
use serde_json;

pub(crate) trait ArithmeticAdditionParams {
    fn params(&self) -> &[(usize, u64)];
}

pub(crate) trait ArithmeticMultiplicationParams {
    fn params(&self) -> (&[(usize, u64)], &[(usize, u64)]);
}

pub(crate) trait ArithmeticMultiexpParams {
    fn params(&self) -> (u64, (usize, u64), &HashMap<usize, u64>);
}


#[derive(Clone, Deserialize, Debug)]
pub(crate) struct G1G2AdditionParams {
    // #[serde(deserialize_with = "parse_tuple_usize_u64")]
    #[serde(rename = "price")]
    lookup_parameters: Vec<(usize, u64)>
}

impl ArithmeticAdditionParams for G1G2AdditionParams {
    fn params(&self) -> &[(usize, u64)] {
        &self.lookup_parameters
    }
}

#[derive(Clone, Deserialize, Debug)]
pub(crate) struct G1G2MultiplicationParams {
    // #[serde(deserialize_with = "parse_tuple_usize_u64")]
    #[serde(rename = "base")]
    base: Vec<(usize, u64)>,

    // #[serde(deserialize_with = "parse_tuple_usize_u64")]
    #[serde(rename = "per_limb")]
    per_limb: Vec<(usize, u64)>
}

impl ArithmeticMultiplicationParams for G1G2MultiplicationParams {
    fn params(&self) -> (&[(usize, u64)], &[(usize, u64)]) {
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

pub(crate) fn parse_tuple_usize_u64<'de, D>(deserializer: D) -> Result<Vec<(usize, u64)>, D::Error>
where
    D: Deserializer<'de>,
{
    use serde_json::Value;
    use serde::de::{Visitor, SeqAccess};

    struct MyVisitor;

    impl<'de> Visitor<'de> for MyVisitor
    {
        type Value = Vec<(usize, u64)>;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("a nonempty sequence of numbers")
        }

        fn visit_seq<S>(self, mut seq: S) -> Result<Self::Value, S::Error>
        where
            S: SeqAccess<'de>,
        {
            let mut results = Vec::with_capacity(seq.size_hint().unwrap_or(100));
            while let Some(value) = seq.next_element::<[Value; 2]>()? {
                let first = value[0].as_str().expect("is a string").parse::<usize>().expect(&format!("should be an integer {:?}", value[0]));
                let second = value[1].as_str().expect("is a string").parse::<u64>().expect(&format!("should be an integer {:?}", value[1]));
                results.push((first, second))
            }

            Ok(results)
        }
    }

    let visitor = MyVisitor;
    let result = deserializer.deserialize_seq(visitor)?;

    Ok(result)
}

pub(crate) fn parse_hashmap_usize_u64<'de, D>(deserializer: D) -> Result<HashMap<usize, u64>, D::Error>
where
    D: Deserializer<'de>,
{
    use serde_json::Value;
    use serde::de::{Visitor, SeqAccess};

    struct MyVisitor;

    impl<'de> Visitor<'de> for MyVisitor
    {
        type Value = HashMap<usize, u64>;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("a nonempty sequence of numbers")
        }

        fn visit_seq<S>(self, mut seq: S) -> Result<Self::Value, S::Error>
        where
            S: SeqAccess<'de>,
        {
            let mut results = HashMap::with_capacity(seq.size_hint().unwrap_or(100));
            while let Some(value) = seq.next_element::<[Value; 2]>()? {
                let first = value[0].as_str().expect("is a string").parse::<usize>().expect(&format!("should be an integer {:?}", value[0]));
                let second = value[1].as_str().expect("is a string").parse::<u64>().expect(&format!("should be an integer {:?}", value[1]));
                results.insert(first, second);
            }

            Ok(results)
        }
    }

    let visitor = MyVisitor;
    let result = deserializer.deserialize_seq(visitor)?;

    Ok(result)
}


pub(crate) fn parse_hashmap_usize_u64_from_ints<'de, D>(deserializer: D) -> Result<HashMap<usize, u64>, D::Error>
where
    D: Deserializer<'de>,
{
    use serde_json::Value;
    use serde::de::{Visitor, SeqAccess};

    struct MyVisitor;

    impl<'de> Visitor<'de> for MyVisitor
    {
        type Value = HashMap<usize, u64>;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            formatter.write_str("a nonempty sequence of numbers")
        }

        fn visit_seq<S>(self, mut seq: S) -> Result<Self::Value, S::Error>
        where
            S: SeqAccess<'de>,
        {
            let mut results = HashMap::with_capacity(seq.size_hint().unwrap_or(100));
            while let Some(value) = seq.next_element::<[Value; 2]>()? {
                let first = value[0].as_u64().expect(&format!("should be an integer {:?}", value[0])) as usize;
                let second = value[1].as_u64().expect(&format!("should be an integer {:?}", value[1]));
                results.insert(first, second);
            }

            Ok(results)
        }
    }

    let visitor = MyVisitor;
    let result = deserializer.deserialize_seq(visitor)?;

    Ok(result)
}

pub(crate) fn parse_usize<'de, D>(deserializer: D) -> Result<usize, D::Error>
where
    D: Deserializer<'de>,
{
    let string_value = &String::deserialize(deserializer)?;
    let value = string_value.parse::<usize>().expect("is integer");

    Ok(value)
}

pub(crate) fn parse_u64<'de, D>(deserializer: D) -> Result<u64, D::Error>
where
    D: Deserializer<'de>,
{
    let string_value = &String::deserialize(deserializer)?;
    let value = string_value.parse::<u64>().expect("is integer");

    Ok(value)
}

pub(crate) fn meter_addition<P: ArithmeticAdditionParams>(modulus_limbs: usize, parameters: &P) -> Result<u64, ApiError> {
    let mut res: Vec<_> = parameters.params().iter().filter(|(limbs, _)| *limbs == modulus_limbs).collect();
    if res.len() != 1 {
        return Err(ApiError::UnknownParameter(format!("Unknown number of limbs = {}", modulus_limbs)));
    }

    let found = res.pop().expect("result exists").1;

    return Ok(found)
}

pub(crate) fn meter_multiplication<P: ArithmeticMultiplicationParams>(modulus_limbs: usize, group_limbs: usize, parameters: &P) -> Result<u64, ApiError> {
    let (one_shot_params, per_limb_params) =  parameters.params();
    let mut one_shot_res: Vec<_> = one_shot_params.iter().filter(|(limbs, _)| *limbs == modulus_limbs).collect();
    if one_shot_res.len() != 1 {
        return Err(ApiError::UnknownParameter(format!("Unknown number of limbs = {}", modulus_limbs)));
    }

    let mut per_limb_coeff: Vec<_> = per_limb_params.iter().filter(|(limbs, _)| *limbs == modulus_limbs).collect();
    if per_limb_coeff.len() != 1 {
        return Err(ApiError::UnknownParameter(format!("Unknown number of limbs = {}", modulus_limbs)));
    }

    let one_shot = one_shot_res.pop().expect("result exists").1;
    let per_limb = per_limb_coeff.pop().expect("result exists").1;

    let mut result = per_limb.checked_mul(group_limbs as u64).ok_or(ApiError::Overflow)?;
    result = result.checked_add(one_shot).ok_or(ApiError::Overflow)?;

    return Ok(result)
}

pub(crate) fn meter_multiexp<P: ArithmeticMultiplicationParams, M: ArithmeticMultiexpParams>(
    modulus_limbs: usize, 
    group_limbs: usize, 
    num_pairs: usize, 
    parameters: &P, 
    multiexp_discounts: &M
) -> Result<u64, ApiError> {
    let per_pair = meter_multiplication(modulus_limbs, group_limbs, parameters)?;

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
        let mul_price = super::meter_multiplication(4, 4, &*super::G1_MULTIPLICATION_PARAMS_INSTANCE).unwrap();

        println!("BN254 addition price = {}", addition_price);
        println!("BN254 multiplication price = {}", mul_price);
        
    }
}