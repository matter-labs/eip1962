use num_bigint::BigUint;
use num_traits::Num;

extern crate serde;
extern crate serde_json;

use serde::{Deserialize, Deserializer};

#[derive(Deserialize, Debug, Clone)]
pub(crate) struct JsonBls12PairingCurveParameters {
    #[serde(deserialize_with = "biguint_with_sign_from_hex_string")]
    pub non_residue: (BigUint, bool),

    #[serde(rename = "is_D_type")]
    #[serde(deserialize_with = "bool_from_string")]
    pub is_d_type: bool,

    #[serde(deserialize_with = "biguint_with_sign_from_hex_string")]
    pub quadratic_non_residue_0: (BigUint, bool),

    #[serde(deserialize_with = "biguint_with_sign_from_hex_string")]
    pub quadratic_non_residue_1: (BigUint, bool),

    #[serde(deserialize_with = "biguint_with_sign_from_hex_string")]
    pub x: (BigUint, bool),

    #[serde(deserialize_with = "biguint_from_hex_string")]
    pub q: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    pub r: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "A")]
    pub a: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "B")]
    pub b: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "A_twist_0")]
    pub a_twist_0: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "A_twist_1")]
    pub a_twist_1: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "B_twist_0")]
    pub b_twist_0: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "B_twist_1")]
    pub b_twist_1: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    pub g1_x: BigUint,
    #[serde(deserialize_with = "biguint_from_hex_string")]
    pub g1_y: BigUint,
    #[serde(deserialize_with = "biguint_from_hex_string")]
    pub g2_x_0: BigUint,
    #[serde(deserialize_with = "biguint_from_hex_string")]
    pub g2_x_1: BigUint,
    #[serde(deserialize_with = "biguint_from_hex_string")]
    pub g2_y_0: BigUint,
    #[serde(deserialize_with = "biguint_from_hex_string")]
    pub g2_y_1: BigUint,

    #[serde(rename = "g1_scalar_mult_test_vectors")]
    pub g1_mul_vectors: Vec<JsonG1PointScalarMultiplicationPair>,

    #[serde(rename = "g2_scalar_mult_test_vectors")]
    pub g2_mul_vectors: Vec<JsonG2PointScalarMultiplicationPair>,
}

#[derive(Deserialize, Debug, Clone)]
pub(crate) struct JsonG1PointScalarMultiplicationPair {
    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "a")]
    pub scalar: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "g_x")]
    pub base_x: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "g_y")]
    pub base_y: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "h_x")]
    pub result_x: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "h_y")]
    pub result_y: BigUint,
}

#[derive(Deserialize, Debug, Clone)]
pub(crate) struct JsonG2PointScalarMultiplicationPair {
    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "a")]
    pub scalar: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "g_x_0")]
    pub base_x_0: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "g_x_1")]
    pub base_x_1: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "g_y_0")]
    pub base_y_0: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "g_y_1")]
    pub base_y_1: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "h_x_0")]
    pub result_x_0: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "h_x_1")]
    pub result_x_1: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "h_y_0")]
    pub result_y_0: BigUint,

    #[serde(deserialize_with = "biguint_from_hex_string")]
    #[serde(rename = "h_y_1")]
    pub result_y_1: BigUint,
}

fn bool_from_string<'de, D>(deserializer: D) -> Result<bool, D::Error>
where
    D: Deserializer<'de>,
{
    match String::deserialize(deserializer)?.as_ref() {
        "True" => Ok(true),
        "False" => Ok(false),
        other => Err(serde::de::Error::invalid_value(
            serde::de::Unexpected::Str(other),
            &"True or False",
        )),
    }
}

fn biguint_from_hex_string<'de, D>(deserializer: D) -> Result<BigUint, D::Error>
where
    D: Deserializer<'de>,
{
    let string_value = strip_0x(&String::deserialize(deserializer)?);
    let value = BigUint::from_str_radix(&string_value, 16).map_err(|_| {
        serde::de::Error::invalid_value(
            serde::de::Unexpected::Str(&string_value),
            &"Not valid hex number",
        )
    })?;

    Ok(value)
}

fn biguint_with_sign_from_hex_string<'de, D>(deserializer: D) -> Result<(BigUint, bool), D::Error>
where
    D: Deserializer<'de>,
{
    let (string_value, is_positive) = strip_0x_and_get_sign(&String::deserialize(deserializer)?);
    let value = BigUint::from_str_radix(&string_value, 16).map_err(|_| {
        serde::de::Error::invalid_value(
            serde::de::Unexpected::Str(&string_value),
            &"Not valid hex number",
        )
    })?;

    Ok((value, is_positive))
}

fn strip_0x(string: &str) -> String {
    let string = string.trim();
    let mut string = string.to_ascii_lowercase().as_bytes().to_vec();
    if string.len() > 2 && string[0..1] == b"0"[..] && string[1..2] == b"x"[..] {
        string = string[2..].to_vec();
    }

    if string.len() > 1 && string[string.len() - 1] == b"l"[0] {
        string = string[..(string.len() - 1)].to_vec();
    }
    
    std::string::String::from_utf8(string).unwrap()
}

fn strip_0x_and_get_sign(string: &str) -> (String, bool) {
    let string = string.trim();
    let mut string = string.to_ascii_lowercase().as_bytes().to_vec();
    let mut positive = true;
    if string.len() > 1 && string[0..1] == b"-"[..] {
        string = string[1..].to_vec();
        positive = false;
    }
    if string.len() > 1 && string[0..1] == b"+"[..] {
        string = string[1..].to_vec();
    }
    if string.len() > 2 && string[0..1] == b"0"[..] && string[1..2] == b"x"[..] {
        string = string[2..].to_vec();
    }
    
    (std::string::String::from_utf8(string).unwrap(), positive)
}

fn strip_0x_and_pad(string: &str) -> String {
    let string = string.trim();
    let mut string = string.to_ascii_lowercase().as_bytes().to_vec();
    if string.len() > 2 && string[0..1] == b"0"[..] && string[1..2] == b"x"[..] {
        let mut string = string[2..].to_vec();
        if string.len() % 2 == 1 {
            string = {
                let mut res = "0".as_bytes().to_vec();
                res.extend(string.into_iter());

                res
            };
        }
        return std::string::String::from_utf8(string).unwrap();
    }
    if string.len() % 2 == 1 {
        string = {
            let mut res = "0".as_bytes().to_vec();
            res.extend(string.into_iter());

            res
        };
    }

    std::string::String::from_utf8(string).unwrap()
}

pub(crate) fn read_dir_and_grab_curves(dir_path: &str) -> Vec<(JsonBls12PairingCurveParameters, String)> {
    use std::io::Read;
    use std::fs::{self};
    use std::path::Path;
    use std::fs::File;

    let dir = Path::new(dir_path);
    assert!(dir.is_dir());
    let mut results = vec![];
    for entry in fs::read_dir(dir).expect("must read the directory") {
        let entry = entry.expect("directory should contain files");
        let path = entry.path();
        if path.is_dir() {
            continue
        } else {
            let extension = path.extension();
            if extension.is_none() {
                continue
            }
            let extension = extension.unwrap();
            if extension != "curve" {
                continue
            }
        }
        let mut buffer = Vec::new();
        let file_name = path.file_name().unwrap().to_str().unwrap().to_owned();
        let mut f = File::open(path).expect("must open file");
        f.read_to_end(&mut buffer).expect("must read bytes from file");
        let c: JsonBls12PairingCurveParameters = serde_json::from_slice(&buffer[..]).expect("must deserialize");
        results.push((c, file_name));
    }
    
    results
}

pub(crate) fn pad_for_len_be(input: Vec<u8>, len: usize) -> Vec<u8> {
    if input.len() < len {
        let mut res = input;
        res.reverse();
        res.resize(len, 0u8);
        res.reverse();
        return res;
    }

    input
}

pub(crate) fn prepend_0x(input: &str) -> String {
    format!("0x{}", input)
}