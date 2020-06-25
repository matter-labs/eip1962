pub struct EIP196Executor;

use crate::engines::bn254::*;
use crate::public_interface::ApiError;

pub const SCALAR_BYTE_LENGTH: usize = 32;

pub const SERIALIZED_FP_BYTE_LENGTH: usize = 32;
pub const SERIALIZED_G1_POINT_BYTE_LENGTH: usize = SERIALIZED_FP_BYTE_LENGTH * 2;

pub const SERIALIZED_FP2_BYTE_LENGTH: usize = SERIALIZED_FP_BYTE_LENGTH * 2;
pub const SERIALIZED_G2_POINT_BYTE_LENGTH: usize = SERIALIZED_FP2_BYTE_LENGTH * 2;

pub const SERIALIZED_PAIRING_RESULT_BYTE_LENGTH: usize = 32;

use crate::public_interface::decode_g1;
use crate::public_interface::decode_g2;

use crate::weierstrass::Group;
use crate::pairings::PairingEngine;

#[cfg(feature = "eip_196_c_api")]
pub mod c_api;

fn pairing_result_false() -> [u8; SERIALIZED_PAIRING_RESULT_BYTE_LENGTH] {
    [0u8; SERIALIZED_PAIRING_RESULT_BYTE_LENGTH]
}

fn pairing_result_true() -> [u8; SERIALIZED_PAIRING_RESULT_BYTE_LENGTH] {
    let mut res = [0u8; SERIALIZED_PAIRING_RESULT_BYTE_LENGTH];
    res[31] = 1u8;

    res
}

impl EIP196Executor {
    pub fn add<'a>(input: &'a [u8]) -> Result<[u8; SERIALIZED_G1_POINT_BYTE_LENGTH], ApiError> {
        if input.len() != SERIALIZED_G1_POINT_BYTE_LENGTH * 2 {
            return Err(ApiError::InputError("invalid input length for G1 addition".to_owned()));
        }

        let (mut p_0, rest) = decode_g1::decode_g1_point_from_xy_oversized(input, SERIALIZED_FP_BYTE_LENGTH, &*BN254_G1_CURVE)?;
        let (p_1, _) = decode_g1::decode_g1_point_from_xy_oversized(rest, SERIALIZED_FP_BYTE_LENGTH, &*BN254_G1_CURVE)?;

        if !p_0.is_on_curve() {
            if !crate::features::in_fuzzing_or_gas_metering() {
                return Err(ApiError::InputError(format!("Point 0 is not on curve, file {}, line {}", file!(), line!())));
            }
        }
        if !p_1.is_on_curve() {
            if !crate::features::in_fuzzing_or_gas_metering() {
                return Err(ApiError::InputError(format!("Point 1 is not on curve, file {}, line {}", file!(), line!())));
            }
        }

        p_0.add_assign(&p_1);

        let mut output = [0u8; SERIALIZED_G1_POINT_BYTE_LENGTH];

        let as_vec = decode_g1::serialize_g1_point(SERIALIZED_FP_BYTE_LENGTH, &p_0)?;

        output.copy_from_slice(&as_vec[..]);

        Ok(output)
    }

    pub fn mul<'a>(input: &'a [u8]) -> Result<[u8; SERIALIZED_G1_POINT_BYTE_LENGTH], ApiError> {
        if input.len() != SERIALIZED_G1_POINT_BYTE_LENGTH + SCALAR_BYTE_LENGTH {
            return Err(ApiError::InputError("invalid input length for G1 multiplication".to_owned()));
        }

        let (p_0, rest) = decode_g1::decode_g1_point_from_xy_oversized(input, SERIALIZED_FP_BYTE_LENGTH, &*BN254_G1_CURVE)?;
        let (scalar, _) = decode_g1::decode_scalar_representation(rest, SCALAR_BYTE_LENGTH)?;

        if !p_0.is_on_curve() {
            if !crate::features::in_fuzzing_or_gas_metering() {
                return Err(ApiError::InputError(format!("Point is not on curve, file {}, line {}", file!(), line!())));
            }
        }

        let p = p_0.mul(&scalar);

        let mut output = [0u8; SERIALIZED_G1_POINT_BYTE_LENGTH];

        let as_vec = decode_g1::serialize_g1_point(SERIALIZED_FP_BYTE_LENGTH, &p)?;

        output.copy_from_slice(&as_vec[..]);

        Ok(output)
    }

    pub fn pair<'a>(input: &'a [u8]) -> Result<[u8; SERIALIZED_PAIRING_RESULT_BYTE_LENGTH], ApiError> {
        if input.len() % (SERIALIZED_G2_POINT_BYTE_LENGTH + SERIALIZED_G1_POINT_BYTE_LENGTH) != 0 {
            return Err(ApiError::InputError("invalid input length for pairing".to_owned()));
        }
        let num_pairs = input.len() / (SERIALIZED_G2_POINT_BYTE_LENGTH + SERIALIZED_G1_POINT_BYTE_LENGTH);

        if num_pairs == 0 {
            return Ok(pairing_result_true());
            // return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
        }

        let mut global_rest = input;

        let mut g1_points = Vec::with_capacity(num_pairs);
        let mut g2_points = Vec::with_capacity(num_pairs);

        for _ in 0..num_pairs {
            let (g1, rest) = decode_g1::decode_g1_point_from_xy_oversized(global_rest, SERIALIZED_FP_BYTE_LENGTH, &*BN254_G1_CURVE)?;

            // g2 encoding in EIP 196/197 is non-standard: Fp2 element c0 + v*c1 where v is non-residue is
            // encoded as (c1, c0) instead of usual (c0, c1)
            let (g2, rest) = {
                use crate::public_interface::decode_utils::split;

                let (g2_encoding_bytes, rest) = split(rest, SERIALIZED_G2_POINT_BYTE_LENGTH, "not enough bytes to read G2 point")?;
                let mut swapped_encoding = [0u8; SERIALIZED_G2_POINT_BYTE_LENGTH];

                // swap for x coordinate
                (&mut swapped_encoding[0..SERIALIZED_FP_BYTE_LENGTH]).copy_from_slice(&g2_encoding_bytes[SERIALIZED_FP_BYTE_LENGTH..(SERIALIZED_FP_BYTE_LENGTH*2)]);
                (&mut swapped_encoding[SERIALIZED_FP_BYTE_LENGTH..(SERIALIZED_FP_BYTE_LENGTH*2)]).copy_from_slice(&g2_encoding_bytes[0..SERIALIZED_FP_BYTE_LENGTH]);

                // swap for y coordinate
                (&mut swapped_encoding[(SERIALIZED_FP_BYTE_LENGTH*2)..(SERIALIZED_FP_BYTE_LENGTH*3)]).copy_from_slice(&g2_encoding_bytes[(SERIALIZED_FP_BYTE_LENGTH*3)..(SERIALIZED_FP_BYTE_LENGTH*4)]);
                (&mut swapped_encoding[(SERIALIZED_FP_BYTE_LENGTH*3)..(SERIALIZED_FP_BYTE_LENGTH*4)]).copy_from_slice(&g2_encoding_bytes[(SERIALIZED_FP_BYTE_LENGTH*2)..(SERIALIZED_FP_BYTE_LENGTH*3)]);


                let (g2, _) = decode_g2::decode_g2_point_from_xy_in_fp2_oversized(&swapped_encoding[..], SERIALIZED_FP_BYTE_LENGTH, &*BN254_G2_CURVE)?;

                (g2, rest)
            };
            

            global_rest = rest;

            if !g1.is_on_curve() {
                if !crate::features::in_fuzzing_or_gas_metering() {
                    return Err(ApiError::InputError("G1 point is not on curve".to_owned()));
                }
            }

            if !g2.is_on_curve() {
                if !crate::features::in_fuzzing_or_gas_metering() {
                    return Err(ApiError::InputError("G2 point is not on curve".to_owned()));
                }
            }
            // "fast" subgroup checks using empirical data
            if g1.wnaf_mul_with_window_size(&BN254_SUBGROUP_ORDER[..], 5).is_zero() == false {
                if !crate::features::in_fuzzing_or_gas_metering() {
                    return Err(ApiError::InputError("G1 point is not in the expected subgroup".to_owned()));
                }
            }

            if g2.wnaf_mul_with_window_size(&BN254_SUBGROUP_ORDER[..], 5).is_zero() == false {
                if !crate::features::in_fuzzing_or_gas_metering() {
                    return Err(ApiError::InputError("G2 point is not in the expected subgroup".to_owned()));
                }
            }

            if !g1.is_zero() && !g2.is_zero() {
                g1_points.push(g1);
                g2_points.push(g2);
            }
        }

        debug_assert!(g1_points.len() == g2_points.len());

        if g1_points.len() == 0 {
            return Ok(pairing_result_true());
        }

        let engine = &*BN254_PAIRING_ENGINE;

        let pairing_result = engine.pair(&g1_points, &g2_points);

        if pairing_result.is_none() {
            return Err(ApiError::UnknownParameter("Pairing engine returned no value".to_owned()));
        }

        use crate::extension_towers::fp12_as_2_over3_over_2::Fp12;
        use crate::traits::ZeroAndOne;

        let one_fp12 = Fp12::one(&*BN254_EXT12_FIELD);
        let pairing_result = pairing_result.unwrap();
        let result = if pairing_result == one_fp12 {
            pairing_result_true()
        } else {
            pairing_result_false()
        };

        Ok(result)
    }
}