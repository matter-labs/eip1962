pub struct EIP2537Executor;

use crate::engines::bls12_381::{self, mapping};
use crate::public_interface::ApiError;

pub const SCALAR_BYTE_LENGTH: usize = 32;

pub const SERIALIZED_FP_BYTE_LENGTH: usize = 64;
pub const SERIALIZED_G1_POINT_BYTE_LENGTH: usize = SERIALIZED_FP_BYTE_LENGTH * 2;

pub const SERIALIZED_FP2_BYTE_LENGTH: usize = SERIALIZED_FP_BYTE_LENGTH * 2;
pub const SERIALIZED_G2_POINT_BYTE_LENGTH: usize = SERIALIZED_FP2_BYTE_LENGTH * 2;

pub const SERIALIZED_PAIRING_RESULT_BYTE_LENGTH: usize = 32;

use crate::public_interface::decode_fp;
use crate::public_interface::decode_g1;
use crate::public_interface::decode_g2;

use crate::weierstrass::Group;
use crate::multiexp::peppinger;
use crate::pairings::PairingEngine;

#[cfg(feature = "eip_2357_c_api")]
mod c_api;
#[cfg(feature = "eip_2357_c_api")]
pub use self::c_api::{c_perform_operation};

fn pairing_result_false() -> [u8; SERIALIZED_PAIRING_RESULT_BYTE_LENGTH] {
    [0u8; SERIALIZED_PAIRING_RESULT_BYTE_LENGTH]
}

fn pairing_result_true() -> [u8; SERIALIZED_PAIRING_RESULT_BYTE_LENGTH] {
    let mut res = [0u8; SERIALIZED_PAIRING_RESULT_BYTE_LENGTH];
    res[31] = 1u8;

    res
}

impl EIP2537Executor {
    pub fn g1_add<'a>(input: &'a [u8]) -> Result<[u8; SERIALIZED_G1_POINT_BYTE_LENGTH], ApiError> {
        if input.len() != SERIALIZED_G1_POINT_BYTE_LENGTH * 2 {
            return Err(ApiError::InputError("invalid input length for G1 addition".to_owned()));
        }

        let (mut p_0, rest) = decode_g1::decode_g1_point_from_xy_oversized(input, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_G1_CURVE)?;
        let (p_1, _) = decode_g1::decode_g1_point_from_xy_oversized(rest, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_G1_CURVE)?;

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

    pub fn g1_mul<'a>(input: &'a [u8]) -> Result<[u8; SERIALIZED_G1_POINT_BYTE_LENGTH], ApiError> {
        if input.len() != SERIALIZED_G1_POINT_BYTE_LENGTH + SCALAR_BYTE_LENGTH {
            return Err(ApiError::InputError("invalid input length for G1 multiplication".to_owned()));
        }

        let (p_0, rest) = decode_g1::decode_g1_point_from_xy_oversized(input, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_G1_CURVE)?;
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

    pub fn g1_multiexp<'a>(input: &'a [u8]) -> Result<[u8; SERIALIZED_G1_POINT_BYTE_LENGTH], ApiError> {
        if input.len() % (SERIALIZED_G1_POINT_BYTE_LENGTH + SCALAR_BYTE_LENGTH) != 0 {
            return Err(ApiError::InputError("invalid input length for G1 multiplication".to_owned()));
        }
        let num_pairs = input.len() / (SERIALIZED_G1_POINT_BYTE_LENGTH + SCALAR_BYTE_LENGTH);

        if num_pairs == 0 {
            return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
        }

        let mut global_rest = input;
        let mut bases = Vec::with_capacity(num_pairs);
        let mut scalars = Vec::with_capacity(num_pairs);

        for _ in 0..num_pairs {
            let (p, local_rest) = decode_g1::decode_g1_point_from_xy_oversized(global_rest, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_G1_CURVE)?;
            let (scalar, local_rest) = decode_g1::decode_scalar_representation(local_rest, SCALAR_BYTE_LENGTH)?;
            if !p.is_on_curve() {
                if !crate::features::in_fuzzing_or_gas_metering() {
                    return Err(ApiError::InputError(format!("Point is not on curve, file {}, line {}", file!(), line!())));
                }
            }
            bases.push(p);
            scalars.push(scalar);
            global_rest = local_rest;
        }

        if bases.len() != scalars.len() || bases.len() == 0 {
            return Err(ApiError::InputError(format!("Multiexp with empty input pairs, file {}, line {}", file!(), line!())));
        } 

        let result = peppinger(&bases, scalars);

        let mut output = [0u8; SERIALIZED_G1_POINT_BYTE_LENGTH];

        let as_vec = decode_g1::serialize_g1_point(SERIALIZED_FP_BYTE_LENGTH, &result)?;

        output.copy_from_slice(&as_vec[..]);

        Ok(output)
    }

    pub fn g2_add<'a>(input: &'a [u8]) -> Result<[u8; SERIALIZED_G2_POINT_BYTE_LENGTH], ApiError> {
        if input.len() != SERIALIZED_G2_POINT_BYTE_LENGTH * 2 {
            return Err(ApiError::InputError("invalid input length for G2 addition".to_owned()));
        }

        let (mut p_0, rest) = decode_g2::decode_g2_point_from_xy_in_fp2_oversized(input, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_G2_CURVE)?;
        let (p_1, _) = decode_g2::decode_g2_point_from_xy_in_fp2_oversized(rest, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_G2_CURVE)?;

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

        let mut output = [0u8; SERIALIZED_G2_POINT_BYTE_LENGTH];

        let as_vec = decode_g2::serialize_g2_point_in_fp2(SERIALIZED_FP_BYTE_LENGTH, &p_0)?;

        output.copy_from_slice(&as_vec[..]);

        Ok(output)
    }

    pub fn g2_mul<'a>(input: &'a [u8]) -> Result<[u8; SERIALIZED_G2_POINT_BYTE_LENGTH], ApiError> {
        if input.len() != SERIALIZED_G2_POINT_BYTE_LENGTH + SCALAR_BYTE_LENGTH {
            return Err(ApiError::InputError("invalid input length for G1 multiplication".to_owned()));
        }

        let (p_0, rest) = decode_g2::decode_g2_point_from_xy_in_fp2_oversized(input, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_G2_CURVE)?;
        let (scalar, _) = decode_g1::decode_scalar_representation(rest, SCALAR_BYTE_LENGTH)?;

        if !p_0.is_on_curve() {
            if !crate::features::in_fuzzing_or_gas_metering() {
                return Err(ApiError::InputError(format!("Point is not on curve, file {}, line {}", file!(), line!())));
            }
        }

        let p = p_0.mul(&scalar);

        let mut output = [0u8; SERIALIZED_G2_POINT_BYTE_LENGTH];

        let as_vec = decode_g2::serialize_g2_point_in_fp2(SERIALIZED_FP_BYTE_LENGTH, &p)?;

        output.copy_from_slice(&as_vec[..]);

        Ok(output)
    }

    pub fn g2_multiexp<'a>(input: &'a [u8]) -> Result<[u8; SERIALIZED_G2_POINT_BYTE_LENGTH], ApiError> {
        if input.len() % (SERIALIZED_G2_POINT_BYTE_LENGTH + SCALAR_BYTE_LENGTH) != 0 {
            return Err(ApiError::InputError("invalid input length for G1 multiplication".to_owned()));
        }
        let num_pairs = input.len() / (SERIALIZED_G2_POINT_BYTE_LENGTH + SCALAR_BYTE_LENGTH);

        if num_pairs == 0 {
            return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
        }

        let mut global_rest = input;
        let mut bases = Vec::with_capacity(num_pairs);
        let mut scalars = Vec::with_capacity(num_pairs);

        for _ in 0..num_pairs {
            let (p, local_rest) = decode_g2::decode_g2_point_from_xy_in_fp2_oversized(global_rest, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_G2_CURVE)?;
            let (scalar, local_rest) = decode_g1::decode_scalar_representation(local_rest, SCALAR_BYTE_LENGTH)?;
            if !p.is_on_curve() {
                if !crate::features::in_fuzzing_or_gas_metering() {
                    return Err(ApiError::InputError(format!("Point is not on curve, file {}, line {}", file!(), line!())));
                }
            }
            bases.push(p);
            scalars.push(scalar);
            global_rest = local_rest;
        }

        if bases.len() != scalars.len() || bases.len() == 0 {
            return Err(ApiError::InputError(format!("Multiexp with empty input pairs, file {}, line {}", file!(), line!())));
        } 

        let result = peppinger(&bases, scalars);

        let mut output = [0u8; SERIALIZED_G2_POINT_BYTE_LENGTH];

        let as_vec = decode_g2::serialize_g2_point_in_fp2(SERIALIZED_FP_BYTE_LENGTH, &result)?;

        output.copy_from_slice(&as_vec[..]);

        Ok(output)
    }

    pub fn pair<'a>(input: &'a [u8]) -> Result<[u8; 32], ApiError> {
        if input.len() % (SERIALIZED_G2_POINT_BYTE_LENGTH + SERIALIZED_G1_POINT_BYTE_LENGTH) != 0 {
            return Err(ApiError::InputError("invalid input length for pairing".to_owned()));
        }
        let num_pairs = input.len() / (SERIALIZED_G2_POINT_BYTE_LENGTH + SERIALIZED_G1_POINT_BYTE_LENGTH);

        if num_pairs == 0 {
            return Err(ApiError::InputError("Invalid number of pairs".to_owned()));
        }

        let mut global_rest = input;

        let mut g1_points = Vec::with_capacity(num_pairs);
        let mut g2_points = Vec::with_capacity(num_pairs);

        for _ in 0..num_pairs {
            let (g1, rest) = decode_g1::decode_g1_point_from_xy_oversized(global_rest, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_G1_CURVE)?;
            let (g2, rest) = decode_g2::decode_g2_point_from_xy_in_fp2_oversized(rest, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_G2_CURVE)?;

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

            if !g1.check_correct_subgroup() {
                if !crate::features::in_fuzzing_or_gas_metering() {
                    return Err(ApiError::InputError("G1 or G2 point is not in the expected subgroup".to_owned()));
                }
            }

            if !g2.check_correct_subgroup() {
                if !crate::features::in_fuzzing_or_gas_metering() {
                    return Err(ApiError::InputError("G1 or G2 point is not in the expected subgroup".to_owned()));
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

        let engine = &bls12_381::BLS12_381_PAIRING_ENGINE;

        let pairing_result = engine.pair(&g1_points, &g2_points);

        if pairing_result.is_none() {
            return Err(ApiError::UnknownParameter("Pairing engine returned no value".to_owned()));
        }

        use crate::extension_towers::fp12_as_2_over3_over_2::Fp12;
        use crate::traits::ZeroAndOne;

        let one_fp12 = Fp12::one(&bls12_381::BLS12_381_EXTENSION_12_FIELD);
        let pairing_result = pairing_result.unwrap();
        let result = if pairing_result == one_fp12 {
            pairing_result_true()
        } else {
            pairing_result_false()
        };

        Ok(result)
    }

    pub fn map_to_curve<'a>(input: &'a [u8]) -> Result<Vec<u8>, ApiError> {
        if input.len() == SERIALIZED_FP_BYTE_LENGTH {
            return Self::map_fp_to_g1(input);
        } else if input.len() == SERIALIZED_FP2_BYTE_LENGTH {
            return Self::map_fp2_to_g2(input);
        } else {
            return Err(ApiError::InputError("invalid input length for field to curve mapping".to_owned()));
        }
    }

    fn map_fp_to_g1<'a>(input: &'a [u8]) -> Result<Vec<u8>, ApiError> {
        let (fe, _) = decode_fp::decode_fp_oversized(input, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_FIELD)?;
        let point = mapping::fp_to_g1(&fe)?;

        let as_vec = decode_g1::serialize_g1_point(SERIALIZED_FP_BYTE_LENGTH, &point)?;

        Ok(as_vec)
    }

    fn map_fp2_to_g2<'a>(input: &'a [u8]) -> Result<Vec<u8>, ApiError> {
        let (fe, _) = decode_fp::decode_fp2_oversized(input, SERIALIZED_FP_BYTE_LENGTH, &bls12_381::BLS12_381_EXTENSION_2_FIELD)?;
        let point = mapping::fp2_to_g2(&fe)?;

        let as_vec = decode_g2::serialize_g2_point_in_fp2(SERIALIZED_FP_BYTE_LENGTH, &point)?;

        Ok(as_vec)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::{Rng};
    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    use indicatif::{ProgressBar, ProgressStyle};

    use csv::Writer;
    use hex;

    use num_bigint::BigUint;
    use num_traits::Num;
    use crate::fp::Fp;

    type Scalar = crate::integers::MaxGroupSizeUint;

    type FpElement = crate::fp::Fp<'static, crate::field::U384Repr, crate::field::PrimeField<crate::field::U384Repr>>;
    type Fp2Element = crate::extension_towers::fp2::Fp2<'static, crate::field::U384Repr, crate::field::PrimeField<crate::field::U384Repr>>;

    type G1 = crate::weierstrass::curve::CurvePoint<'static, crate::weierstrass::CurveOverFpParameters<'static, crate::field::U384Repr, crate::field::PrimeField<crate::field::U384Repr>>>;
    type G2 = crate::weierstrass::curve::CurvePoint<'static, crate::weierstrass::CurveOverFp2Parameters<'static, crate::field::U384Repr, crate::field::PrimeField<crate::field::U384Repr>>>;

    fn make_random_fp_with_encoding<R: Rng>(rng: &mut R, modulus: &BigUint) -> (FpElement, Vec<u8>) {
        let mut buff = vec![0u8; 48*3];
        rng.fill_bytes(&mut buff);

        let num = BigUint::from_bytes_be(&buff);
        let num = num % modulus.clone();

        let x = Fp::from_be_bytes(&bls12_381::BLS12_381_FIELD, &num.to_bytes_be(), true).unwrap();

        let as_vec = decode_fp::serialize_fp_fixed_len(SERIALIZED_FP_BYTE_LENGTH, &x).unwrap();

        assert!(as_vec.len() == SERIALIZED_FP_BYTE_LENGTH);
        assert_eq!(&as_vec[..16], &[0u8; 16]);

        (x, as_vec)
    }

    fn make_random_fp2_with_encoding<R: Rng>(rng: &mut R, modulus: &BigUint) -> (Fp2Element, Vec<u8>) {
        let mut buff = vec![0u8; 48*3];

        let mut encoding = Vec::with_capacity(SERIALIZED_FP2_BYTE_LENGTH);

        rng.fill_bytes(&mut buff);
        let num = BigUint::from_bytes_be(&buff);
        let num = num % modulus.clone();

        let c0 = Fp::from_be_bytes(&bls12_381::BLS12_381_FIELD, &num.to_bytes_be(), true).unwrap();

        let as_vec = decode_fp::serialize_fp_fixed_len(SERIALIZED_FP_BYTE_LENGTH, &c0).unwrap();
        encoding.extend(as_vec);

        rng.fill_bytes(&mut buff);
        let num = BigUint::from_bytes_be(&buff);
        let num = num % modulus.clone();

        let c1 = Fp::from_be_bytes(&bls12_381::BLS12_381_FIELD, &num.to_bytes_be(), true).unwrap();

        let as_vec = decode_fp::serialize_fp_fixed_len(SERIALIZED_FP_BYTE_LENGTH, &c1).unwrap();
        encoding.extend(as_vec);

        assert!(encoding.len() == SERIALIZED_FP2_BYTE_LENGTH);
        assert_eq!(&encoding[..16], &[0u8; 16]);
        assert_eq!(&encoding[64..80], &[0u8; 16]);

        let mut fe = bls12_381::BLS12_381_FP2_ZERO.clone();
        fe.c0 = c0;
        fe.c1 = c1;

        (fe, encoding)
    }
    
    fn make_random_g1_with_encoding<R: Rng>(rng: &mut R) -> (G1, Vec<u8>) {
        let (scalar, _) = make_random_scalar_with_encoding(rng);
        let p = bls12_381::BLS12_381_G1_GENERATOR.mul(&scalar);

        let as_vec = decode_g1::serialize_g1_point(SERIALIZED_FP_BYTE_LENGTH, &p).unwrap();

        assert!(as_vec.len() == SERIALIZED_G1_POINT_BYTE_LENGTH);
        assert_eq!(&as_vec[..16], &[0u8; 16]);
        assert_eq!(&as_vec[64..80], &[0u8; 16]);

        (p, as_vec)
    }

    fn make_random_g2_with_encoding<R: Rng>(rng: &mut R) -> (G2, Vec<u8>) {
        let (scalar, _) = make_random_scalar_with_encoding(rng);

        let p = bls12_381::BLS12_381_G2_GENERATOR.mul(&scalar);

        let as_vec = decode_g2::serialize_g2_point_in_fp2(SERIALIZED_FP_BYTE_LENGTH, &p).unwrap();

        assert!(as_vec.len() == SERIALIZED_G2_POINT_BYTE_LENGTH);

        (p, as_vec)
    }

    fn make_random_scalar_with_encoding<R: Rng>(rng: &mut R) -> (Scalar, Vec<u8>) {
        let mut buff = vec![0u8; SCALAR_BYTE_LENGTH];
        rng.fill_bytes(&mut buff);

        let (scalar, _) = decode_g1::decode_scalar_representation(&buff, SCALAR_BYTE_LENGTH).unwrap();

        (scalar, buff)
    }

    fn make_csv_writer(path: &str) -> Option<Writer<std::fs::File>> {
        if WRITE_VECTORS {
            let mut writer = Writer::from_path(path).expect("must open a test file");
            writer.write_record(&["input", "result"]).expect("must write header");

            Some(writer)
        } else {
            None
        }
    }

    const NUM_TESTS: usize = 100;
    const MULTIEXP_INPUT: usize = 16;
    const WRITE_VECTORS: bool = true;

    #[test]
    fn test_g1_add() {
        let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        pb.set_length(NUM_TESTS as u64);

        let mut writer = make_csv_writer("src/test/test_vectors/eip2537/g1_add.csv");

        for _ in 0..NUM_TESTS {
            let mut encoding = Vec::with_capacity(SERIALIZED_G1_POINT_BYTE_LENGTH * 2);

            let (mut p0, e) = make_random_g1_with_encoding(&mut rng);
            encoding.extend(e);

            let (p1, e) = make_random_g1_with_encoding(&mut rng);
            encoding.extend(e);

            p0.add_assign(&p1);

            let expected = decode_g1::serialize_g1_point(SERIALIZED_FP_BYTE_LENGTH, &p0).unwrap();
            assert!(expected.len() == SERIALIZED_G1_POINT_BYTE_LENGTH);

            let api_result = EIP2537Executor::g1_add(&encoding).unwrap();

            assert_eq!(&expected[..], &api_result[..]);

            if let Some(writer) = writer.as_mut() {
                writer.write_record(
                    &[
                        &hex::encode(&encoding[..]), 
                        &hex::encode(&api_result[..])
                    ],
                ).expect("must write a test vector");
            }

            pb.inc(1);
        }

        pb.finish_with_message("Completed");
    }

    #[test]
    fn test_g1_point_mul() {
        let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        pb.set_length(NUM_TESTS as u64);

        let mut writer = make_csv_writer("src/test/test_vectors/eip2537/g1_mul.csv");

        for _ in 0..NUM_TESTS {
            let mut encoding = Vec::with_capacity(SERIALIZED_G1_POINT_BYTE_LENGTH + SCALAR_BYTE_LENGTH);

            let (p0, e) = make_random_g1_with_encoding(&mut rng);
            encoding.extend(e);

            let (scalar, e) = make_random_scalar_with_encoding(&mut rng);
            encoding.extend(e);

            let p = p0.mul(&scalar);

            let expected = decode_g1::serialize_g1_point(SERIALIZED_FP_BYTE_LENGTH, &p).unwrap();
            assert!(expected.len() == SERIALIZED_G1_POINT_BYTE_LENGTH);

            let api_result = EIP2537Executor::g1_mul(&encoding).unwrap();

            assert_eq!(&expected[..], &api_result[..]);

            if let Some(writer) = writer.as_mut() {
                writer.write_record(
                    &[
                        &hex::encode(&encoding[..]), 
                        &hex::encode(&api_result[..])
                    ],
                ).expect("must write a test vector");
            }

            pb.inc(1);
        }

        pb.finish_with_message("Completed");
    }

    #[test]
    fn test_g1_multiexp() {
        let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        pb.set_length(NUM_TESTS as u64);

        let mut writer = make_csv_writer("src/test/test_vectors/eip2537/g1_multiexp.csv");

        for _ in 0..NUM_TESTS {
            let mut encoding = Vec::with_capacity((SERIALIZED_G1_POINT_BYTE_LENGTH + SCALAR_BYTE_LENGTH) * MULTIEXP_INPUT);
            let mut points = Vec::with_capacity(MULTIEXP_INPUT);
            let mut scalars = Vec::with_capacity(MULTIEXP_INPUT);
            for _ in 0..MULTIEXP_INPUT {
                let (p0, e) = make_random_g1_with_encoding(&mut rng);
                encoding.extend(e);

                let (scalar, e) = make_random_scalar_with_encoding(&mut rng);
                encoding.extend(e);

                points.push(p0);
                scalars.push(scalar);
            }

            let p = peppinger(&points, scalars);

            let expected = decode_g1::serialize_g1_point(SERIALIZED_FP_BYTE_LENGTH, &p).unwrap();
            assert!(expected.len() == SERIALIZED_G1_POINT_BYTE_LENGTH);

            let api_result = EIP2537Executor::g1_multiexp(&encoding).unwrap();

            assert_eq!(&expected[..], &api_result[..]);

            if let Some(writer) = writer.as_mut() {
                writer.write_record(
                    &[
                        &hex::encode(&encoding[..]), 
                        &hex::encode(&api_result[..])
                    ],
                ).expect("must write a test vector");
            }

            pb.inc(1);
        }

        pb.finish_with_message("Completed");
    }

    #[test]
    fn test_g2_add() {
        let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        pb.set_length(NUM_TESTS as u64);

        let mut writer = make_csv_writer("src/test/test_vectors/eip2537/g2_add.csv");

        for _ in 0..NUM_TESTS {
            let mut encoding = Vec::with_capacity(SERIALIZED_G2_POINT_BYTE_LENGTH * 2);

            let (mut p0, e) = make_random_g2_with_encoding(&mut rng);
            encoding.extend(e);

            let (p1, e) = make_random_g2_with_encoding(&mut rng);
            encoding.extend(e);

            p0.add_assign(&p1);

            let expected = decode_g2::serialize_g2_point_in_fp2(SERIALIZED_FP_BYTE_LENGTH, &p0).unwrap();
            assert!(expected.len() == SERIALIZED_G2_POINT_BYTE_LENGTH);

            let api_result = EIP2537Executor::g2_add(&encoding).unwrap();

            assert_eq!(&expected[..], &api_result[..]);

            if let Some(writer) = writer.as_mut() {
                writer.write_record(
                    &[
                        &hex::encode(&encoding[..]), 
                        &hex::encode(&api_result[..])
                    ],
                ).expect("must write a test vector");
            }

            pb.inc(1);
        }

        pb.finish_with_message("Completed");
    }

    #[test]
    fn test_g2_point_mul() {
        let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        pb.set_length(NUM_TESTS as u64);

        let mut writer = make_csv_writer("src/test/test_vectors/eip2537/g2_mul.csv");

        for _ in 0..NUM_TESTS {
            let mut encoding = Vec::with_capacity(SERIALIZED_G2_POINT_BYTE_LENGTH + SCALAR_BYTE_LENGTH);

            let (p0, e) = make_random_g2_with_encoding(&mut rng);
            encoding.extend(e);

            let (scalar, e) = make_random_scalar_with_encoding(&mut rng);
            encoding.extend(e);

            let p = p0.mul(&scalar);

            let expected = decode_g2::serialize_g2_point_in_fp2(SERIALIZED_FP_BYTE_LENGTH, &p).unwrap();
            assert!(expected.len() == SERIALIZED_G2_POINT_BYTE_LENGTH);

            let api_result = EIP2537Executor::g2_mul(&encoding).unwrap();

            assert_eq!(&expected[..], &api_result[..]);

            if let Some(writer) = writer.as_mut() {
                writer.write_record(
                    &[
                        &hex::encode(&encoding[..]), 
                        &hex::encode(&api_result[..])
                    ],
                ).expect("must write a test vector");
            }

            pb.inc(1);
        }

        pb.finish_with_message("Completed");
    }

    #[test]
    fn test_g2_multiexp() {
        let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        pb.set_length(NUM_TESTS as u64);

        let mut writer = make_csv_writer("src/test/test_vectors/eip2537/g2_multiexp.csv");

        for _ in 0..NUM_TESTS {
            let mut encoding = Vec::with_capacity((SERIALIZED_G2_POINT_BYTE_LENGTH + SCALAR_BYTE_LENGTH) * MULTIEXP_INPUT);
            let mut points = Vec::with_capacity(MULTIEXP_INPUT);
            let mut scalars = Vec::with_capacity(MULTIEXP_INPUT);
            for _ in 0..MULTIEXP_INPUT {
                let (p0, e) = make_random_g2_with_encoding(&mut rng);
                encoding.extend(e);

                let (scalar, e) = make_random_scalar_with_encoding(&mut rng);
                encoding.extend(e);

                points.push(p0);
                scalars.push(scalar);
            }

            let p = peppinger(&points, scalars);

            let expected = decode_g2::serialize_g2_point_in_fp2(SERIALIZED_FP_BYTE_LENGTH, &p).unwrap();
            assert!(expected.len() == SERIALIZED_G2_POINT_BYTE_LENGTH);

            let api_result = EIP2537Executor::g2_multiexp(&encoding).unwrap();

            assert_eq!(&expected[..], &api_result[..]);

            if let Some(writer) = writer.as_mut() {
                writer.write_record(
                    &[
                        &hex::encode(&encoding[..]), 
                        &hex::encode(&api_result[..])
                    ],
                ).expect("must write a test vector");
            }

            pb.inc(1);
        }

        pb.finish_with_message("Completed");
    }

    #[test]
    fn generate_fp_to_g1_mapping_vectors() {
        let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        pb.set_length(NUM_TESTS as u64);

        let mut writer = make_csv_writer("src/test/test_vectors/eip2537/fp_to_g1.csv");
        assert!(writer.is_some());
        let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();

        for _ in 0..NUM_TESTS {
            let (_, input) = make_random_fp_with_encoding(&mut rng, &modulus);

            let api_result = EIP2537Executor::map_to_curve(&input).unwrap();
            assert!(api_result.len() == SERIALIZED_G1_POINT_BYTE_LENGTH);

            if let Some(writer) = writer.as_mut() {
                writer.write_record(
                    &[
                        &hex::encode(&input[..]), 
                        &hex::encode(&api_result[..])
                    ],
                ).expect("must write a test vector");
            }

            pb.inc(1);
        }

        pb.finish_with_message("Completed");
    }

    #[test]
    fn generate_fp2_to_g2_mapping_vectors() {
        let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        pb.set_length(NUM_TESTS as u64);

        let mut writer = make_csv_writer("src/test/test_vectors/eip2537/fp2_to_g2.csv");
        assert!(writer.is_some());
        let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();

        for _ in 0..NUM_TESTS {
            let (_, input) = make_random_fp2_with_encoding(&mut rng, &modulus);

            let api_result = EIP2537Executor::map_to_curve(&input).unwrap();
            assert!(api_result.len() == SERIALIZED_G2_POINT_BYTE_LENGTH);

            if let Some(writer) = writer.as_mut() {
                writer.write_record(
                    &[
                        &hex::encode(&input[..]), 
                        &hex::encode(&api_result[..])
                    ],
                ).expect("must write a test vector");
            }

            pb.inc(1);
        }

        pb.finish_with_message("Completed");
    }
}