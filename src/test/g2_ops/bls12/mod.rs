use crate::test::parsers::*;
use crate::public_interface::constants::*;

use super::*;

const EXTENSION_DEGREE: usize = 2;

pub(crate) fn assemble_single_curve_params(curve: JsonBls12PairingCurveParameters) -> (Vec<u8>, usize, usize) {
    // - Lengths of modulus (in bytes)
    // - Field modulus
    // - Extension degree
    // - Non-redisue
    // - Curve A
    // - Curve B

    // first determine the length of the modulus
    let modulus = curve.q;
    let modulus_length = modulus.clone().to_bytes_be().len();

    let modulus_len_encoded = vec![modulus_length as u8];
    let modulus_encoded = pad_for_len_be(modulus.clone().to_bytes_be(), modulus_length);

    let encoded_extension_degree = vec![EXTENSION_DEGREE as u8];

    let fp2_nonres_encoded = {
        let (mut nonres, is_positive) = curve.non_residue;
        if !is_positive {
            nonres = modulus.clone() - nonres;
        }
        pad_for_len_be(nonres.to_bytes_be(), modulus_length)
    };

    let mut a_encoded = pad_for_len_be(curve.a_twist_0.to_bytes_be(), modulus_length);
    a_encoded.extend(pad_for_len_be(curve.a_twist_1.to_bytes_be(), modulus_length));
    let mut b_encoded = pad_for_len_be(curve.b_twist_0.to_bytes_be(), modulus_length);
    b_encoded.extend(pad_for_len_be(curve.b_twist_1.to_bytes_be(), modulus_length));

    // now we make two random scalars and do scalar multiplications in G1 and G2 to get pairs that should
    // at the end of the day pair to identity element

    let group_size = curve.r;
    let group_size_encoded = group_size.clone().to_bytes_be();
    let group_size_length = group_size_encoded.len();
    let group_len_encoded = vec![group_size_length as u8];

    let mut calldata = vec![];
    calldata.extend(modulus_len_encoded.into_iter());
    calldata.extend(modulus_encoded.into_iter());
    calldata.extend(encoded_extension_degree.into_iter());
    calldata.extend(fp2_nonres_encoded.into_iter());
    calldata.extend(a_encoded.into_iter());
    calldata.extend(b_encoded.into_iter());
    calldata.extend(group_len_encoded.into_iter());
    calldata.extend(group_size_encoded.into_iter());

    (calldata, modulus_length, group_size_length)
}

fn assemble_single_point_scalar_pair(
    pair: JsonG2PointScalarMultiplicationPair,
    modulus_len: usize,
    group_len: usize,
) -> (Vec<u8>, Vec<u8>) {
    // - X,
    // - Y,
    // - Scalar

    // - Result X
    // - Result Y

    // first determine the length of the modulus
    let mut x_encoded = pad_for_len_be(pair.base_x_0.to_bytes_be(), modulus_len);
    x_encoded.extend(pad_for_len_be(pair.base_x_1.to_bytes_be(), modulus_len));
    let mut y_encoded = pad_for_len_be(pair.base_y_0.to_bytes_be(), modulus_len);
    y_encoded.extend(pad_for_len_be(pair.base_y_1.to_bytes_be(), modulus_len));

    let scalar_encoded = pad_for_len_be(pair.scalar.to_bytes_be(), group_len);

    let mut calldata = vec![];
    calldata.extend(x_encoded.into_iter());
    calldata.extend(y_encoded.into_iter());
    calldata.extend(scalar_encoded.into_iter());

    let mut result_x_encoded = pad_for_len_be(pair.result_x_0.to_bytes_be(), modulus_len);
    result_x_encoded.extend(pad_for_len_be(pair.result_x_1.to_bytes_be(), modulus_len));
    let mut result_y_encoded = pad_for_len_be(pair.result_y_0.to_bytes_be(), modulus_len);
    result_y_encoded.extend(pad_for_len_be(pair.result_y_1.to_bytes_be(), modulus_len));

    let mut result = vec![];
    result.extend(result_x_encoded.into_iter());
    result.extend(result_y_encoded.into_iter());

    (calldata, result)
}

#[test]
fn test_g2_mul_from_vectors() {
    let curves = read_dir_and_grab_curves::<JsonBls12PairingCurveParameters>("src/test/test_vectors/bls12/");
    assert!(curves.len() != 0);
    for (curve, _) in curves.into_iter() {
        let (calldata, modulus_len, group_len) = assemble_single_curve_params(curve.clone());
        for pair in curve.g2_mul_vectors.into_iter() {
            let (points_data, expected_result) = assemble_single_point_scalar_pair(pair, modulus_len, group_len);

            let mut calldata = calldata.clone();
            calldata.extend(points_data);

            let result = call_g2_engine_mul(&calldata[..]);
            if result.is_err() {
                panic!("{}", result.err().expect("guaranteed to exist"));
            }
            assert!(result.is_ok());

            let result = result.expect("guaranteed to exist");
            assert!(result == expected_result);
        }
    }
}

extern crate hex;
extern crate csv;

use hex::{encode};
use csv::{Writer};

#[test]
#[ignore]
fn dump_g2_mul_vectors() {
    let curves = read_dir_and_grab_curves::<JsonBls12PairingCurveParameters>("src/test/test_vectors/bls12/");
    assert!(curves.len() != 0);
    let mut writer = Writer::from_path("src/test/test_vectors/bls12/g2_mul.csv").expect("must open a test file");
    writer.write_record(&["input", "result"]).expect("must write header");
    for (curve, _) in curves.into_iter() {
        let (calldata, modulus_len, group_len) = assemble_single_curve_params(curve.clone());
        for pair in curve.g2_mul_vectors.into_iter() {
            let (points_data, expected_result) = assemble_single_point_scalar_pair(pair, modulus_len, group_len);
            let mut input_data = vec![OPERATION_G2_MUL];
            input_data.extend(calldata.clone());
            input_data.extend(points_data);

            writer.write_record(&[
                prepend_0x(&encode(&input_data[..])), 
                prepend_0x(&encode(&expected_result[..]))],
            ).expect("must write a record");
        }
    }
    writer.flush().expect("must finalize writing");
}


// use rust_test::Bencher;

// #[bench]
// fn bench_single(b: &mut Bencher) {
//     let calldata = assemble_single();
//     b.iter(|| {
//         call_bls12_engine(&calldata[..]).expect("must use");
//     });
// }

