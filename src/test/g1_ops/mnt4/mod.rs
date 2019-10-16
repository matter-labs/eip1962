use crate::public_interface::constants::*;

use crate::test::parsers::*;

use super::*; 

pub(crate) fn assemble_single_curve_params(curve: JsonMnt4PairingCurveParameters) -> (Vec<u8>, usize, usize) {
    // - Lengths of modulus (in bytes)
    // - Field modulus
    // - Curve A
    // - Curve B

    // first determine the length of the modulus
    let modulus = curve.q;
    let modulus_length = modulus.clone().to_bytes_be().len();

    let modulus_len_encoded = vec![modulus_length as u8];
    let modulus_encoded = pad_for_len_be(modulus.clone().to_bytes_be(), modulus_length);

    let a_encoded = apply_sign(curve.a, &modulus);
    let a_encoded = pad_for_len_be(a_encoded.to_bytes_be(), modulus_length);

    let b_encoded = apply_sign(curve.b, &modulus);
    let b_encoded = pad_for_len_be(b_encoded.to_bytes_be(), modulus_length);

    let group_size = curve.r;
    let group_size_encoded = group_size.clone().to_bytes_be();
    let group_size_length = group_size_encoded.len();
    let group_len_encoded = vec![group_size_length as u8];

    let mut calldata = vec![];
    calldata.extend(modulus_len_encoded.into_iter());
    calldata.extend(modulus_encoded.into_iter());
    calldata.extend(a_encoded.into_iter());
    calldata.extend(b_encoded.into_iter());
    calldata.extend(group_len_encoded.into_iter());
    calldata.extend(group_size_encoded.into_iter());

    (calldata, modulus_length, group_size_length)
}

pub(crate) fn assemble_single_point_scalar_pair(
    pair: JsonG1PointScalarMultiplicationPair,
    modulus_len: usize,
    group_len: usize,
) -> (Vec<u8>, Vec<u8>) {
    // - X,
    // - Y,
    // - Scalar

    // - Result X
    // - Result Y

    // first determine the length of the modulus
    let x_encoded = pad_for_len_be(pair.base_x.to_bytes_be(), modulus_len);
    let y_encoded = pad_for_len_be(pair.base_y.to_bytes_be(), modulus_len);

    let scalar_encoded = pad_for_len_be(pair.scalar.to_bytes_be(), group_len);

    let mut calldata = vec![];
    calldata.extend(x_encoded.into_iter());
    calldata.extend(y_encoded.into_iter());
    calldata.extend(scalar_encoded.into_iter());

    let result_x_encoded = pad_for_len_be(pair.result_x.to_bytes_be(), modulus_len);
    let result_y_encoded = pad_for_len_be(pair.result_y.to_bytes_be(), modulus_len);

    let mut result = vec![];
    result.extend(result_x_encoded.into_iter());
    result.extend(result_y_encoded.into_iter());

    (calldata, result)
}

pub(crate) fn assemble_single_points_addition_pair(
    pair: JsonG1PointScalarMultiplicationPair,
    modulus_len: usize,
    _group_len: usize,
) -> (Vec<u8>, Vec<u8>) {
    // - X1,
    // - Y1,
    // - X2,
    // - Y2,

    // - Result X
    // - Result Y
    let mut calldata = vec![];
    let x_encoded = pad_for_len_be(pair.base_x.to_bytes_be(), modulus_len);
    let y_encoded = pad_for_len_be(pair.base_y.to_bytes_be(), modulus_len);
    calldata.extend(x_encoded.into_iter());
    calldata.extend(y_encoded.into_iter());

    let result_x_encoded = pad_for_len_be(pair.result_x.to_bytes_be(), modulus_len);
    let result_y_encoded = pad_for_len_be(pair.result_y.to_bytes_be(), modulus_len);

    let mut result = vec![];
    result.extend(result_x_encoded.into_iter());
    result.extend(result_y_encoded.into_iter());

    (calldata, result)
}

// #[test]
// fn test_g1_mul_from_vectors() {
//     let curves = read_dir_and_grab_curves::<JsonBnPairingCurveParameters>("src/test/test_vectors/bn/");
//     assert!(curves.len() != 0);
//     for (curve, _) in curves.into_iter() {
//         let (calldata, modulus_len, group_len) = assemble_single_curve_params(curve.clone());
//         for pair in curve.g1_mul_vectors.into_iter() {
//             let (points_data, expected_result) = assemble_single_point_scalar_pair(pair, modulus_len, group_len);

//             let mut calldata = calldata.clone();
//             calldata.extend(points_data);

//             let result = call_g1_engine_mul(&calldata[..]);
//             if result.is_err() {
//                 panic!("{}", result.err().unwrap());
//             }
//             assert!(result.is_ok());

//             let result = result.unwrap();
//             assert!(result == expected_result);
//         }
//     }
// }

// extern crate hex;
// extern crate csv;

// use hex::{encode};
// use csv::{Writer};

// #[test]
// fn dump_g1_mul_vectors() {
//     let curves = read_dir_and_grab_curves::<JsonBnPairingCurveParameters>("src/test/test_vectors/bn/");
//     assert!(curves.len() != 0);
//     let mut writer = Writer::from_path("src/test/test_vectors/bn/g1_mul.csv").expect("must open a test file");
//     writer.write_record(&["input", "result"]).expect("must write header");
//     for (curve, _) in curves.into_iter() {
//         let (calldata, modulus_len, group_len) = assemble_single_curve_params(curve.clone());
//         for pair in curve.g1_mul_vectors.into_iter() {
//             let (points_data, expected_result) = assemble_single_point_scalar_pair(pair, modulus_len, group_len);
//             let mut input_data = vec![OPERATION_G1_MUL];
//             input_data.extend(calldata.clone());
//             input_data.extend(points_data);

//             writer.write_record(&[
//                 prepend_0x(&encode(&input_data[..])), 
//                 prepend_0x(&encode(&expected_result[..]))],
//             ).expect("must write a record");
//         }
//     }
//     writer.flush().expect("must finalize writing");
// }

// use rust_test::Bencher;

// #[bench]
// fn bench_single(b: &mut Bencher) {
//     let calldata = assemble_single();
//     b.iter(|| {
//         call_bls12_engine(&calldata[..]).expect("must use");
//     });
// }

