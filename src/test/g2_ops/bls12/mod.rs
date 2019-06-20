extern crate hex;

use crate::public_interface::constants::*;

use num_bigint::BigUint;
use crate::test::parsers::*;

use super::*;

const EXTENSION_DEGREE: usize = 2;

fn assemble_single_curve_params(curve: JsonCurveParameters) -> Vec<u8> {
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

    let mut a_encoded = pad_for_len_be(curve.a_0.to_bytes_be(), modulus_length);
    a_encoded.extend(pad_for_len_be(curve.a_1.to_bytes_be(), modulus_length));
    let mut b_encoded = pad_for_len_be(curve.b_0.to_bytes_be(), modulus_length);
    b_encoded.extend(pad_for_len_be(curve.b_1.to_bytes_be(), modulus_length));

    let twist_type = if curve.is_d_type { vec![TWIST_TYPE_D] } else { vec![TWIST_TYPE_M] };


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

    calldata
}

#[test]
fn test_from_vectors() {
    let curves = read_dir_and_grab_curves("src/test/test_vectors/bls12/");
    assert!(curves.len() != 0);
    for curve in curves.into_iter() {
        let calldata = assemble_single_curve_params(curve);
        let result = call_bls12_engine(&calldata[..]);
        assert!(result.is_ok());

        let result = result.unwrap()[0];
        assert!(result == 1u8);
    }
}

use rust_test::Bencher;

#[bench]
fn bench_single(b: &mut Bencher) {
    let calldata = assemble_single();
    b.iter(|| {
        call_bls12_engine(&calldata[..]).expect("must use");
    });
}

