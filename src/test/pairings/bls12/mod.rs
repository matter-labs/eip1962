use crate::public_interface::constants::*;
use crate::public_interface::{PublicG1Api, G1Api};

use num_bigint::BigUint;

use crate::test::parsers::*;
use super::call_pairing_engine;

pub(crate) fn assemble_single_curve_params(curve: JsonBls12PairingCurveParameters) -> Vec<u8> {
    // - Curve type
    // - Lengths of modulus (in bytes)
    // - Field modulus
    // - Curve A
    // - Curve B
    // - non-residue for Fp2
    // - non-residue for Fp6
    // - twist type M/D
    // - parameter X
    // - sign of X
    // - number of pairs
    // - list of encoded pairs

    // first determine the length of the modulus
    let modulus = curve.q;
    let modulus_length = modulus.clone().to_bytes_be().len();

    let curve_type = vec![BLS12];
    let modulus_len_encoded = vec![modulus_length as u8];
    let modulus_encoded = pad_for_len_be(modulus.clone().to_bytes_be(), modulus_length);

    let a_encoded = pad_for_len_be(curve.a.to_bytes_be(), modulus_length);
    let b_encoded = pad_for_len_be(curve.b.to_bytes_be(), modulus_length);

    let fp2_nonres_encoded = {
        let (mut nonres, is_positive) = curve.non_residue;
        if !is_positive {
            nonres = modulus.clone() - nonres;
        }
        pad_for_len_be(nonres.to_bytes_be(), modulus_length)
    };

    let fp6_nonres_encoded_c0 = {
        let (mut nonres, is_positive) = curve.quadratic_non_residue_0;
        if !is_positive {
            nonres = modulus.clone() - nonres;
        }
        pad_for_len_be(nonres.to_bytes_be(), modulus_length)
    };

    let fp6_nonres_encoded_c1 = {
        let (mut nonres, is_positive) = curve.quadratic_non_residue_1;
        if !is_positive {
            nonres = modulus.clone() - nonres;
        }
        pad_for_len_be(nonres.to_bytes_be(), modulus_length)
    };

    let twist_type = if curve.is_d_type { vec![TWIST_TYPE_D] } else { vec![TWIST_TYPE_M] };

    let (x_decoded, x_is_positive) = curve.x;
    let x_sign = if x_is_positive { vec![0u8] } else { vec![1u8] };
    let x_encoded = x_decoded.to_bytes_be();
    let x_length = vec![x_encoded.len() as u8];

    // now we make two random scalars and do scalar multiplications in G1 and G2 to get pairs that should
    // at the end of the day pair to identity element

    let group_size = curve.r;
    let group_size_encoded = group_size.clone().to_bytes_be();
    let group_size_length = group_size_encoded.len();
    let group_len_encoded = vec![group_size_length as u8];

    // first parse generators
    // g1 generator
    let g1_x = curve.g1_x;
    let g1_x = pad_for_len_be(g1_x.to_bytes_be(), modulus_length);
    let g1_y = curve.g1_y;
    let g1_y = pad_for_len_be(g1_y.to_bytes_be(), modulus_length);

    // g2 generator
    let g2_x_0 = curve.g2_x_0;
    let g2_x_1 = curve.g2_x_1;

    let g2_y_0 = curve.g2_y_0;
    let g2_y_1 = curve.g2_y_1;

    let num_pairs = vec![2u8];

    let mut g1_encodings = vec![];
    let mut g2_encodings = vec![];

    {
        let mut g2_generator_encoding = vec![];
        g2_generator_encoding.extend(pad_for_len_be(g2_x_0.to_bytes_be(), modulus_length));
        g2_generator_encoding.extend(pad_for_len_be(g2_x_1.to_bytes_be(), modulus_length));
        g2_generator_encoding.extend(pad_for_len_be(g2_y_0.to_bytes_be(), modulus_length));
        g2_generator_encoding.extend(pad_for_len_be(g2_y_1.to_bytes_be(), modulus_length));

        g2_encodings.push(g2_generator_encoding.clone());
        g2_encodings.push(g2_generator_encoding.clone());
    }

    // for multiplications we use the public API itself - just construct the corresponding G1
    // multiplication API. Leave G2 as generators for now

    use rand::{Rng, SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

    {
        // - Multiplication API signature
        // - Lengths of modulus (in bytes)
        // - Field modulus
        // - Curve A
        // - Curve B
        // - Length of a scalar field (curve order) (in bytes)
        // - Curve order
        // - X
        // - Y
        // - Scalar
        
        let random_scalar_bytes: Vec<u8> = (0..group_size_length).map(|_| rng.gen()).collect();
        let random_scalar = BigUint::from_bytes_be(&random_scalar_bytes[..]);
        let random_scalar = random_scalar % &group_size;

        let negated_random_scalar = group_size.clone() - &random_scalar;

        let mut common_api_bytes = vec![];
        common_api_bytes.extend_from_slice(&modulus_len_encoded[..]);
        common_api_bytes.extend_from_slice(&modulus_encoded[..]);
        common_api_bytes.extend_from_slice(&a_encoded[..]);
        common_api_bytes.extend_from_slice(&b_encoded[..]);
        common_api_bytes.extend_from_slice(&group_len_encoded[..]);
        common_api_bytes.extend_from_slice(&group_size_encoded[..]);

        let mut mul_calldata = common_api_bytes.clone();
        mul_calldata.extend_from_slice(&g1_x[..]);
        mul_calldata.extend_from_slice(&g1_y[..]);
        mul_calldata.extend(pad_for_len_be(random_scalar.to_bytes_be(), group_size_length));

        let first_g1_encoded = PublicG1Api::mul_point(&mul_calldata[..]).expect("must multiply");

        g1_encodings.push(first_g1_encoded);

        let mut mul_calldata = common_api_bytes.clone();
        mul_calldata.extend_from_slice(&g1_x[..]);
        mul_calldata.extend_from_slice(&g1_y[..]);
        mul_calldata.extend(pad_for_len_be(negated_random_scalar.to_bytes_be(), group_size_length));

        let second_g1_encoded = PublicG1Api::mul_point(&mul_calldata[..]).expect("must multiply");

        g1_encodings.push(second_g1_encoded);
    }

    let mut calldata = vec![];
    calldata.extend(curve_type.into_iter());
    calldata.extend(modulus_len_encoded.into_iter());
    calldata.extend(modulus_encoded.into_iter());
    calldata.extend(a_encoded.into_iter());
    calldata.extend(b_encoded.into_iter());
    calldata.extend(group_len_encoded.into_iter());
    calldata.extend(group_size_encoded.into_iter());
    calldata.extend(fp2_nonres_encoded.into_iter());
    calldata.extend(fp6_nonres_encoded_c0.into_iter());
    calldata.extend(fp6_nonres_encoded_c1.into_iter());
    calldata.extend(twist_type.into_iter());
    calldata.extend(x_length.into_iter());
    calldata.extend(x_encoded.into_iter());
    calldata.extend(x_sign.into_iter());
    calldata.extend(num_pairs.into_iter());
    for (g1, g2) in g1_encodings.into_iter().zip(g2_encodings.into_iter()) {
        calldata.extend(g1.into_iter());
        calldata.extend(g2.into_iter());
    }

    calldata
}

// #[test]
// fn test_single() {
//     let calldata = assemble_single();
//     let result = call_pairing_engine(&calldata[..]);
//     assert!(result.is_ok());

//     let result = result.unwrap()[0];
//     assert!(result == 1u8);
// }

#[test]
fn test_bls12_pairings_from_vectors() {
    let curves = read_dir_and_grab_curves("src/test/test_vectors/bls12/");
    assert!(curves.len() != 0);
    for (curve, _) in curves.into_iter() {
        let calldata = assemble_single_curve_params(curve);
        let result = call_pairing_engine(&calldata[..]);
        assert!(result.is_ok());

        let result = result.unwrap()[0];
        assert!(result == 1u8);
    }
}

extern crate hex;
extern crate csv;

use hex::{encode};
use csv::{Writer};

#[test]
fn dump_pairing_vectors() {
    let curves = read_dir_and_grab_curves("src/test/test_vectors/bls12/");
    assert!(curves.len() != 0);
    let mut writer = Writer::from_path("src/test/test_vectors/bls12/pairing.csv").expect("must open a test file");
    writer.write_record(&["input", "result"]).expect("must write header");
    for (curve, _) in curves.into_iter() {
        let mut input_data = vec![OPERATION_G1_MUL];
        let calldata = assemble_single_curve_params(curve.clone());
        input_data.extend(calldata);
        let expected_result = vec![1u8];
        writer.write_record(&[
            prepend_0x(&encode(&input_data[..])), 
            prepend_0x(&encode(&expected_result[..]))],
        ).expect("must write a record");
    }
    writer.flush().expect("must finalize writing");
}

#[test]
fn dump_fuzzing_vectors() {
    use std::io::Write;
    use std::fs::File;
    let curves = read_dir_and_grab_curves("src/test/test_vectors/bls12/");
    assert!(curves.len() != 0);
    let mut f = File::create("src/test/test_vectors/bls12/pairing.csv").unwrap();
    let mut writer = Writer::from_path("src/test/test_vectors/bls12/pairing.csv").expect("must open a test file");
    writer.write_record(&["input", "result"]).expect("must write header");
    for (curve, _) in curves.into_iter() {
        let mut input_data = vec![OPERATION_G1_MUL];
        let calldata = assemble_single_curve_params(curve.clone());
        input_data.extend(calldata);
        let expected_result = vec![1u8];
        writer.write_record(&[
            prepend_0x(&encode(&input_data[..])), 
            prepend_0x(&encode(&expected_result[..]))],
        ).expect("must write a record");
    }
    writer.flush().expect("must finalize writing");
}

// use rust_test::Bencher;

// #[bench]
// fn bench_single(b: &mut Bencher) {
//     let calldata = assemble_single();
//     b.iter(|| {
//         call_pairing_engine(&calldata[..]).expect("must use");
//     });
// }

fn assemble_single() -> Vec<u8> {
    /// - Curve type
    /// - Lengths of modulus (in bytes)
    /// - Field modulus
    /// - Curve A
    /// - Curve B
    // - non-residue for Fp2
    // - non-residue for Fp6
    // - twist type M/D
    // - parameter X
    // - sign of X
    // - number of pairs
    // - list of encoded pairs

    use num_traits::FromPrimitive;
    use num_integer::Integer;
    use num_traits::Zero;
    use num_traits::Num;
    let modulus_length = 48;
    let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
    let curve_type = vec![1u8];
    let modulus_len_encoded = vec![48u8];
    let modulus_encoded = pad_for_len_be(modulus.clone().to_bytes_be(), modulus_length);
    let a_encoded = pad_for_len_be(BigUint::from(0u64).to_bytes_be(), modulus_length);
    let b_encoded = pad_for_len_be(BigUint::from(4u64).to_bytes_be(), modulus_length);
    let fp2_nonres_encoded = pad_for_len_be((modulus.clone() - BigUint::from(1u64)).to_bytes_be(), modulus_length);
    let fp6_nonres_encoded_c0 = pad_for_len_be(BigUint::from(1u64).to_bytes_be(), modulus_length);
    let fp6_nonres_encoded_c1 = pad_for_len_be(BigUint::from(1u64).to_bytes_be(), modulus_length);
    let twist_type = vec![1u8]; // M
    let x_length = vec![8u8];
    let x_encoded = pad_for_len_be(BigUint::from(0xd201000000010000 as u64).to_bytes_be(), 8); 
    let x_sign = vec![1u8]; // negative x
    let num_pairs = vec![2u8];
    // first pair
    let p_x = BigUint::from_str_radix("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507", 10).unwrap().to_bytes_be();
    let p_y = BigUint::from_str_radix("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569", 10).unwrap().to_bytes_be();
    let mut g1_0_encoding: Vec<u8> = vec![];
    g1_0_encoding.extend(pad_for_len_be(p_x.clone(), modulus_length).into_iter());
    g1_0_encoding.extend(pad_for_len_be(p_y, modulus_length).into_iter());

    let q_x_0 = BigUint::from_str_radix("352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160", 10).unwrap().to_bytes_be();
    let q_x_1 = BigUint::from_str_radix("3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758", 10).unwrap().to_bytes_be();
    let q_y_0 = BigUint::from_str_radix("1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905", 10).unwrap().to_bytes_be();
    let q_y_1 = BigUint::from_str_radix("927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582", 10).unwrap().to_bytes_be();
    let mut g2_0_encoding = vec![];
    g2_0_encoding.extend(pad_for_len_be(q_x_0.clone(), modulus_length).into_iter());
    g2_0_encoding.extend(pad_for_len_be(q_x_1.clone(), modulus_length).into_iter());
    g2_0_encoding.extend(pad_for_len_be(q_y_0.clone(), modulus_length).into_iter());
    g2_0_encoding.extend(pad_for_len_be(q_y_1.clone(), modulus_length).into_iter());

    // second pair 
    let y = modulus.clone() - BigUint::from_str_radix("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569", 10).unwrap();

    let mut g1_1_encoding: Vec<u8> = vec![];
    g1_1_encoding.extend(pad_for_len_be(p_x.clone(), modulus_length).into_iter());
    g1_1_encoding.extend(pad_for_len_be(y.to_bytes_be(), modulus_length).into_iter());

    let g2_1_encoding = g2_0_encoding.clone();

    let mut calldata = vec![];
    calldata.extend(curve_type.into_iter());
    calldata.extend(modulus_len_encoded.into_iter());
    calldata.extend(modulus_encoded.into_iter());
    calldata.extend(a_encoded.into_iter());
    calldata.extend(b_encoded.into_iter());
    calldata.extend(fp2_nonres_encoded.into_iter());
    calldata.extend(fp6_nonres_encoded_c0.into_iter());
    calldata.extend(fp6_nonres_encoded_c1.into_iter());
    calldata.extend(twist_type.into_iter());
    calldata.extend(x_length.into_iter());
    calldata.extend(x_encoded.into_iter());
    calldata.extend(x_sign.into_iter());
    calldata.extend(num_pairs.into_iter());
    calldata.extend(g1_0_encoding.into_iter());
    calldata.extend(g2_0_encoding.into_iter());
    calldata.extend(g1_1_encoding.into_iter());
    calldata.extend(g2_1_encoding.into_iter());

    calldata
}



fn strip_0x(string: &str) -> String {
    let string = string.trim();
    let mut string = string.to_ascii_lowercase().as_bytes().to_vec();
    if string.len() > 2 && string[0..1] == b"0"[..] && string[1..2] == b"x"[..] {
        string = string[2..].to_vec();
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