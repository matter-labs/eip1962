extern crate hex;
extern crate csv;

use crate::{API, PrecompileAPI};
use hex::{decode};
use csv::{Reader};

// #[test]
// fn test_mul_from_csv() {
//     let mut reader = Reader::from_path("src/test/pairings/bls12/test_vectors.csv").expect("must open a test file");
//     for line in reader.records() {
//         let record = line.expect("line must decode");
//         let mut it = record.iter().map(|el| strip_0x_and_pad(el));
//         let err = it.next().unwrap();
//         let field = it.next().unwrap();
//         let a = it.next().unwrap();
//         let b = it.next().unwrap();
//         let order = it.next().unwrap();
//         let point_x = it.next().unwrap();
//         let point_y = it.next().unwrap();
//         let scalar = it.next().unwrap();
//         let result_x = it.next().unwrap();
//         let result_y = it.next().unwrap();
//         // now convert to bytes for API call
//         let err = if err == "1" {
//             true
//         } else {
//             false
//         };
//         let mut encoding = vec![];
//         let field = decode(field).expect("must decode hex");
//         let field_byte_len = field.len() as u8;
//         encoding.push(field_byte_len);
//         encoding.extend(field.into_iter());

//         let a = decode(a).expect("must decode hex");
//         let a = pad_for_len_be(a, field_byte_len as usize);
//         encoding.extend(a.into_iter());
//         let b = decode(b).expect("must decode hex");
//         let b = pad_for_len_be(b, field_byte_len as usize);
//         encoding.extend(b.into_iter());

//         let order = decode(order).expect("must decode hex");
//         let order_byte_len = order.len() as u8;
//         // println!("Order len = {}", order_byte_len);
//         encoding.push(order_byte_len);
//         encoding.extend(order.into_iter());

//         // x
//         let mut coord = decode(point_x).expect("must decode hex");
//         let coord = pad_for_len_be(coord, field_byte_len as usize);
//         encoding.extend(coord.into_iter());

//         // y
//         let mut coord = decode(point_y).expect("must decode hex");
//         let coord = pad_for_len_be(coord, field_byte_len as usize);
//         encoding.extend(coord.into_iter());

//         // scalar 

//         let mut scalar = decode(scalar).expect("must decode hex");
//         let scalar = pad_for_len_be(scalar, order_byte_len as usize);
//         encoding.extend(scalar.into_iter());

//         // println!("Encoding = {}", encode(encoding.clone()));

//         let mut result = vec![];

//         // result_x
//         let mut coord = decode(result_x).expect("must decode hex");
//         let coord = pad_for_len_be(coord, field_byte_len as usize);
//         result.extend(coord.into_iter());

//         // result_y 
//         let mut coord = decode(result_y).expect("must decode hex");
//         let coord = pad_for_len_be(coord, field_byte_len as usize);
//         result.extend(coord.into_iter());

//         // println!("Result encoding = {}", encode(result.clone()));

//         let mul_result = API::mul_point(&encoding[..]).expect("must multiply");

//         // println!("Multiplication result encoding = {}", encode(mul_result.clone()));

//         if result != mul_result {
//             // println!("{:?} != {:?}", result, mul_result);
//             assert!(err);
//         }

//         // assert_eq!(result, mul_result);
//     }

// }

#[test]
fn test_single() {
    let calldata = assemble_single();
    let result = call_bls12_engine(&calldata[..]);
    assert!(result.is_ok());

    let result = result.unwrap()[0];
    assert!(result == 1u8);
}

use rust_test::Bencher;

#[bench]
fn bench_single(b: &mut Bencher) {
    let calldata = assemble_single();
    b.iter(|| {
        call_bls12_engine(&calldata[..]);
    });
}

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

    use num_bigint::BigUint;
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

fn call_bls12_engine(bytes: &[u8]) -> Result<Vec<u8>, ()> {
    use crate::public_interface::{PairingApi, PublicPairingApi};
    PublicPairingApi::pair(&bytes)
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

fn pad_for_len_be(input: Vec<u8>, len: usize) -> Vec<u8> {
    if input.len() < len {
        let mut res = input;
        res.reverse();
        res.resize(len, 0u8);
        res.reverse();
        return res;
    }

    input
}