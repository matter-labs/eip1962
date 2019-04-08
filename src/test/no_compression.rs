extern crate hex;
extern crate csv;

use crate::{ApiImplementation, PrecompileAPI};
use crate::field::{U256Repr, U320Repr};
use hex::{decode};
use csv::{Reader};

#[test]
fn test_mul_from_csv() {
    let mut reader = Reader::from_path("src/test/no_compression_mul.csv").expect("must open a test file");
    for line in reader.records() {
        let record = line.expect("line must decode");
        let mut it = record.iter().map(|el| strip_0x_and_pad(el));
        let err = it.next().unwrap();
        let field = it.next().unwrap();
        let a = it.next().unwrap();
        let b = it.next().unwrap();
        let order = it.next().unwrap();
        let point_x = it.next().unwrap();
        let point_y = it.next().unwrap();
        let scalar = it.next().unwrap();
        let result_x = it.next().unwrap();
        let result_y = it.next().unwrap();
        // now convert to bytes for API call
        let err = if err == "1" {
            true
        } else {
            false
        };
        let mut encoding = vec![];
        let field = decode(field).expect("must decode hex");
        let field_byte_len = field.len() as u8;
        encoding.push(field_byte_len);
        encoding.extend(field.into_iter());

        let a = decode(a).expect("must decode hex");
        let a = pad_for_len_be(a, field_byte_len as usize);
        encoding.extend(a.into_iter());
        let b = decode(b).expect("must decode hex");
        let b = pad_for_len_be(b, field_byte_len as usize);
        encoding.extend(b.into_iter());

        let order = decode(order).expect("must decode hex");
        let order_byte_len = order.len() as u8;
        // println!("Order len = {}", order_byte_len);
        encoding.push(order_byte_len);
        encoding.extend(order.into_iter());

        // x
        let mut coord = decode(point_x).expect("must decode hex");
        let coord = pad_for_len_be(coord, field_byte_len as usize);
        encoding.extend(coord.into_iter());

        // y
        let mut coord = decode(point_y).expect("must decode hex");
        let coord = pad_for_len_be(coord, field_byte_len as usize);
        encoding.extend(coord.into_iter());

        // scalar 

        let mut scalar = decode(scalar).expect("must decode hex");
        let scalar = pad_for_len_be(scalar, order_byte_len as usize);
        encoding.extend(scalar.into_iter());

        // println!("Encoding = {}", encode(encoding.clone()));

        let mut result = vec![];

        // result_x
        let mut coord = decode(result_x).expect("must decode hex");
        let coord = pad_for_len_be(coord, field_byte_len as usize);
        result.extend(coord.into_iter());

        // result_y 
        let mut coord = decode(result_y).expect("must decode hex");
        let coord = pad_for_len_be(coord, field_byte_len as usize);
        result.extend(coord.into_iter());

        // println!("Result encoding = {}", encode(result.clone()));

        let mul_result = ApiImplementation::<U256Repr, U256Repr>::mul_point(&encoding[..]).expect("must multiply");

        // println!("Multiplication result encoding = {}", encode(mul_result.clone()));

        if result != mul_result {
            assert!(err);
            // println!("{:?} != {:?}", result, mul_result);
        }

        // assert_eq!(result, mul_result);
    }

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