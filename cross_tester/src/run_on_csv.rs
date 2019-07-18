extern crate hex;
extern crate csv;

use self::hex::{decode};
use self::csv::{Reader};

use super::*;

use std::path::Path;

fn cross_check_run_from_csv<P: AsRef<Path>>(path: P) {
    let mut reader = Reader::from_path(path).expect("must open a test file");
    for line in reader.records() {
        let record = line.expect("line must decode");
        let mut it = record.iter().map(|el| strip_0x_and_pad(el));
        let input = decode(strip_0x_and_pad(&it.next().expect("some input"))).expect("must decode");
        let _output = decode(strip_0x_and_pad(&it.next().expect("some output"))).expect("must decode");
        run(&input);
    }
}

#[test]
fn cross_check_g1_muls() {
    let path = "../src/test/test_vectors/bls12/g1_mul.csv";
    cross_check_run_from_csv(path);
}

#[test]
fn cross_check_g2_muls() {
    let path = "../src/test/test_vectors/bls12/g2_mul.csv";
    cross_check_run_from_csv(path);
}

#[test]
fn cross_check_pairings() {
    let path = "../src/test/test_vectors/bls12/pairing.csv";
    cross_check_run_from_csv(path);
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