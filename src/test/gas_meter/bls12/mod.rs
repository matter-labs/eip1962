use crate::field::calculate_num_limbs;
use crate::field::biguint_to_u64_vec;
use crate::public_interface::API;
use crate::public_interface::constants::*;
use crate::public_interface::sane_limits::*;

use crate::test::parsers::*;
use crate::test::pairings::bls12::*;

use num_bigint::BigUint;

fn make_x_bit_length_and_hamming_weight(bit_length: usize, hamming_weight: usize) -> BigUint {
    assert!(bit_length > 0);
    assert!(bit_length >= hamming_weight);
    if bit_length == hamming_weight {
        let mut x = BigUint::from(1u64);
        x <<= bit_length;
        x -= BigUint::from(1u64);
        return x;
    }
    let mut x = BigUint::from(1u64);
    x <<= bit_length - 1;
    for i in 1..hamming_weight {
        let mut tmp = BigUint::from(1u64);
        tmp <<= i - 1;
        x += tmp;
    }

    x
}


struct Bls12Report {
    x_bit_length: usize,
    x_hamming_weight: usize,
    modulus_limbs: usize,
    num_pairs: usize,
    x_is_negative: bool,
    run_microseconds: u64,
}

extern crate csv;
use std::path::Path;

use csv::{Writer};

fn write_reports<P: AsRef<Path>>(reports: Vec<Bls12Report>, path: P) {
    assert!(reports.len() != 0);
    let mut writer = Writer::from_path(path).expect("must open a test file");
    writer.write_record(&["x_bit_length", "x_hamming_weight", "modulus_limbs", "num_pairs", "x_is_negative", "run_microseconds"]).expect("must write header");
    for report in reports.into_iter() {
        let x_is_negative = if report.x_is_negative {
            "1"
        } else {
            "0"
        };
        writer.write_record(&[
            report.x_bit_length.to_string(),
            report.x_hamming_weight.to_string(),
            report.modulus_limbs.to_string(),
            report.num_pairs.to_string(),
            x_is_negative.to_owned(),
            report.run_microseconds.to_string()
            ]
        ).expect("must write a record");
    }
    writer.flush().expect("must finalize writing");
}


fn process_curve(curve: JsonBls12PairingCurveParameters) -> Vec<Bls12Report> {
    let max_bits = MAX_BLS12_X_BIT_LENGTH;
    let max_bits = 64;
    let max_hamming = MAX_BLS12_X_HAMMING;
    let max_num_pairs = 8;

    let mut reports = vec![];

    use std::time::Instant;

    for bits in (1..=max_bits).step_by(1) {
        for hamming in (1..=bits).step_by(2) {
            for num_pairs in (2..=max_num_pairs).step_by(2) {
                let new_x = make_x_bit_length_and_hamming_weight(bits, hamming);

                for x_is_negative in vec![false, true] {
                    let mut new_curve = curve.clone();
                    new_curve.x = (new_x.clone(), x_is_negative);
                    let limbs = calculate_num_limbs(&new_curve.q).expect("must work");
                    let mut input_data = vec![OPERATION_PAIRING];
                    let calldata = assemble_single_curve_params(new_curve, num_pairs);
                    input_data.extend(calldata);
                    let now = Instant::now();
                    let res = API::run(&input_data);
                    let elapsed = now.elapsed();
                    assert!(res.is_ok());
                    let report = Bls12Report {
                        x_bit_length: bits,
                        x_hamming_weight: hamming,
                        modulus_limbs: limbs,
                        num_pairs: num_pairs,
                        x_is_negative: x_is_negative,
                        run_microseconds: elapsed.as_micros() as u64,
                    };

                    reports.push(report);
                }
            }
        }
    }

    reports
}

#[test]
fn test_bench_bls12_pairings() {
    let curves = read_dir_and_grab_curves::<JsonBls12PairingCurveParameters>("src/test/test_vectors/bls12/");
    let curves = vec![curves[0].clone()];
    let mut total_results = vec![];
    for (curve, _) in curves.into_iter() {
        let subresult = process_curve(curve);
        total_results.extend(subresult);
    }

    write_reports(total_results, "src/test/gas_meter/bls12/reports.csv");
}

    

