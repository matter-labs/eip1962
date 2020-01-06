use crate::test::*;
use crate::public_interface::API;
use crate::public_interface::constants::*;
use crate::public_interface::sane_limits::*;
use crate::public_interface::decode_utils::*;

use crate::test::parsers::*;
use crate::test::pairings::bn::*;

use super::*;

pub(crate) struct BnReport {
    six_u_plus_two_bit_length: usize,
    six_u_plus_two_hamming: usize,
    modulus_limbs: usize,
    num_pairs: usize,
    group_limbs: usize,
    x_is_negative: bool,
    x_bit_length: usize,
    x_hamming_weight: usize,
    run_microseconds: u64,
}

extern crate csv;
use std::path::Path;

use csv::{Writer};
use std::fs::File;

pub(crate) struct BnReportWriter {
    writer: Writer<File>
}

impl BnReportWriter {
    pub(crate) fn new_for_path<P: AsRef<Path>>(path: P) -> Self {
        let mut writer = Writer::from_path(path).expect("must open a test file");
        writer.write_record(&["six_u_plus_two_bit_length", 
                            "six_u_plus_two_hamming",
                            "modulus_limbs", 
                            "group_limbs",
                            "num_pairs", 
                            "x_is_negative", 
                            "x_bit_length", 
                            "x_hamming_weight", 
                            "run_microseconds"
                        ]).expect("must write header");
        writer.flush().expect("must finalize writing");

        Self {
            writer
        }
    }

    pub fn write_report(&mut self, report: BnReport) {
        let x_is_negative = if report.x_is_negative {
            "1"
        } else {
            "0"
        };
        self.writer.write_record(&[
            report.six_u_plus_two_bit_length.to_string(),
            report.six_u_plus_two_hamming.to_string(),
            report.modulus_limbs.to_string(),
            report.group_limbs.to_string(),
            report.num_pairs.to_string(),
            x_is_negative.to_owned(),
            report.x_bit_length.to_string(),
            report.x_hamming_weight.to_string(),
            report.run_microseconds.to_string()
            ]
        ).expect("must write a record");

        self.writer.flush().expect("must write to disk");
    } 
}

pub(crate) fn process_for_curve_and_bit_sizes(curve: JsonBnPairingCurveParameters, bits: usize, hamming: usize, num_pairs: usize) -> Vec<BnReport> {
    use std::time::Instant;
    
    let mut reports = vec![];
    
    let new_x = make_x_bit_length_and_hamming_weight(bits, hamming);
    for x_is_negative in vec![true] {
    // for x_is_negative in vec![false, true] {
        let mut new_curve = curve.clone();
        new_curve.x = (new_x.clone(), x_is_negative);
        let (_six_u_plus_two, six_u_plus_two_bit_length, six_u_plus_two_hamming) = six_u_plus_two(&new_x, !x_is_negative);
        let limbs = crate::test::calculate_num_limbs(&new_curve.q).expect("must work");
        let group_order_limbs = crate::test::num_units_for_group_order(&new_curve.r).expect("must work");
        let mut input_data = vec![OPERATION_PAIRING];
        let calldata = assemble_single_curve_params(new_curve, num_pairs);
        if calldata.is_err() {
            // panic!("Bn curve encoding error = {}", calldata.err().unwrap());
            // println!("Bn curve encoding error = {}", calldata.err().unwrap());
            continue
        };
        let calldata = calldata.unwrap();
        input_data.extend(calldata);
        let now = Instant::now();
        let res = API::run(&input_data);
        let elapsed = now.elapsed();
        if res.is_ok() {
            let report = BnReport {
                six_u_plus_two_bit_length: six_u_plus_two_bit_length,
                six_u_plus_two_hamming: six_u_plus_two_hamming,
                modulus_limbs: limbs,
                num_pairs: num_pairs,
                group_limbs: group_order_limbs,
                x_is_negative: x_is_negative,
                x_bit_length: bits,
                x_hamming_weight: hamming,
                run_microseconds: elapsed.as_micros() as u64,
            };

            reports.push(report);
        } else {
            println!("BN error {:?}", res.err().unwrap());
        }
    }

    reports
}

fn process_curve(curve: JsonBnPairingCurveParameters) -> Vec<BnReport> {
    let max_bits = MAX_BN_U_BIT_LENGTH;
    let max_bits = 64;
    let max_hamming = MAX_BN_SIX_U_PLUS_TWO_HAMMING;
    let max_num_pairs = 8;

    let mut reports = vec![];

    for bits in (1..=max_bits).step_by(1) {
        for hamming in (1..=bits).step_by(2) {
            for num_pairs in (2..=max_num_pairs).step_by(2) {
                let subreports = process_for_curve_and_bit_sizes(
                    curve.clone(), bits, hamming, num_pairs
                );
                reports.extend(subreports);
            }
        }
    }

    reports
}

// #[test]
// #[ignore]
// fn test_bench_bn_pairings() {
//     let curves = read_dir_and_grab_curves::<JsonBnPairingCurveParameters>("src/test/test_vectors/bn/");
//     let curves = vec![curves[0].clone()];
//     let mut total_results = vec![];
//     for (curve, _) in curves.into_iter() {
//         let subresult = process_curve(curve);
//         total_results.extend(subresult);
//     }

//     write_reports(total_results, "src/test/gas_meter/bn/reports.csv");
// }

    

