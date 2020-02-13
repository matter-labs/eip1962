use crate::test::*;
use crate::public_interface::API;
use crate::public_interface::constants::*;
use crate::public_interface::sane_limits::*;
use crate::public_interface::decode_utils::*;

use crate::test::parsers::*;
use crate::test::pairings::mnt6::*;

use super::*;

#[derive(Clone, Debug)]
pub(crate) struct Mnt6Report {
    pub(crate) modulus_limbs: usize,
    pub(crate) num_pairs: usize,
    pub(crate) group_order_limbs: usize,
    pub(crate) x_is_negative: bool,
    pub(crate) x_bit_length: usize,
    pub(crate) x_hamming_weight: usize,
    pub(crate) exp_w0_bit_length: usize,
    pub(crate) exp_w0_hamming: usize,
    pub(crate) exp_w0_is_negative: bool,
    pub(crate) exp_w1_bit_length: usize,
    pub(crate) exp_w1_hamming: usize,
    pub(crate) run_microseconds: u64,
}

extern crate csv;
use std::path::Path;

use csv::{Writer};
use std::fs::File;

pub(crate) struct Mnt6ReportWriter {
    writer: Writer<File>
}

impl Mnt6ReportWriter {
    pub(crate) fn new_for_path<P: AsRef<Path>>(path: P) -> Self {
        let mut writer = Writer::from_path(path).expect("must open a test file");
        writer.write_record(&[
                            "modulus_limbs", 
                            "group_limbs",
                            "num_pairs", 
                            "x_is_negative", 
                            "x_bit_length", 
                            "x_hamming_weight",
                            "exp_w0_bit_length",
                            "exp_w0_hamming",
                            "exp_w0_is_negative",
                            "exp_w1_bit_length",
                            "exp_w1_hamming",
                            "run_microseconds"
                        ]).expect("must write header");
        writer.flush().expect("must finalize writing");

        Self {
            writer
        }
    }

    pub fn write_report(&mut self, report: Mnt6Report) {
        let x_is_negative = if report.x_is_negative {
            "1"
        } else {
            "0"
        };

        let exp_w0_is_negative = if report.exp_w0_is_negative {
            "1"
        } else {
            "0"
        };

        self.writer.write_record(&[
            report.modulus_limbs.to_string(),
            report.group_order_limbs.to_string(),
            report.num_pairs.to_string(),
            x_is_negative.to_owned(),
            report.x_bit_length.to_string(),
            report.x_hamming_weight.to_string(),
            report.exp_w0_bit_length.to_string(),
            report.exp_w0_hamming.to_string(),
            exp_w0_is_negative.to_owned(),
            report.exp_w1_bit_length.to_string(),
            report.exp_w1_hamming.to_string(),
            report.run_microseconds.to_string(),
            ]
        ).expect("must write a record");

        self.writer.flush().expect("must write to disk");
    } 
}

pub(crate) fn process_for_curve_and_bit_sizes(
    curve: JsonMnt6PairingCurveParameters, 
    bits: usize, 
    hamming: usize, 
    w_0_bits: usize,
    w_0_hamming: usize,
    w_1_bits: usize,
    w_1_hamming: usize,
    num_pairs: usize
) -> Vec<(Mnt6Report, Vec<u8>, Vec<u8>)> {
    use std::time::Instant;
    
    let mut reports = vec![];
    
    let new_x = make_x_bit_length_and_hamming_weight(bits, hamming);
    let new_w0 = make_x_bit_length_and_hamming_weight(w_0_bits, w_0_hamming);
    let new_w1 = make_x_bit_length_and_hamming_weight(w_1_bits, w_1_hamming);
    let exp_w0_is_negative = true;
    for x_is_negative in vec![true] {
    // for x_is_negative in vec![false, true] {
        let mut new_curve = curve.clone();
        new_curve.x = (new_x.clone(), x_is_negative);
        new_curve.exp_w0 = (new_w0.clone(), exp_w0_is_negative);
        new_curve.exp_w1 = new_w1.clone();
        let limbs = crate::test::calculate_num_limbs(&new_curve.q).expect("must work");
        let group_order_limbs = crate::test::num_units_for_group_order(&new_curve.r).expect("must work");
        let mut input_data = vec![OPERATION_PAIRING];
        let calldata = assemble_single_curve_params(new_curve, num_pairs, false);
        if calldata.is_err() {
            continue
        };
        let calldata = calldata.unwrap();
        input_data.extend(calldata);
        let now = Instant::now();
        let res = API::run(&input_data);
        let elapsed = now.elapsed();
        if let Ok(res_data) = res {
            let report = Mnt6Report {
                modulus_limbs: limbs,
                group_order_limbs, 
                num_pairs: num_pairs,
                x_is_negative: x_is_negative,
                x_bit_length: bits,
                x_hamming_weight: hamming,
                exp_w0_bit_length: w_0_bits,
                exp_w0_hamming: w_0_hamming,
                exp_w0_is_negative: exp_w0_is_negative,
                exp_w1_bit_length: w_1_bits,
                exp_w1_hamming: w_1_hamming,
                run_microseconds: elapsed.as_micros() as u64,
            };

            reports.push((report, res_data, input_data));
        } else {
            println!("MNT6 error {:?}", res.err().unwrap());
        }
    }

    reports
}

// pub(crate) fn estimate_gas_meter_difference(
//     curve: JsonMnt6PairingCurveParameters, 
//     bits: usize, 
//     hamming: usize, 
//     w_0_bits: usize,
//     w_0_hamming: usize,
//     w_1_bits: usize,
//     w_1_hamming: usize,
//     num_pairs: usize
// ) -> Vec<i64> {
//     use std::time::Instant;

//     let gas_factor = 15u64;
//     let mut reports = vec![];
    
//     let new_x = make_x_bit_length_and_hamming_weight(bits, hamming);
//     let new_w0 = make_x_bit_length_and_hamming_weight(w_0_bits, w_0_hamming);
//     let new_w1 = make_x_bit_length_and_hamming_weight(w_1_bits, w_1_hamming);
//     let exp_w0_is_negative = true;
//     for x_is_negative in vec![true] {
//     // for x_is_negative in vec![false, true] {
//         let mut new_curve = curve.clone();
//         new_curve.x = (new_x.clone(), x_is_negative);
//         new_curve.exp_w0 = (new_w0.clone(), exp_w0_is_negative);
//         new_curve.exp_w1 = new_w1.clone();
//         let mut input_data = vec![OPERATION_PAIRING];
//         let calldata = assemble_single_curve_params(new_curve, num_pairs, false);
//         if calldata.is_err() {
//             continue
//         };
//         let calldata = calldata.unwrap();
//         input_data.extend(calldata);
//         let now = Instant::now();
//         let res = API::run(&input_data);
//         let elapsed = now.elapsed();
//         if res.is_ok() {
//             let gas_estimated = crate::gas_meter::GasMeter::meter(&input_data);
//             if gas_estimated.is_ok() {
//                 let running_gas = (elapsed.as_micros() as u64) * gas_factor;
//                 let difference = (gas_estimated.unwrap() as i64) - (running_gas as i64);

//                 reports.push(difference);
//             } else {
//                 println!("MNT6 gas estimation error {:?}", gas_estimated.err().unwrap());
//                 println!("Data = {}", hex::encode(&input_data));
//             }
//         } else {
//             println!("MNT6 error {:?}", res.err().unwrap());
//         }
//     }

//     reports
// }

