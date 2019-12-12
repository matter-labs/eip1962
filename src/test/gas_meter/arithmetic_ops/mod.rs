use crate::test::*;
use crate::public_interface::API;
use crate::public_interface::constants::*;
use crate::public_interface::sane_limits::*;
use crate::public_interface::decode_utils::*;

use crate::test::parsers::*;
use crate::test::g1_ops::mnt4 as g1_mnt4;
use crate::test::g1_ops::mnt6 as g1_mnt6;

use crate::test::g2_ops::mnt4 as g2_mnt4;
use crate::test::g2_ops::mnt6 as g2_mnt6;

use super::*;

#[derive(Clone, Debug)]
pub(crate) struct ArithmeticReport {
    pub modulus_limbs: usize,
    pub group_limbs: usize,
    pub num_mul_pairs: usize,
    pub a_is_zero: bool,
    pub ext_degree: usize,
    pub run_microseconds_add: u64,
    pub run_microseconds_mul: u64,
    pub run_microseconds_multiexp: u64,
}

pub(crate) fn encode_g1_point(point: (BigUint, BigUint), modulus_length: usize) -> Vec<u8> {
    let (g1_x, g1_y) = point;
    let g1_x = pad_for_len_be(g1_x.to_bytes_be(), modulus_length);
    let g1_y = pad_for_len_be(g1_y.to_bytes_be(), modulus_length);

    let mut res = vec![];
    res.extend(g1_x);
    res.extend(g1_y);

    res                
}

pub(crate) fn encode_g2_point_ext2(point: ((BigUint, BigUint), (BigUint, BigUint)), modulus_length: usize) -> Vec<u8> {
    let ((x_0, x_1), (y_0, y_1)) = point;
    let x_0 = pad_for_len_be(x_0.to_bytes_be(), modulus_length);
    let x_1 = pad_for_len_be(x_1.to_bytes_be(), modulus_length);
    let y_0 = pad_for_len_be(y_0.to_bytes_be(), modulus_length);
    let y_1 = pad_for_len_be(y_1.to_bytes_be(), modulus_length);

    let mut res = vec![];
    res.extend(x_0);
    res.extend(x_1);
    res.extend(y_0);
    res.extend(y_1);

    res                
}

pub(crate) fn encode_g2_point_ext3(point: ((BigUint, BigUint, BigUint), (BigUint, BigUint, BigUint)), modulus_length: usize) -> Vec<u8> {
    let ((x_0, x_1, x_2), (y_0, y_1, y_2)) = point;
    let x_0 = pad_for_len_be(x_0.to_bytes_be(), modulus_length);
    let x_1 = pad_for_len_be(x_1.to_bytes_be(), modulus_length);
    let x_2 = pad_for_len_be(x_2.to_bytes_be(), modulus_length);
    let y_0 = pad_for_len_be(y_0.to_bytes_be(), modulus_length);
    let y_1 = pad_for_len_be(y_1.to_bytes_be(), modulus_length);
    let y_2 = pad_for_len_be(y_2.to_bytes_be(), modulus_length);

    let mut res = vec![];
    res.extend(x_0);
    res.extend(x_1);
    res.extend(x_2);
    res.extend(y_0);
    res.extend(y_1);
    res.extend(y_2);

    res                
}

extern crate csv;
use std::path::Path;

use csv::{Writer};
use std::fs::File;

pub(crate) struct ArithmeticReportWriter {
    writer: Writer<File>
}

#[derive(Clone, Debug)]
pub(crate) struct MaxReportFilter {
    current_max: Option<ArithmeticReport>
}

impl MaxReportFilter {
    pub(crate) fn new() -> Self {
        Self {
            current_max: None
        }
    }

    pub(crate) fn filter(&mut self, report: ArithmeticReport) {
        if self.current_max.is_none() {
            self.current_max = Some(report);
        } else {
            let mut current = self.current_max.take().unwrap();
            assert_eq!(current.modulus_limbs, report.modulus_limbs, "current = {:?}, other = {:?}", current, report);
            assert_eq!(current.group_limbs, report.group_limbs, "current = {:?}, other = {:?}", current, report);
            assert_eq!(current.num_mul_pairs, report.num_mul_pairs, "current = {:?}, other = {:?}", current, report);
            assert_eq!(current.a_is_zero, report.a_is_zero, "current = {:?}, other = {:?}", current, report);
            assert_eq!(current.ext_degree, report.ext_degree, "current = {:?}, other = {:?}", current, report);

            if current.run_microseconds_add < report.run_microseconds_add {
                current.run_microseconds_add = report.run_microseconds_add;
            }

            if current.run_microseconds_mul < report.run_microseconds_mul {
                current.run_microseconds_mul = report.run_microseconds_mul;
            }

            if current.run_microseconds_multiexp < report.run_microseconds_multiexp {
                current.run_microseconds_multiexp = report.run_microseconds_multiexp;
            }

            self.current_max = Some(current);
        }
    }

    pub(crate) fn get(self) -> Option<ArithmeticReport> {
        self.current_max
    }
}

impl ArithmeticReportWriter {
    pub(crate) fn new_for_path<P: AsRef<Path>>(path: P) -> Self {
        let mut writer = Writer::from_path(path).expect("must open a test file");
        writer.write_record(&[
            "modulus_limbs", 
            "group_limbs",
            "num_mul_pairs", 
            "a_is_zero", 
            "ext_degree", 
            "run_microseconds_add",
            "run_microseconds_mul",
            "run_microseconds_multiexp"
        ]).expect("must write header");
        writer.flush().expect("must finalize writing");

        Self {
            writer
        }
    }

    pub fn write_report(&mut self, report: ArithmeticReport) {
        let a_is_zero = if report.a_is_zero {
            "1"
        } else {
            "0"
        };

        self.writer.write_record(&[
            report.modulus_limbs.to_string(),
            report.group_limbs.to_string(),
            report.num_mul_pairs.to_string(),
            a_is_zero.to_owned(),
            report.ext_degree.to_string(),
            report.run_microseconds_add.to_string(),
            report.run_microseconds_mul.to_string(),
            report.run_microseconds_multiexp.to_string()
            ]
        ).expect("must write a record");

        self.writer.flush().expect("must write to disk");
    } 
}

pub(crate) fn process_for_ext2(
    curve: JsonMnt4PairingCurveParameters, 
    g1_worst_case_pair: JsonG1PointScalarMultiplicationPair,
    g2_worst_case_pair: JsonG2PointScalarMultiplicationPair
) -> Vec<ArithmeticReport> {
    use std::time::Instant;
    
    let mut reports = vec![];

    let (curve_a, _) = curve.a.clone();
    
    let a_is_zero = if curve_a.is_zero() && curve.a_twist_0.is_zero() && curve.a_twist_1.is_zero() {
        true
    } else {
        assert!(!curve_a.is_zero() && !curve.a_twist_0.is_zero() && !curve.a_twist_1.is_zero());
        false
    };

    let limbs = calculate_num_limbs(&curve.q).expect("must work");
    let group_order_limbs = crate::test::num_units_for_group_order(&curve.r).expect("must work");
    let (common_g1_data, modulus_length, group_length) = g1_mnt4::assemble_single_curve_params(curve.clone());
    let (common_g2_data, _, _) = g2_mnt4::assemble_single_curve_params(curve.clone());

    let num_mul_pairs_g1 = curve.g1_mul_vectors.len();
    let num_mul_pairs_g2 = curve.g2_mul_vectors.len();

    assert!(num_mul_pairs_g1 >= 2);
    assert!(num_mul_pairs_g2 >= 2);

    let addition_timing_g1 = {
        let mut input_data = vec![OPERATION_G1_ADD];
        input_data.extend(common_g1_data.clone());
        let p0 = encode_g1_point((curve.g1_x.clone(), curve.g1_y.clone()), modulus_length);
        let p1 = encode_g1_point((g1_worst_case_pair.base_x.clone(), g1_worst_case_pair.base_y.clone()), modulus_length);
        input_data.extend(p0);
        input_data.extend(p1);

        let now = Instant::now();
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let multiplication_timing_g1 = {
        let mut input_data = vec![OPERATION_G1_MUL];
        input_data.extend(common_g1_data.clone());
        let (p, _) = g1_mnt4::assemble_single_point_scalar_pair(g1_worst_case_pair, modulus_length, group_length);
        input_data.extend(p);

        let now = Instant::now();
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let multiexp_timing_g1 = {
        let mut input_data = vec![OPERATION_G1_MULTIEXP];
        input_data.extend(common_g1_data.clone());
        input_data.extend(vec![curve.g1_mul_vectors.len() as u8]);
        for pair in curve.g1_mul_vectors.into_iter(){
            let (p, _) = g1_mnt4::assemble_single_point_scalar_pair(pair, modulus_length, group_length);
            input_data.extend(p);
        }

        let now = Instant::now();
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let addition_timing_g2 = {
        let mut input_data = vec![OPERATION_G2_ADD];
        input_data.extend(common_g2_data.clone());
        let p0 = encode_g2_point_ext2(( (curve.g2_x_0.clone(), curve.g2_x_1.clone()), (curve.g2_y_0.clone(), curve.g2_y_1.clone()) ), modulus_length);
        let p1 = encode_g2_point_ext2(( (g2_worst_case_pair.base_x_0.clone(), g2_worst_case_pair.base_x_1.clone()), (g2_worst_case_pair.base_y_0.clone(), g2_worst_case_pair.base_y_1.clone()) ), modulus_length);
        input_data.extend(p0);
        input_data.extend(p1);

        let now = Instant::now();
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let multiplication_timing_g2 = {
        let mut input_data = vec![OPERATION_G2_MUL];
        input_data.extend(common_g2_data.clone());
        let (p, _) = g2_mnt4::assemble_single_point_scalar_pair(g2_worst_case_pair, modulus_length, group_length);
        input_data.extend(p);

        let now = Instant::now();
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let multiexp_timing_g2 = {
        let mut input_data = vec![OPERATION_G2_MULTIEXP];
        input_data.extend(common_g2_data.clone());
        input_data.extend(vec![curve.g2_mul_vectors.len() as u8]);
        for pair in curve.g2_mul_vectors.into_iter(){
            let (p, _) = g2_mnt4::assemble_single_point_scalar_pair(pair, modulus_length, group_length);
            input_data.extend(p);
        }

        let now = Instant::now();
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let report_g1 = ArithmeticReport {
        modulus_limbs: limbs,
        group_limbs: group_order_limbs,
        num_mul_pairs: num_mul_pairs_g1,
        a_is_zero: a_is_zero,
        ext_degree: 1,
        run_microseconds_add: addition_timing_g1.as_micros() as u64,
        run_microseconds_mul: multiplication_timing_g1.as_micros() as u64,
        run_microseconds_multiexp: multiexp_timing_g1.as_micros() as u64,
    };

    reports.push(report_g1);

    let report_g2 = ArithmeticReport {
        modulus_limbs: limbs,
        group_limbs: group_order_limbs,
        num_mul_pairs: num_mul_pairs_g2,
        a_is_zero: a_is_zero,
        ext_degree: 2,
        run_microseconds_add: addition_timing_g2.as_micros() as u64,
        run_microseconds_mul: multiplication_timing_g2.as_micros() as u64,
        run_microseconds_multiexp: multiexp_timing_g2.as_micros() as u64,
    };

    reports.push(report_g2);
    
    reports
}

pub(crate) fn process_for_ext3(
    curve: JsonMnt6PairingCurveParameters, 
    g1_worst_case_pair: JsonG1PointScalarMultiplicationPair,
    g2_worst_case_pair: JsonG2Ext3PointScalarMultiplicationPair
) -> Vec<ArithmeticReport> {
    use std::time::Instant;
    
    let mut reports = vec![];

    let (curve_a, _) = curve.a.clone();
    
    let a_is_zero = if curve_a.is_zero() && curve.a_twist_0.is_zero() && curve.a_twist_1.is_zero() && curve.a_twist_2.is_zero(){
        true
    } else {
        assert!(!curve_a.is_zero() && !curve.a_twist_0.is_zero() && !curve.a_twist_1.is_zero() && !curve.a_twist_2.is_zero());
        false
    };

    let limbs = calculate_num_limbs(&curve.q).expect("must work");
    let group_order_limbs = crate::test::num_units_for_group_order(&curve.r).expect("must work");
    let (common_g1_data, modulus_length, group_length) = g1_mnt6::assemble_single_curve_params(curve.clone());
    let (common_g2_data, _, _) = g2_mnt6::assemble_single_curve_params(curve.clone());

    let num_mul_pairs_g1 = curve.g1_mul_vectors.len();
    let num_mul_pairs_g2 = curve.g2_mul_vectors.len();

    assert!(num_mul_pairs_g1 >= 2);
    assert!(num_mul_pairs_g2 >= 2);

    let addition_timing_g1 = {
        let mut input_data = vec![OPERATION_G1_ADD];
        input_data.extend(common_g1_data.clone());
        let p0 = encode_g1_point((curve.g1_x.clone(), curve.g1_y.clone()), modulus_length);
        let p1 = encode_g1_point((g1_worst_case_pair.base_x.clone(), g1_worst_case_pair.base_y.clone()), modulus_length);
        input_data.extend(p0);
        input_data.extend(p1);

        let now = Instant::now();
        // let _ = API::run(&input_data);
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let multiplication_timing_g1 = {
        let mut input_data = vec![OPERATION_G1_MUL];
        input_data.extend(common_g1_data.clone());
        let (p, _) = g1_mnt6::assemble_single_point_scalar_pair(g1_worst_case_pair, modulus_length, group_length);
        input_data.extend(p);

        let now = Instant::now();
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let multiexp_timing_g1 = {
        let mut input_data = vec![OPERATION_G1_MULTIEXP];
        input_data.extend(common_g1_data.clone());
        input_data.extend(vec![curve.g1_mul_vectors.len() as u8]);
        for pair in curve.g1_mul_vectors.into_iter(){
            let (p, _) = g1_mnt6::assemble_single_point_scalar_pair(pair, modulus_length, group_length);
            input_data.extend(p);
        }

        let now = Instant::now();
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let addition_timing_g2 = {
        let mut input_data = vec![OPERATION_G2_ADD];
        input_data.extend(common_g2_data.clone());
        let p0 = encode_g2_point_ext3(( (curve.g2_x_0.clone(), curve.g2_x_1.clone(), curve.g2_x_2.clone()), (curve.g2_y_0.clone(), curve.g2_y_1.clone(), curve.g2_y_2.clone()) ), modulus_length);
        let p1 = encode_g2_point_ext3(( (g2_worst_case_pair.base_x_0.clone(), g2_worst_case_pair.base_x_1.clone(), g2_worst_case_pair.base_x_2.clone()), (g2_worst_case_pair.base_y_0.clone(), g2_worst_case_pair.base_y_1.clone(), g2_worst_case_pair.base_y_2.clone()) ), modulus_length);
        input_data.extend(p0);
        input_data.extend(p1);

        let now = Instant::now();
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let multiplication_timing_g2 = {
        let mut input_data = vec![OPERATION_G2_MUL];
        input_data.extend(common_g2_data.clone());
        let (p, _) = g2_mnt6::assemble_single_point_scalar_pair(g2_worst_case_pair, modulus_length, group_length);
        input_data.extend(p);

        let now = Instant::now();
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let multiexp_timing_g2 = {
        let mut input_data = vec![OPERATION_G2_MULTIEXP];
        input_data.extend(common_g2_data.clone());
        input_data.extend(vec![curve.g2_mul_vectors.len() as u8]);
        for pair in curve.g2_mul_vectors.into_iter(){
            let (p, _) = g2_mnt6::assemble_single_point_scalar_pair(pair, modulus_length, group_length);
            input_data.extend(p);
        }

        let now = Instant::now();
        let _ = API::run(&input_data).unwrap();
        let elapsed = now.elapsed();

        elapsed
    };

    let report_g1 = ArithmeticReport {
        modulus_limbs: limbs,
        group_limbs: group_order_limbs,
        num_mul_pairs: num_mul_pairs_g1,
        a_is_zero: a_is_zero,
        ext_degree: 1,
        run_microseconds_add: addition_timing_g1.as_micros() as u64,
        run_microseconds_mul: multiplication_timing_g1.as_micros() as u64,
        run_microseconds_multiexp: multiexp_timing_g1.as_micros() as u64,
    };

    reports.push(report_g1);

    let report_g2 = ArithmeticReport {
        modulus_limbs: limbs,
        group_limbs: group_order_limbs,
        num_mul_pairs: num_mul_pairs_g2,
        a_is_zero: a_is_zero,
        ext_degree: 3,
        run_microseconds_add: addition_timing_g2.as_micros() as u64,
        run_microseconds_mul: multiplication_timing_g2.as_micros() as u64,
        run_microseconds_multiexp: multiexp_timing_g2.as_micros() as u64,
    };

    reports.push(report_g2);
    
    reports
}