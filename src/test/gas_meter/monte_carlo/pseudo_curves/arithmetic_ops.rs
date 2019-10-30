use crate::test::gas_meter::bls12;
use crate::test::gas_meter::bn;
use crate::test::gas_meter::mnt4;
use crate::test::gas_meter::mnt6;
use crate::test::gas_meter::arithmetic_ops;

use crate::public_interface::API;
use crate::public_interface::constants::*;
use crate::public_interface::sane_limits::*;

use crate::test::parsers::*;

use super::gen_params;

use crate::public_interface::constants::*;

use rand::{Rng, thread_rng};
use rand::distributions::Distribution;
use rand::distributions::Uniform;

use num_bigint::BigUint;
use num_traits::Num;
use num_integer::Integer;

fn make_a_zero_ext2(curve: &mut JsonMnt4PairingCurveParameters) {
    curve.a = (BigUint::from(0u64), false);
    curve.a_twist_0 = BigUint::from(0u64);
    curve.a_twist_1 = BigUint::from(0u64);
}

fn make_a_zero_ext3(curve: &mut JsonMnt6PairingCurveParameters) {
    curve.a = (BigUint::from(0u64), false);
    curve.a_twist_0 = BigUint::from(0u64);
    curve.a_twist_1 = BigUint::from(0u64);
    curve.a_twist_2 = BigUint::from(0u64);
}

fn trim_multiexp_ext_2(curve: &mut JsonMnt4PairingCurveParameters, len: usize) {
    curve.g1_mul_vectors.truncate(len);
    curve.g2_mul_vectors.truncate(len);
}

fn trim_multiexp_ext_3(curve: &mut JsonMnt6PairingCurveParameters, len: usize) {
    curve.g1_mul_vectors.truncate(len);
    curve.g2_mul_vectors.truncate(len);
}

#[test]
#[ignore]
fn run_arithmetic_ops_pseudo_curves_monte_carlo() {
    assert!(std::option_env!("GAS_METERING").is_some());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    // low number of samples is ok cause we test many things redundantly
    const SAMPLES: usize = 10_000;

    let mut writer = arithmetic_ops::ArithmeticReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/monte_carlo_arith_{}.csv", SAMPLES));

    let limbs_rng = Uniform::new_inclusive(4, 16);
    let group_limbs_rng = Uniform::new_inclusive(1, 16);

    use pbr::ProgressBar;

    let mut pb = ProgressBar::new(SAMPLES as u64);
    let mut multiexp_len = vec![2, 4, 8, 16, 32, 64, 128];
    multiexp_len.reverse();

    let mut samples_processed = 0;

    while samples_processed < SAMPLES {
        let mut got_results = false;
        let num_limbs = limbs_rng.sample(&mut rng);
        let num_group_limbs = group_limbs_rng.sample(&mut rng);
        let (mut curve_ext3, g1_worst_case, g2_worst_case) = gen_params::random_mul_params_a_non_zero_ext3(num_limbs, num_group_limbs, 128, &mut rng);
        let mut curve_ext3_a_zero = curve_ext3.clone();
        make_a_zero_ext3(&mut curve_ext3_a_zero);
        for len in multiexp_len.iter() {
            trim_multiexp_ext_3(&mut curve_ext3, *len);
            let reports = arithmetic_ops::process_for_ext3(curve_ext3.clone(), g1_worst_case.clone(), g2_worst_case.clone());
            for r in reports.into_iter() {
                got_results = true;
                writer.write_report(r);
            }
        }
        if !got_results {
            continue;
        }
        got_results = false;
        for len in multiexp_len.iter() {
            trim_multiexp_ext_3(&mut curve_ext3_a_zero, *len);
            let reports = arithmetic_ops::process_for_ext3(curve_ext3_a_zero.clone(), g1_worst_case.clone(), g2_worst_case.clone());
            for r in reports.into_iter() {
                got_results = true;
                writer.write_report(r);
            }
        }
        if !got_results {
            continue;
        }
        got_results = false;
        let (mut curve_ext2, g1_worst_case, g2_worst_case) = gen_params::random_mul_params_a_non_zero_ext2(num_limbs, num_group_limbs, 128, &mut rng);
        let mut curve_ext2_a_zero = curve_ext2.clone();
        make_a_zero_ext2(&mut curve_ext2_a_zero);
        for len in multiexp_len.iter() {
            trim_multiexp_ext_2(&mut curve_ext2, *len);
            let reports = arithmetic_ops::process_for_ext2(curve_ext2.clone(), g1_worst_case.clone(), g2_worst_case.clone());
            for r in reports.into_iter() {
                got_results = true;
                writer.write_report(r);
            }
        }
        if !got_results {
            continue;
        }
        got_results = false;
        for len in multiexp_len.iter() {
            trim_multiexp_ext_2(&mut curve_ext2_a_zero, *len);
            let reports = arithmetic_ops::process_for_ext2(curve_ext2_a_zero.clone(), g1_worst_case.clone(), g2_worst_case.clone());
            for r in reports.into_iter() {
                got_results = true;
                writer.write_report(r);
            }
        }
        
        if got_results {
            samples_processed += 1;
            pb.inc();
        }
    }


    pb.finish_print("done");
}


#[test]
// #[ignore]
fn run_single_curve_arithmetic_ops() {
    assert!(std::option_env!("GAS_METERING").is_some());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    // low number of samples is ok cause we test many things redundantly
    const SAMPLES: usize = 10_000;

    let limbs_rng = Uniform::new_inclusive(4, 16);
    let group_limbs_rng = Uniform::new_inclusive(1, 16);

    use pbr::ProgressBar;

    let mut pb = ProgressBar::new(SAMPLES as u64);
    let mut multiexp_len = vec![2, 4, 8, 16, 32, 64, 128];
    multiexp_len.reverse();

    let num_limbs = limbs_rng.sample(&mut rng);
    let num_group_limbs = group_limbs_rng.sample(&mut rng);
    let (mut curve_ext3, g1_worst_case, g2_worst_case) = gen_params::random_mul_params_a_non_zero_ext3(num_limbs, num_group_limbs, 128, &mut rng);
    let mut curve_ext3_a_zero = curve_ext3.clone();
    make_a_zero_ext3(&mut curve_ext3_a_zero);
    for len in multiexp_len.iter() {
        trim_multiexp_ext_3(&mut curve_ext3, *len);
        let reports = arithmetic_ops::process_for_ext3(curve_ext3.clone(), g1_worst_case.clone(), g2_worst_case.clone());
    }
    
    pb.finish_print("done");
}