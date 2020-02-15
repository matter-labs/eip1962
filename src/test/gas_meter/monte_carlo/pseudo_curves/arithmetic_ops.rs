use crate::test::gas_meter::bls12;
use crate::test::gas_meter::bn;
use crate::test::gas_meter::mnt4;
use crate::test::gas_meter::mnt6;
use crate::test::gas_meter::arithmetic_ops;
use crate::test::*;

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
fn run_deterministic_parallel_search_no_filtering() {
    assert!(crate::features::in_gas_metering());

    use crate::public_interface::constants::*;

    use std::thread;

    use std::sync::mpsc::{channel, TryRecvError};

    use rayon::prelude::*;

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

    let samples = std::env::var("NUM_SAMPLES").expect("`NUM_SAMPLES` variable must be set");
    let samples : usize  = samples.parse().expect("`NUM_SAMPLES` variable must be an unsigned integer");
    assert!(samples > 0);

    let mut writer = arithmetic_ops::ArithmeticReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/monte_carlo_arith_deterministic_parallel_unfiltered_{}.csv", samples));

    // let mut multiexp_len = vec![0, 2, 4, 8, 16, 32, 64, 128];
    let mut multiexp_len = vec![2, 4, 8, 16, 32, 64, 128];
    multiexp_len.reverse();

    use indicatif::{ProgressBar, ProgressStyle};

    let pb = ProgressBar::new(1u64);

    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
        .progress_chars("##-"));

    let mut parameters_space = vec![];
    let (tx, rx) = channel();
    for num_limbs in NUM_LIMBS_MIN..=NUM_LIMBS_MAX {
        for num_group_limbs in NUM_GROUP_LIMBS_MIN..=NUM_GROUP_LIMBS_MAX {
            parameters_space.push((num_limbs, num_group_limbs, rng.clone(), pb.clone(), tx.clone()));
        }
    }

    drop(tx);

    pb.set_length((parameters_space.len() * samples) as u64);

    let handler = thread::spawn(move || {
        parameters_space.into_par_iter().for_each(|(num_limbs, num_group_limbs, mut rng, pb, tx)| {
            for _ in 0..samples {
                let (mut curve_ext3, g1_worst_case, g2_worst_case) = gen_params::random_mul_params_a_non_zero_ext3(num_limbs, num_group_limbs, 128, &mut rng);
                let mut curve_ext3_a_zero = curve_ext3.clone();
                make_a_zero_ext3(&mut curve_ext3_a_zero);
                for len in multiexp_len.iter() {
                    trim_multiexp_ext_3(&mut curve_ext3, *len);
                    let reports = arithmetic_ops::process_for_ext3(curve_ext3.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                    for report in reports.into_iter() {
                        tx.send(report).unwrap();
                    }
                }
                for len in multiexp_len.iter() {
                    trim_multiexp_ext_3(&mut curve_ext3_a_zero, *len);
                    let reports = arithmetic_ops::process_for_ext3(curve_ext3_a_zero.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                    for report in reports.into_iter() {
                        tx.send(report).unwrap();
                    }
                }
                let (mut curve_ext2, g1_worst_case, g2_worst_case) = gen_params::random_mul_params_a_non_zero_ext2(num_limbs, num_group_limbs, 128, &mut rng);
                let mut curve_ext2_a_zero = curve_ext2.clone();
                make_a_zero_ext2(&mut curve_ext2_a_zero);
                for len in multiexp_len.iter() {
                    trim_multiexp_ext_2(&mut curve_ext2, *len);
                    let reports = arithmetic_ops::process_for_ext2(curve_ext2.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                    for report in reports.into_iter() {
                        tx.send(report).unwrap();
                    }
                }
                for len in multiexp_len.iter() {
                    trim_multiexp_ext_2(&mut curve_ext2_a_zero, *len);
                    let reports = arithmetic_ops::process_for_ext2(curve_ext2_a_zero.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                    for report in reports.into_iter() {
                        tx.send(report).unwrap();
                    }
                }
            
                pb.inc(1);
            }

        });
    });

    loop {
        let subres = rx.try_recv();
        match subres {
            Ok(r) => {
                writer.write_report(r);
            },
            Err(TryRecvError::Empty) => {
                std::thread::sleep(std::time::Duration::from_millis(1000u64));
            },
            Err(TryRecvError::Disconnected) => {
                handler.join().unwrap();
                break;
            }
        }
    }

    pb.finish();
}

fn run(bytes: &[u8]) {
    let _ = API::run(&bytes);
}