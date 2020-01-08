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
fn run_arithmetic_ops_pseudo_curves_monte_carlo() {
    assert!(crate::features::in_gas_metering());

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
    assert!(crate::features::in_gas_metering());

    use std::fs::File;

    use crate::public_interface::decode_utils::*;
    use crate::test::parsers::*;
    use crate::test::g1_ops::mnt4 as g1_mnt4;
    use crate::test::g1_ops::mnt6 as g1_mnt6;

    use crate::test::g2_ops::mnt4 as g2_mnt4;
    use crate::test::g2_ops::mnt6 as g2_mnt6;
    use crate::test::gas_meter::arithmetic_ops::*;

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    // low number of samples is ok cause we test many things redundantly
    const SAMPLES: usize = 10_000;

    let limbs_rng = Uniform::new_inclusive(4, 16);
    let group_limbs_rng = Uniform::new_inclusive(1, 16);

    let mut multiexp_len = vec![2, 4, 8, 16, 32, 64, 128];
    multiexp_len.reverse();

    let num_limbs = limbs_rng.sample(&mut rng);
    let num_group_limbs = group_limbs_rng.sample(&mut rng);
    let (mut curve_ext3, g1_worst_case_pair, g2_worst_case_pair) = gen_params::random_mul_params_a_non_zero_ext3(num_limbs, num_group_limbs, 128, &mut rng);
    let mut curve_ext3_a_zero = curve_ext3.clone();
    make_a_zero_ext3(&mut curve_ext3_a_zero);

    let curve = curve_ext3;

    use std::time::Instant;

    let limbs = calculate_num_limbs(&curve.q).expect("must work");
    let group_order_limbs = crate::test::num_units_for_group_order(&curve.r).expect("must work");
    let (common_g1_data, modulus_length, group_length) = g1_mnt6::assemble_single_curve_params(curve.clone());
    let (common_g2_data, _, _) = g2_mnt6::assemble_single_curve_params(curve.clone());

    let num_mul_pairs_g1 = curve.g1_mul_vectors.len();
    let num_mul_pairs_g2 = curve.g2_mul_vectors.len();

    assert!(num_mul_pairs_g1 >= 2);
    assert!(num_mul_pairs_g2 >= 2);

    let addition_timing_g2 = {
        let mut input_data = vec![OPERATION_G2_ADD];
        input_data.extend(common_g2_data.clone());
        let p0 = encode_g2_point_ext3(( (curve.g2_x_0.clone(), curve.g2_x_1.clone(), curve.g2_x_2.clone()), (curve.g2_y_0.clone(), curve.g2_y_1.clone(), curve.g2_y_2.clone()) ), modulus_length);
        let p1 = encode_g2_point_ext3(( (g2_worst_case_pair.base_x_0.clone(), g2_worst_case_pair.base_x_1.clone(), g2_worst_case_pair.base_x_2.clone()), (g2_worst_case_pair.base_y_0.clone(), g2_worst_case_pair.base_y_1.clone(), g2_worst_case_pair.base_y_2.clone()) ), modulus_length);
        input_data.extend(p0);
        input_data.extend(p1);

        let now = Instant::now();
        run(&input_data);
        let elapsed = now.elapsed();

        elapsed
    };
}


#[test]
#[ignore]
fn run_single_curve_and_fields_construction() {
    assert!(crate::features::in_gas_metering());

    use crate::public_interface::decode_utils::*;
    use crate::test::g1_ops::mnt6 as g1_mnt6;

    use crate::test::g2_ops::mnt6 as g2_mnt6;
    use crate::test::gas_meter::arithmetic_ops::*;

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    // low number of samples is ok cause we test many things redundantly
    const SAMPLES: u64 = 1000;

    let mut multiexp_len = vec![2, 4, 8, 16, 32, 64, 128];
    multiexp_len.reverse();

    for num_limbs in 4..=16 {
        let mut results = vec![];
        for num_group_limbs in 1..=16 {
    // for num_limbs in 16..=16 {
        // for num_group_limbs in 16..=16 {
            let (curve_ext3, g1_worst_case_pair, g2_worst_case_pair) = gen_params::random_mul_params_a_non_zero_ext3(num_limbs, num_group_limbs, 128, &mut rng);
            let mut curve_ext3_a_zero = curve_ext3.clone();
            make_a_zero_ext3(&mut curve_ext3_a_zero);

            let curve = curve_ext3;

            use std::time::Instant;

            let limbs = calculate_num_limbs(&curve.q).expect("must work");
            assert!(limbs == num_limbs);
            let group_order_limbs = crate::test::num_units_for_group_order(&curve.r).expect("must work");
            assert!(num_group_limbs == group_order_limbs);
            let (common_g1_data, modulus_length, group_length) = g1_mnt6::assemble_single_curve_params(curve.clone());
            let (common_g2_data, _, _) = g2_mnt6::assemble_single_curve_params(curve.clone());

            let num_mul_pairs_g1 = curve.g1_mul_vectors.len();
            let num_mul_pairs_g2 = curve.g2_mul_vectors.len();

            assert!(num_mul_pairs_g1 >= 2);
            assert!(num_mul_pairs_g2 >= 2);

            let addition_timing_g2 = {
                let mut input_data = vec![OPERATION_G2_ADD];
                input_data.extend(common_g2_data.clone());
                let p0 = encode_g2_point_ext3(( (curve.g2_x_0.clone(), curve.g2_x_1.clone(), curve.g2_x_2.clone()), (curve.g2_y_0.clone(), curve.g2_y_1.clone(), curve.g2_y_2.clone()) ), modulus_length);
                let p1 = encode_g2_point_ext3(( (g2_worst_case_pair.base_x_0.clone(), g2_worst_case_pair.base_x_1.clone(), g2_worst_case_pair.base_x_2.clone()), (g2_worst_case_pair.base_y_0.clone(), g2_worst_case_pair.base_y_1.clone(), g2_worst_case_pair.base_y_2.clone()) ), modulus_length);
                input_data.extend(p0);
                input_data.extend(p1);
                input_data.extend(vec![0u8]);
                // if num_limbs == 16 && num_group_limbs == 13 {
                //     println!("{}", hex::encode(&input_data));
                // }
                
                let mut total = 0;

                for _ in 0..SAMPLES {
                    let now = Instant::now();
                    run(&input_data);
                    let elapsed = now.elapsed();

                    total += elapsed.as_micros();
                }

                (total as u64) / SAMPLES
            };

            results.push(addition_timing_g2);
        }

        let max = *results.iter().max().unwrap();

        println!("G2 Ext3: For {} limbs max curve construction taken {} microseconds", num_limbs, max);
    }
}

#[test]
#[ignore]
fn run_deterministic_search_over_parameter_space_for_g1_and_g2() {
    assert!(crate::features::in_gas_metering());

    use crate::public_interface::constants::*;

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

    let mut writer = arithmetic_ops::ArithmeticReportWriter::new_for_path("src/test/gas_meter/pseudo_curves/monte_carlo_arith_deterministic.csv");

    const RUNS_PER_PARTICULAR_CURVE: usize = 5;
    const RUNS_PER_PARAMETERS_COMBINATION: usize = 40;

    // let mut multiexp_len = vec![0, 2, 4, 8, 16, 32, 64, 128];
    let mut multiexp_len = vec![2, 4, 8, 16, 32, 64, 128];
    multiexp_len.reverse();

    use indicatif::{ProgressBar, ProgressStyle};

    let pb = ProgressBar::new((RUNS_PER_PARAMETERS_COMBINATION * RUNS_PER_PARTICULAR_CURVE * 13 * 16) as u64);

    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
        .progress_chars("##-"));

    for num_limbs in NUM_LIMBS_MIN..=NUM_LIMBS_MAX {
        for num_group_limbs in NUM_GROUP_LIMBS_MIN..=NUM_GROUP_LIMBS_MAX {

            let mut filter_g1 = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];
            let mut filter_g1_a_is_zero = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];
            let mut filter_g2_ext_2 = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];
            let mut filter_g2_ext_3 = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];
            let mut filter_g2_ext_2_a_is_zero = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];
            let mut filter_g2_ext_3_a_is_zero = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];

            for _ in 0..RUNS_PER_PARAMETERS_COMBINATION {
            // while samples_processed < RUNS_PER_PARAMETERS_COMBINATION {
                let (mut curve_ext3, g1_worst_case, g2_worst_case) = gen_params::random_mul_params_a_non_zero_ext3(num_limbs, num_group_limbs, 128, &mut rng);
                let mut curve_ext3_a_zero = curve_ext3.clone();
                make_a_zero_ext3(&mut curve_ext3_a_zero);
                for ((len, filter_g2), filter_g1) in multiexp_len.iter().zip(filter_g2_ext_3.iter_mut()).zip(filter_g1.iter_mut()) {
                    trim_multiexp_ext_3(&mut curve_ext3, *len);
                    for _ in 0..RUNS_PER_PARTICULAR_CURVE {
                        let mut reports = arithmetic_ops::process_for_ext3(curve_ext3.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                        if reports.len() != 2 {
                            break;
                        }
                        // println!("{:?}", filter_g2);
                        filter_g2.filter(reports.pop().unwrap());
                        filter_g1.filter(reports.pop().unwrap());
                    }
                }
                for ((len, filter_g2), filter_g1) in multiexp_len.iter().zip(filter_g2_ext_3_a_is_zero.iter_mut()).zip(filter_g1_a_is_zero.iter_mut()) {
                    trim_multiexp_ext_3(&mut curve_ext3_a_zero, *len);
                    for _ in 0..RUNS_PER_PARTICULAR_CURVE {
                        let mut reports = arithmetic_ops::process_for_ext3(curve_ext3_a_zero.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                        if reports.len() != 2 {
                            break;
                        }
                        // println!("{:?}", filter_g2);
                        filter_g2.filter(reports.pop().unwrap());
                        filter_g1.filter(reports.pop().unwrap());
                    }
                }
                let (mut curve_ext2, g1_worst_case, g2_worst_case) = gen_params::random_mul_params_a_non_zero_ext2(num_limbs, num_group_limbs, 128, &mut rng);
                let mut curve_ext2_a_zero = curve_ext2.clone();
                make_a_zero_ext2(&mut curve_ext2_a_zero);
                for ((len, filter_g2), filter_g1) in multiexp_len.iter().zip(filter_g2_ext_2.iter_mut()).zip(filter_g1.iter_mut()) {
                    trim_multiexp_ext_2(&mut curve_ext2, *len);
                    for _ in 0..RUNS_PER_PARTICULAR_CURVE {
                        let mut reports = arithmetic_ops::process_for_ext2(curve_ext2.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                        if reports.len() != 2 {
                            break;
                        }
                        filter_g2.filter(reports.pop().unwrap());
                        filter_g1.filter(reports.pop().unwrap());
                    }
                }
                for ((len, filter_g2), filter_g1) in multiexp_len.iter().zip(filter_g2_ext_2_a_is_zero.iter_mut()).zip(filter_g1_a_is_zero.iter_mut()) {
                    trim_multiexp_ext_2(&mut curve_ext2_a_zero, *len);
                    for _ in 0..RUNS_PER_PARTICULAR_CURVE {
                        let mut reports = arithmetic_ops::process_for_ext2(curve_ext2_a_zero.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                        if reports.len() != 2 {
                            break;
                        }
                        filter_g2.filter(reports.pop().unwrap());
                        filter_g1.filter(reports.pop().unwrap());
                    }
                }
                
                pb.inc(1);
            }

            for f in vec![filter_g1, filter_g1_a_is_zero, filter_g2_ext_2, filter_g2_ext_2_a_is_zero, filter_g2_ext_3, filter_g2_ext_3_a_is_zero].into_iter() {
                for f in f.into_iter() {
                    let r = f.get().unwrap();
                    writer.write_report(r);
                }
            }
        }
    }

    pb.finish();
}



#[test]
#[ignore]
fn run_deterministic_parallel_search_over_parameter_space_for_g1_and_g2() {
    assert!(crate::features::in_gas_metering());

    use crate::public_interface::constants::*;

    use std::thread;

    use std::sync::mpsc::{channel, TryRecvError};

    use rayon::prelude::*;

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

    let mut writer = arithmetic_ops::ArithmeticReportWriter::new_for_path("src/test/gas_meter/pseudo_curves/monte_carlo_arith_deterministic_parallel.csv");

    const RUNS_PER_PARTICULAR_CURVE: usize = 3;
    const RUNS_PER_PARAMETERS_COMBINATION: usize = 15;

    // let mut multiexp_len = vec![0, 2, 4, 8, 16, 32, 64, 128];
    let mut multiexp_len = vec![2, 4, 8, 16, 32, 64, 128];
    multiexp_len.reverse();

    use indicatif::{ProgressBar, ProgressStyle};

    let pb = ProgressBar::new((RUNS_PER_PARAMETERS_COMBINATION * 13 * 16) as u64);

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

    pb.set_length((parameters_space.len() * RUNS_PER_PARAMETERS_COMBINATION) as u64);

    let handler = thread::spawn(move || {
        parameters_space.into_par_iter().for_each(|(num_limbs, num_group_limbs, mut rng, pb, tx)| {
            let mut filter_g1 = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];
            let mut filter_g1_a_is_zero = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];
            let mut filter_g2_ext_2 = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];
            let mut filter_g2_ext_3 = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];
            let mut filter_g2_ext_2_a_is_zero = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];
            let mut filter_g2_ext_3_a_is_zero = vec![arithmetic_ops::MaxReportFilter::new(); multiexp_len.len()];

            for _ in 0..RUNS_PER_PARAMETERS_COMBINATION {
                let (mut curve_ext3, g1_worst_case, g2_worst_case) = gen_params::random_mul_params_a_non_zero_ext3(num_limbs, num_group_limbs, 128, &mut rng);
                let mut curve_ext3_a_zero = curve_ext3.clone();
                make_a_zero_ext3(&mut curve_ext3_a_zero);
                for ((len, filter_g2), filter_g1) in multiexp_len.iter().zip(filter_g2_ext_3.iter_mut()).zip(filter_g1.iter_mut()) {
                    trim_multiexp_ext_3(&mut curve_ext3, *len);
                    for _ in 0..RUNS_PER_PARTICULAR_CURVE {
                        let mut reports = arithmetic_ops::process_for_ext3(curve_ext3.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                        if reports.len() != 2 {
                            break;
                        }
                        filter_g2.filter(reports.pop().unwrap());
                        filter_g1.filter(reports.pop().unwrap());
                    }
                }
                for ((len, filter_g2), filter_g1) in multiexp_len.iter().zip(filter_g2_ext_3_a_is_zero.iter_mut()).zip(filter_g1_a_is_zero.iter_mut()) {
                    trim_multiexp_ext_3(&mut curve_ext3_a_zero, *len);
                    for _ in 0..RUNS_PER_PARTICULAR_CURVE {
                        let mut reports = arithmetic_ops::process_for_ext3(curve_ext3_a_zero.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                        if reports.len() != 2 {
                            break;
                        }
                        filter_g2.filter(reports.pop().unwrap());
                        filter_g1.filter(reports.pop().unwrap());
                    }
                }
                let (mut curve_ext2, g1_worst_case, g2_worst_case) = gen_params::random_mul_params_a_non_zero_ext2(num_limbs, num_group_limbs, 128, &mut rng);
                let mut curve_ext2_a_zero = curve_ext2.clone();
                make_a_zero_ext2(&mut curve_ext2_a_zero);
                for ((len, filter_g2), filter_g1) in multiexp_len.iter().zip(filter_g2_ext_2.iter_mut()).zip(filter_g1.iter_mut()) {
                    trim_multiexp_ext_2(&mut curve_ext2, *len);
                    for _ in 0..RUNS_PER_PARTICULAR_CURVE {
                        let mut reports = arithmetic_ops::process_for_ext2(curve_ext2.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                        if reports.len() != 2 {
                            break;
                        }
                        filter_g2.filter(reports.pop().unwrap());
                        filter_g1.filter(reports.pop().unwrap());
                    }
                }
                for ((len, filter_g2), filter_g1) in multiexp_len.iter().zip(filter_g2_ext_2_a_is_zero.iter_mut()).zip(filter_g1_a_is_zero.iter_mut()) {
                    trim_multiexp_ext_2(&mut curve_ext2_a_zero, *len);
                    for _ in 0..RUNS_PER_PARTICULAR_CURVE {
                        let mut reports = arithmetic_ops::process_for_ext2(curve_ext2_a_zero.clone(), g1_worst_case.clone(), g2_worst_case.clone());
                        if reports.len() != 2 {
                            break;
                        }
                        filter_g2.filter(reports.pop().unwrap());
                        filter_g1.filter(reports.pop().unwrap());
                    }
                }
            
                pb.inc(1);
            }

            let subresult = vec![filter_g1, filter_g1_a_is_zero, filter_g2_ext_2, filter_g2_ext_2_a_is_zero, filter_g2_ext_3, filter_g2_ext_3_a_is_zero];

            tx.send(subresult).unwrap();
        });
    });

    loop {
        let subres = rx.try_recv();
        match subres {
            Ok(subres) => {
                for f in subres.into_iter() {
                    for f in f.into_iter() {
                        let r = f.get().unwrap();
                        writer.write_report(r);
                    }
                }
            },
            Err(TryRecvError::Empty) => {
                std::thread::sleep(std::time::Duration::from_millis(2000u64));
            },
            Err(TryRecvError::Disconnected) => {
                handler.join().unwrap();
                break;
            }
        }
    }

    pb.finish();
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

    let mut writer = arithmetic_ops::ArithmeticReportWriter::new_for_path("src/test/gas_meter/pseudo_curves/monte_carlo_arith_deterministic_parallel_unfiltered_2.csv");

    const RUNS_PER_PARAMETERS_COMBINATION: usize = 200;

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

    pb.set_length((parameters_space.len() * RUNS_PER_PARAMETERS_COMBINATION) as u64);

    let handler = thread::spawn(move || {
        parameters_space.into_par_iter().for_each(|(num_limbs, num_group_limbs, mut rng, pb, tx)| {
            for _ in 0..RUNS_PER_PARAMETERS_COMBINATION {
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


#[test]
// #[ignore]
fn flame_one_off_field_and_curve_constuction_g2_ext3() {
    assert!(crate::features::in_gas_metering());
    const SAMPLES: usize = 10_000;
    let input_data = hex::decode("04804ab81aafb4512bf2cac1e919dcb7d0e4e00ff21efb5f7c32fbabc8e15d9a47803a98c40e9e1993dec072a0fabf9e894608db26e87be0bbaa0e21f0573d3502599adfa20739dfeb722b7c4124430e22ec5656e9f9c15cd7d164d2d87eea9bf086842298b1a2c1c6de7873b5e21b9f98601e1de014f50143fca69724315a366011033ebbc237d9cb858d87c5297656dd1cea633636d457521a1794ec65996077e75c66293f56f0be1321ec175e1493a471815603ad98582d20bb2e57fe3d3486b2f7c6be8c077b49bdad567923dc1232e497740b76f74ca7942722c9f87474e64addab4b8780e79578576feccf543f1b2b57ffe0635e560ca2bc61aa553794597e743ae6114021b46301172be9e99f805e0b6a3f502445f2e58d4d67787b0792e65d5a78807d688fe9681aaf511825c17760002a3d225ec83ef7b99794cafdbd252e00039fa656529e67bc678cb1e35fa08200b752d3315bd98af3e80266adb51509f598cf0a8d30f98f3d1a427ae7ff3387d23016519699890fdb949d64654de11d001cf24a8dc049a7eae28aecdc046140118b1699b813495756926e4852f60a181daeb655f0d7f7e703103d5b18ad94fb1ba8c20f5a3ab9c313222c75b974d525f114f3c52b0e0f970b21a702e10a8ff95f42f68cff77ce7f6b1febb4c123d1c07388d512cba544c3f23797d4d3c58985913ada374c752fb8d7f69a3d015cfe6f4a6c89becebdf84aa94a4d76b3e50e148312320077fdfc27d9a63d16fda1d77d0845e5d5567f2b1cf50f9a7ec10231b338babdd543c00b95067b8dcb95c7c9b1acde6e8b429ced804c6fda22e375c9a61a9f8b8073d48892b93bf8bfaf7856b2c89e5dd38410a40a6bf0fa71adbea67770cd4f548f4b4b3c5aa92bd39bc8049c440080b8a82a9717120ece870fcc37186d6d0da0aed1610eeef028f9700cd95fe3647494623b98949f06892a3235a4bd190dc08a0874d6df932c09bda59398214afabaab01ada0c7d26c3933abd4d9d6900ccdaa6d74c7685bda6090703804886c0c6f4a507d581922ce7c8f05a01226c7c42c14ef8944cedb4320a7b2cd81141939f966d8cf188c1b7b8329c4cf75334209eb9c0873c7d66891aee4d6693879aade5dbd3aa5ef509e5b5ac09aec967aecc56d298e7386042b5e12164f922ddb636e289c6b6405ab2694622f6fc46a9877faee47dde848c7b697b0fcf37d373280043a021a0ab482dd30c1fa9312907365065ca18e2c9a83f8a32e70e423d9833f02c262f02603ab33049d61c8f04215b8d817690133b482eab6da3e472a27c4660007a09b241374eb6eff5dc9195e99c3046710e30b61b59cd781a2715a0a4c8ee3e4acf3577998307c2cf3dfa9c9dcd6b7f35830ef2114b953ebf244fd8a2e815fbee82e84c9bc44d6a4a9d156804cbf70eb988805f609f1d12c8638c22d4668ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff2d89455495a5a05ab98bfe48a7249f35d8cfeb4dfd47381240a527ae3f165162a280d771ca2600a10f9a9e328d9294d8adcb8b15a0891db5ee339573ebc7e90f18b37444bf627d1dbd28adc73c793a8690f637675308a9911b32ef3f435f496aeca712b26dca206bf852192d5fa38a0dc548802a5a9c9e47c388f7a0d9573ff83136b720c372f2707d27ca661b405f983c53db0352c6e984386cb3760f0b2a2d234e29d63ae1909b507e8dd25cd8cba7dcdfd97e4a95eb91c584cb2a329817765fb903a5b4fac93a6a0024cc0096f4057d3b0a77e7f9c5aa3639ef2192fdcf5dacd377800304e2c28731dc56c60865c09a5709813913577bf96467a406593eef34bef42ddfe1688aefa760dddf648026461758868e421a46bc6765e3ffa8d1c67496daa3fe31349f6db32dd0e1f5b801a84dba1b9f3d858e5a892c2ad98eb16f96fe78bd35d1b08253877a6f2abb24d82c15d4dd3252b02e23189435dbb4a9c867dcccb9a07b98cac4580698b218a486287c028d871bc3f585365a3d9fb4275c10f5e47bdbc9e72c4868de55769f9111e64e2a8dd8a7c0bb2a58cf6cb73cd26391271abf079f9f7881a6c31e92486deefdc9dee42b47fa68f7a34bd72fe8f746ccd9c21cdcf5343aa6e57bb444926e042b448ed05bc4d0dd3ca3cf948d2823c6977ef2a3e5b943c7ddb61d5016a94e24258c9f9b9b81dacf5ce9152db84cb2610e05205917e6777f5e9f84e911d156409ab42edb7a15879b3fc7dbe1aeb8cc6c4723e5ceb90c222f6b2bcf2aa3438007483185416a034962e7fa1570cdd35d64fd38c286cc057d4055ca62423f0a1fab46eb2c2e4f245000d6f4f111ca8e7a98bdd7801d0ec57020e0d726026bb854fdbc601a8c15d7ca89056d366683f1143900e0602b4db49d5696bcfcc58f790a3beb1a5af32dbbfbda49bdb1fc7751a82f367875498981d8e42ba27456574cd745f1fe65adcf0fe8436d3e95b2dd7dc92c7bb7db7546cd11c022bb582f1e792d5dccdda7b57772dba3bb9535d86141d02e947c953a394d25a8afbbe642b240fc035f276b5484dea40dd7bcb23eabef08a6437c0c2ec4ed38f54719690e069727e5850184a89c8fe09aa61e05dc873d4305839c438dd87e411338a245e75214f6fb2b28c06a152c6f720214f1b8cad3eef3ec3c04987da400931be361a33e0494265cbfce9f94c37b099b80aea5e4d59d2394ba43c94ba7f31e565688cfb6887702760534d604297e551bffc9d88cee45663c901160db2492c159b32f0969168f861ecb87f97867b32ceb7cb5de4e30b952e3c6a70f00fdca5e22c68877df5b1861d2f64ef829df9f82dde095c04e2f416616a254ea88496e90fc5eaaffca770c945a10a7cc37f57f4671bc4fd6a8d3e3237a0e9935d8b5956c60a36579c2a4ee09b31875f2fa7cd4b1837a19fc8cca656505987910f74fe5ba97809256d5ff164a0e62ba552e3a14218ffd8a1d7827cc6d614d7beb8658224b48292246d05e2bec73f4944f82c66204e0506b9c51849514747c1c80175ef451fb6bfac8a3599052e5ceb650ce5ad8a6224920c5bf2ce159c6757af9c58d267b78cf8ef354ae265b6725787773c1b1c305ae3516c01e64f0370b99f07673e078fbc8fca8dd46df2df182c185fe1bbcf720682031365ee800ced7a1d1f2d9905f275c4d6c9481c8537438d0a7e58ba450db96a10dd8579d038eb3082b21e5706391a121acd9616a74d34e5c237ede71caca7af38705b74f3ddeec4276f2d865fe61858a5a31df4513adb60a66695039b877feafb9b19245453f9f5240f3e017365a975a92cdbaf8c4664616d8d0b9e90b98c7231538a396880f1cb22867ee3f81c6130e5380160262c5b87fa853492a2b151c1c01ffff3749c65c4edfc3e73fb59c6dbf82d299d114306f519f2e85771baefadf80af458c1a3901eb3a9eb4ba55c4e1b52510c9bce5bf39f4b067f6e0905bb6a203c0a2c88907f4d0227b145a29e84b8f13f1a48fe18ce07d2277eb733b8d1b4ea1aae30dc6e370c36af14e9b8e6556985557a4df767487c014e225e041db68851e6ddaf93188821e3df9566f04a85c608b980e520fb95d11c10c84df4c83f2993d8ee8cf5a260e22d5764ed947485c1634f07db3ab91244c40d6aee5741b2cc5d032e3c55900");
    let input_data = input_data.unwrap();
    for _ in 0..SAMPLES {
        let _ = run(&input_data);
    }
}


fn run(bytes: &[u8]) {
    let _ = API::run(&bytes);
}