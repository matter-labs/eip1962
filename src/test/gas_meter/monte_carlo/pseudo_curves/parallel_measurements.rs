use crate::test::gas_meter::bls12;
use crate::test::gas_meter::bn;
use crate::test::gas_meter::mnt4;
use crate::test::gas_meter::mnt6;

use crate::public_interface::API;
use crate::public_interface::constants::*;
use crate::public_interface::sane_limits::*;

use crate::test::parsers::*;

use crate::public_interface::constants::*;

use rand::{Rng, thread_rng};
use rand::distributions::Distribution;
use rand::distributions::Uniform;

use super::gen_params;

#[test]
#[ignore]
fn parallel_measure_one_off_pairing_costs() {
    assert!(crate::features::in_gas_metering());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

    let samples = std::env::var("NUM_SAMPLES").expect("`NUM_SAMPLES` variable must be set");
    let samples : usize  = samples.parse().expect("`NUM_SAMPLES` variable must be an unsigned integer");
    assert!(samples > 0);
    // const SAMPLES: usize = 1_000;
    // const SAMPLES: usize = 2_000;

    use std::thread;

    use std::sync::mpsc::{channel, TryRecvError};

    use rayon::prelude::*;

    let mut bls12_writer = bls12::Bls12ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bls12/one_off_parallel_{}.csv", samples));
    let mut bn_writer = bn::BnReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bn/one_off_parallel_{}.csv", samples));
    let mut mnt4_writer = mnt4::Mnt4ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt4/one_off_parallel_{}.csv", samples));
    let mut mnt6_writer = mnt6::Mnt6ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt6/one_off_parallel_{}.csv", samples));

    let (bls_tx, bls_rx) = channel();
    let (bn_tx, bn_rx) = channel();
    let (mnt4_tx, mnt4_rx) = channel();
    let (mnt6_tx, mnt6_rx) = channel();

    use indicatif::{ProgressBar, ProgressStyle};

    let pb = ProgressBar::new(1u64);

    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
        .progress_chars("##-"));

    let pairs: [usize; 1] = [0];

    // Without subgroup check pairing DOES NOT depend on the group size,
    // BUT for this test and to check it explicitly at least for some part
    // we do a full greedy search

    let mut parameters_space = vec![];
    for num_limbs in NUM_LIMBS_MIN..=NUM_LIMBS_MAX {
        for num_group_limbs in NUM_GROUP_LIMBS_MIN..=NUM_GROUP_LIMBS_MAX {
            parameters_space.push((num_limbs, num_group_limbs, rng.clone(), pb.clone(), (bls_tx.clone(), bn_tx.clone(), mnt4_tx.clone(), mnt6_tx.clone())));
            // parameters_space.push((num_limbs, num_group_limbs, rng.clone(), pb.clone(), (mnt4_tx.clone(), mnt6_tx.clone())));
        }
    }

    drop(bls_tx);
    drop(bn_tx);
    drop(mnt4_tx);
    drop(mnt6_tx);

    pb.set_length((parameters_space.len() * samples) as u64);

    // we are ok with group size beign random here because we DO NOT perform a subgroup checks during gas metering

    let handler = thread::spawn(move || {
        parameters_space.into_par_iter().for_each( |(num_limbs, num_group_limbs, mut rng, pb, (bls_tx, bn_tx, mnt4_tx, mnt6_tx))| {
        // parameters_space.into_par_iter().for_each( |(num_limbs, num_group_limbs, mut rng, pb, (mnt4_tx, mnt6_tx))| {
            for _ in 0..samples {
                {
                    let x_bits = 1;
                    let x_hamming = 1;
                    let curve = gen_params::random_bls12_params(num_limbs, num_group_limbs, &mut rng);
                    for num_pairs in pairs.iter() {
                        let reports = bls12::process_for_curve_and_bit_sizes(curve.clone(), x_bits, x_hamming, *num_pairs);
                        for (r, res_vec) in reports.into_iter() {
                            assert_eq!(res_vec.len(), 1);
                            assert_eq!(res_vec[0], 1u8);
                            bls_tx.send(r).unwrap();
                        }
                    }    
                }
                {
                    // for BN situation is a bit different (cause 6u+2 != 1 always, so miller loop is always non-empty),
                    // so we trim the number of pairs to eliminate the loop and final exponentiation
                    let u_bits = 1;
                    let u_hamming = 1;
                    let curve = gen_params::random_bn_params(num_limbs, num_group_limbs, &mut rng);
                    for num_pairs in pairs.iter() {
                        let reports = bn::process_for_curve_and_bit_sizes(curve.clone(), u_bits, u_hamming, *num_pairs);
                        for (r, res_vec) in reports.into_iter() {
                            assert_eq!(res_vec.len(), 1);
                            assert_eq!(res_vec[0], 1u8);
                            bn_tx.send(r).unwrap();
                        }
                    }
                }
                {
                    let ate_bits = 1;
                    let ate_hamming = 1;

                    let w0_bits = 1;
                    let w0_hamming = 1;

                    let w1_bits = 1;
                    let w1_hamming = 1;

                    let curve = gen_params::random_mnt4_params(num_limbs, num_group_limbs, &mut rng);
                    for num_pairs in pairs.iter() {
                        let reports = mnt4::process_for_curve_and_bit_sizes(
                            curve.clone(), 
                            ate_bits, 
                            ate_hamming,
                            w0_bits,
                            w0_hamming,
                            w1_bits,
                            w1_hamming,                        
                            *num_pairs);
                        for (r, res_vec, _) in reports.into_iter() {
                            assert_eq!(res_vec.len(), 1);
                            assert_eq!(res_vec[0], 1u8);
                            mnt4_tx.send(r).unwrap();
                        }
                    }
                }
                {
                    let ate_bits = 1;
                    let ate_hamming = 1;

                    let w0_bits = 1;
                    let w0_hamming = 1;

                    let w1_bits = 1;
                    let w1_hamming = 1;

                    let curve = gen_params::random_mnt6_params(num_limbs, num_group_limbs, &mut rng);
                    for num_pairs in pairs.iter() {
                        let reports = mnt6::process_for_curve_and_bit_sizes(
                            curve.clone(), 
                            ate_bits, 
                            ate_hamming,
                            w0_bits,
                            w0_hamming,
                            w1_bits,
                            w1_hamming,                        
                            *num_pairs);
                        for (r, res_vec, _) in reports.into_iter() {
                            assert_eq!(res_vec.len(), 1);
                            assert_eq!(res_vec[0], 1u8);
                            mnt6_tx.send(r).unwrap();
                        }
                    }
                }
                pb.inc(1);
            }
        });
    });

    loop {
        let mut all_empty = false;
        let mut all_disconnected = false;
        {
            let subres = bls_rx.try_recv();
            match subres {
                Ok(subres) => {
                    bls12_writer.write_report(subres);
                },
                Err(TryRecvError::Empty) => {
                    all_empty = true;
                },
                Err(TryRecvError::Disconnected) => {
                    all_disconnected = true;
                }
            }
        }
        {
            let subres = bn_rx.try_recv();
            match subres {
                Ok(subres) => {
                    bn_writer.write_report(subres);
                },
                Err(TryRecvError::Empty) => {
                    all_empty = all_empty & true;
                },
                Err(TryRecvError::Disconnected) => {
                    all_disconnected = all_disconnected & true;
                }
            }
        }
        {
            let subres = mnt4_rx.try_recv();
            match subres {
                Ok(subres) => {
                    mnt4_writer.write_report(subres);
                },
                Err(TryRecvError::Empty) => {
                    // all_empty = true;
                    all_empty = all_empty & true;
                },
                Err(TryRecvError::Disconnected) => {
                    // all_disconnected = true;
                    all_disconnected = all_disconnected & true;
                }
            }
        }
        {
            let subres = mnt6_rx.try_recv();
            match subres {
                Ok(subres) => {
                    mnt6_writer.write_report(subres);
                },
                Err(TryRecvError::Empty) => {
                    all_empty = all_empty & true;
                },
                Err(TryRecvError::Disconnected) => {
                    all_disconnected = all_disconnected & true;
                }
            }
        }

        if all_empty {
            std::thread::sleep(std::time::Duration::from_millis(2000u64));
        }

        if all_disconnected {
            pb.println("Joining threads");
            // println!("Joining threads");
            handler.join().unwrap();
            break;
        }
    }

    pb.finish_with_message("Done");
}

#[test]
#[ignore]
fn parallel_measure_bls12_bn_pairing_costs() {
    assert!(crate::features::in_gas_metering());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    let samples = std::env::var("NUM_SAMPLES").expect("`NUM_SAMPLES` variable must be set");
    let samples : usize  = samples.parse().expect("`NUM_SAMPLES` variable must be an unsigned integer");
    assert!(samples > 0);

    use std::thread;

    use std::sync::mpsc::{channel, TryRecvError};

    use rayon::prelude::*;

    let mut bls12_writer = bls12::Bls12ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bls12/miller_loop_and_final_exp_parallel_{}.csv", samples));
    let mut bn_writer = bn::BnReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bn/miller_loop_and_final_exp_parallel_{}.csv", samples));

    let (bls_tx, bls_rx) = channel();
    let (bn_tx, bn_rx) = channel();

    use indicatif::{ProgressBar, ProgressStyle};

    let pb = ProgressBar::new(1u64);

    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
        .progress_chars("##-"));

    // let pairs: [usize; 6] = [2, 4, 6, 8, 12, 16];
    let pairs: [usize; 4] = [2, 4, 8, 16];

    // Without subgroup check pairing DOES NOT depend on the group size,
    // so don't even use it for a sake of smaller noise and just use the maximum one
    let mut parameters_space = vec![];
    for num_limbs in NUM_LIMBS_MIN..=NUM_LIMBS_MAX {
        parameters_space.push((num_limbs, NUM_GROUP_LIMBS_MAX, rng.clone(), pb.clone(), (bls_tx.clone(), bn_tx.clone())));
    }

    drop(bls_tx);
    drop(bn_tx);

    pb.set_length((parameters_space.len() * samples) as u64);

    let handler = thread::spawn(move || {
        parameters_space.into_par_iter().for_each( |(num_limbs, num_group_limbs, mut rng, pb, (bls_tx, bn_tx))| {
            let x_bits_rng = Uniform::new_inclusive(1, MAX_BLS12_X_BIT_LENGTH);
            let u_bits_rng = Uniform::new_inclusive(1, MAX_BN_U_BIT_LENGTH);
            for _ in 0..samples {
                {
                    let x_bits = x_bits_rng.sample(&mut rng);
                    let x_hamming = Uniform::new_inclusive(1, x_bits);
                    let x_hamming = x_hamming.sample(&mut rng);
                    let curve = gen_params::random_bls12_params(num_limbs, num_group_limbs, &mut rng);
                    for num_pairs in pairs.iter() {
                        let reports = bls12::process_for_curve_and_bit_sizes(curve.clone(), x_bits, x_hamming, *num_pairs);
                        for (r, _) in reports.into_iter() {
                            bls_tx.send(r).unwrap();
                        }
                    }     
                }
                {
                    let u_bits = u_bits_rng.sample(&mut rng);
                    let u_hamming = Uniform::new_inclusive(1, u_bits);
                    let u_hamming = u_hamming.sample(&mut rng);
                    let curve = gen_params::random_bn_params(num_limbs, num_group_limbs, &mut rng);
                    for num_pairs in pairs.iter() {
                        let reports = bn::process_for_curve_and_bit_sizes(curve.clone(), u_bits, u_hamming, *num_pairs);
                        for (r, _) in reports.into_iter() {
                            bn_tx.send(r).unwrap();
                        }
                    }
                }
                pb.inc(1);
            }
        });
    });

    loop {
        let mut all_empty = false;
        let mut all_disconnected = false;
        {
            let subres = bls_rx.try_recv();
            match subres {
                Ok(subres) => {
                    bls12_writer.write_report(subres);
                },
                Err(TryRecvError::Empty) => {
                    all_empty = true;
                },
                Err(TryRecvError::Disconnected) => {
                    all_disconnected = true;
                }
            }
        }
        {
            let subres = bn_rx.try_recv();
            match subres {
                Ok(subres) => {
                    bn_writer.write_report(subres);
                },
                Err(TryRecvError::Empty) => {
                    all_empty = all_empty & true;
                },
                Err(TryRecvError::Disconnected) => {
                    all_disconnected = all_disconnected & true;
                }
            }
        }

        if all_empty {
            std::thread::sleep(std::time::Duration::from_millis(2000u64));
        }

        if all_disconnected {
            pb.println("Joining threads");
            handler.join().unwrap();
            break;
        }
    }

    pb.finish_with_message("Done");
}

#[test]
#[ignore]
fn parallel_measure_miller_loop_pairing_costs_mnt() {
    assert!(crate::features::in_gas_metering());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    let samples = std::env::var("NUM_SAMPLES").expect("`NUM_SAMPLES` variable must be set");
    let samples : usize  = samples.parse().expect("`NUM_SAMPLES` variable must be an unsigned integer");
    assert!(samples > 0);

    use std::thread;

    use std::sync::mpsc::{channel, TryRecvError};

    use rayon::prelude::*;

    let mut mnt4_writer = mnt4::Mnt4ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt4/miller_loop_parallel_{}.csv", samples));
    let mut mnt6_writer = mnt6::Mnt6ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt6/miller_loop_parallel_{}.csv", samples));

    let (mnt4_tx, mnt4_rx) = channel();
    let (mnt6_tx, mnt6_rx) = channel();

    use indicatif::{ProgressBar, ProgressStyle};

    let pb = ProgressBar::new(1u64);

    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
        .progress_chars("##-"));

    // let pairs: [usize; 6] = [2, 4, 6, 8, 12, 16];
    let pairs: [usize; 4] = [2, 4, 8, 16];

    // Without subgroup check pairing DOES NOT depend on the group size,
    // so don't even use it for a sake of smaller noise and just use the maximum one

    let mut parameters_space = vec![];
    for num_limbs in NUM_LIMBS_MIN..=NUM_LIMBS_MAX {
        parameters_space.push((num_limbs, NUM_GROUP_LIMBS_MAX, rng.clone(), pb.clone(), (mnt4_tx.clone(), mnt6_tx.clone())));
    }

    drop(mnt4_tx);
    drop(mnt6_tx);

    pb.set_length((parameters_space.len() * samples) as u64);

    let handler = thread::spawn(move || {
        parameters_space.into_par_iter().for_each( |(num_limbs, num_group_limbs, mut rng, pb, (mnt4_tx, mnt6_tx))| {
            let ate_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_ATE_LOOP_COUNT);
            for _ in 0..samples {
                {
                    let ate_bits = ate_rng.sample(&mut rng);
                    let ate_hamming = Uniform::new_inclusive(1, ate_bits);
                    let ate_hamming = ate_hamming.sample(&mut rng);

                    let w0_bits = 1;
                    let w0_hamming = 1;

                    let w1_bits = 1;
                    let w1_hamming = 1;

                    let curve = gen_params::random_mnt4_params(num_limbs, num_group_limbs, &mut rng);
                    for num_pairs in pairs.iter() {
                        let reports = mnt4::process_for_curve_and_bit_sizes(
                            curve.clone(), 
                            ate_bits, 
                            ate_hamming,
                            w0_bits,
                            w0_hamming,
                            w1_bits,
                            w1_hamming,                        
                            *num_pairs);
                        for (r, _, _) in reports.into_iter() {
                            mnt4_tx.send(r).unwrap();
                        }
                    }
                }
                {
                    let ate_bits = ate_rng.sample(&mut rng);
                    let ate_hamming = Uniform::new_inclusive(1, ate_bits);
                    let ate_hamming = ate_hamming.sample(&mut rng);

                    let w0_bits = 1;
                    let w0_hamming = 1;

                    let w1_bits = 1;
                    let w1_hamming = 1;

                    let curve = gen_params::random_mnt6_params(num_limbs, num_group_limbs, &mut rng);
                    for num_pairs in pairs.iter() {
                        let reports = mnt6::process_for_curve_and_bit_sizes(
                            curve.clone(), 
                            ate_bits, 
                            ate_hamming,
                            w0_bits,
                            w0_hamming,
                            w1_bits,
                            w1_hamming,                        
                            *num_pairs);
                        for (r, _, _) in reports.into_iter() {
                            mnt6_tx.send(r).unwrap();
                        }
                    }
                }
                pb.inc(1);
            }
        });
    });

    loop {
        let mut all_empty = false;
        let mut all_disconnected = false;
        {
            let subres = mnt4_rx.try_recv();
            match subres {
                Ok(subres) => {
                    mnt4_writer.write_report(subres);
                },
                Err(TryRecvError::Empty) => {
                    all_empty = true;
                },
                Err(TryRecvError::Disconnected) => {
                    all_disconnected = true;
                }
            }
        }
        {
            let subres = mnt6_rx.try_recv();
            match subres {
                Ok(subres) => {
                    mnt6_writer.write_report(subres);
                },
                Err(TryRecvError::Empty) => {
                    all_empty = all_empty & true;
                },
                Err(TryRecvError::Disconnected) => {
                    all_disconnected = all_disconnected & true;
                }
            }
        }

        if all_empty {
            std::thread::sleep(std::time::Duration::from_millis(2000u64));
        }

        if all_disconnected {
            pb.println("Joining threads");
            handler.join().unwrap();
            break;
        }
    }

    pb.finish_with_message("Done");
}

#[test]
#[ignore]
fn parallel_measure_final_exp_pairing_costs_mnt() {
    assert!(crate::features::in_gas_metering());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    let num_different_bit_length = std::env::var("NUM_BIT_LENGTH").expect("`NUM_BIT_LENGTH` variable must be set");
    let num_different_bit_length : usize  = num_different_bit_length.parse().expect("`NUM_BIT_LENGTH` variable must be an unsigned integer");
    assert!(num_different_bit_length > 0);

    let num_hammings_per_bit_length = std::env::var("NUM_HAMMINGS_PER_BIT_LENGTH").expect("`NUM_HAMMINGS_PER_BIT_LENGTH` variable must be set");
    let num_hammings_per_bit_length : usize  = num_hammings_per_bit_length.parse().expect("`NUM_HAMMINGS_PER_BIT_LENGTH` variable must be an unsigned integer");
    assert!(num_hammings_per_bit_length > 0);

    let samples: usize = num_different_bit_length * num_hammings_per_bit_length;

    use std::thread;

    use std::sync::mpsc::{channel, TryRecvError};

    use rayon::prelude::*;

    let mut mnt4_writer = mnt4::Mnt4ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt4/final_exp_parallel_{}.csv", samples));
    let mut mnt6_writer = mnt6::Mnt6ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt6/final_exp_parallel_{}.csv", samples));

    let (mnt4_tx, mnt4_rx) = channel();
    let (mnt6_tx, mnt6_rx) = channel();

    use indicatif::{ProgressBar, ProgressStyle};

    let pb = ProgressBar::new(1u64);

    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
        .progress_chars("##-"));

    let pairs: [usize; 1] = [2];

    // Without subgroup check pairing DOES NOT depend on the group size,
    // so don't even use it for a sake of smaller noise and just use the maximum one

    let mut parameters_space = vec![];
    for num_limbs in NUM_LIMBS_MIN..=NUM_LIMBS_MAX {
        parameters_space.push((num_limbs, NUM_GROUP_LIMBS_MAX, rng.clone(), pb.clone(), (mnt4_tx.clone(), mnt6_tx.clone())));
    }

    drop(mnt4_tx);
    drop(mnt6_tx);

    pb.set_length((parameters_space.len() * samples * 2) as u64);

    let handler = thread::spawn(move || {
        parameters_space.into_par_iter().for_each( |(num_limbs, num_group_limbs, mut rng, pb, (mnt4_tx, mnt6_tx))| {
            // let ate_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_ATE_LOOP_COUNT);
            let w0_bits_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH);
            let w1_bits_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_FINAL_EXP_W1_BIT_LENGTH);
            for _ in 0..num_different_bit_length {
                {
                let ate_bits = 1;
                let ate_hamming = 1;

                let w0_bits = w0_bits_rng.sample(&mut rng);
                let w0_hamming = Uniform::new_inclusive(1, w0_bits);
                
                let w1_bits = w1_bits_rng.sample(&mut rng);
                let w1_hamming = Uniform::new_inclusive(1, w1_bits);

                for _ in 0..num_hammings_per_bit_length {

                    let w0_hamming = w0_hamming.sample(&mut rng);
                    let w1_hamming = w1_hamming.sample(&mut rng);

                        let curve = gen_params::random_mnt4_params(num_limbs, num_group_limbs, &mut rng);
                        for num_pairs in pairs.iter() {
                            let reports = mnt4::process_for_curve_and_bit_sizes(
                                curve.clone(), 
                                ate_bits, 
                                ate_hamming,
                                w0_bits,
                                w0_hamming,
                                w1_bits,
                                w1_hamming,                        
                                *num_pairs);
                            for (r, _, _) in reports.into_iter() {
                                mnt4_tx.send(r).unwrap();
                            }
                        }
                    }
                    pb.inc(1);
                }
                {
                    let ate_bits = 1;
                    let ate_hamming = 1;
    
                    let w0_bits = w0_bits_rng.sample(&mut rng);
                    let w0_hamming = Uniform::new_inclusive(1, w0_bits);
    
                    let w1_bits = w1_bits_rng.sample(&mut rng);
                    let w1_hamming = Uniform::new_inclusive(1, w1_bits);

                    for _ in 0..num_hammings_per_bit_length {

                        let w0_hamming = w0_hamming.sample(&mut rng);
                        let w1_hamming = w1_hamming.sample(&mut rng);

                        let curve = gen_params::random_mnt6_params(num_limbs, num_group_limbs, &mut rng);
                        for num_pairs in pairs.iter() {
                            let reports = mnt6::process_for_curve_and_bit_sizes(
                                curve.clone(), 
                                ate_bits, 
                                ate_hamming,
                                w0_bits,
                                w0_hamming,
                                w1_bits,
                                w1_hamming,                        
                                *num_pairs);
                            for (r, _, _) in reports.into_iter() {
                                mnt6_tx.send(r).unwrap();
                            }
                        }

                        pb.inc(1);
                    }
                }
                
            }
        });
    });

    loop {
        let mut all_empty = false;
        let mut all_disconnected = false;
        {
            let subres = mnt4_rx.try_recv();
            match subres {
                Ok(subres) => {
                    mnt4_writer.write_report(subres);
                },
                Err(TryRecvError::Empty) => {
                    all_empty = true;
                },
                Err(TryRecvError::Disconnected) => {
                    all_disconnected = true;
                }
            }
        }
        {
            let subres = mnt6_rx.try_recv();
            match subres {
                Ok(subres) => {
                    mnt6_writer.write_report(subres);
                },
                Err(TryRecvError::Empty) => {
                    all_empty = all_empty & true;
                },
                Err(TryRecvError::Disconnected) => {
                    all_disconnected = all_disconnected & true;
                }
            }
        }

        if all_empty {
            std::thread::sleep(std::time::Duration::from_millis(2000u64));
        }

        if all_disconnected {
            pb.println("Joining threads");
            handler.join().unwrap();
            break;
        }
    }

    pb.finish_with_message("Done");
}

// #[test]
// #[ignore]
// fn parallel_check_correspondance_of_gas_metering_mnt() {
//     assert!(crate::features::in_gas_metering());

//     use rand::{SeedableRng};
//     use rand_xorshift::XorShiftRng;

//     let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
//     // const SAMPLES: usize = 1_000;
//     const NUM_BIT_LENGTHS: usize = 10;
//     const NUM_HAMMING_PER_BIT_LENGTH: usize = 1;
//     const SAMPLES: usize = NUM_BIT_LENGTHS * NUM_HAMMING_PER_BIT_LENGTH;

//     use std::thread;

//     use std::sync::mpsc::{channel, TryRecvError};

//     use rayon::prelude::*;

//     let (mnt4_tx, mnt4_rx) = channel();
//     let (mnt6_tx, mnt6_rx) = channel();

//     use indicatif::{ProgressBar, ProgressStyle};

//     let pb = ProgressBar::new(1u64);

//     pb.set_style(ProgressStyle::default_bar()
//         .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
//         .progress_chars("##-"));

//     let pairs: [usize; 4] = [2, 4, 8, 16];

//     let mut parameters_space = vec![];
//     for num_limbs in NUM_LIMBS_MIN..=NUM_LIMBS_MAX {
//         for num_group_limbs in NUM_GROUP_LIMBS_MIN..=NUM_GROUP_LIMBS_MAX {
//             parameters_space.push((num_limbs, num_group_limbs, rng.clone(), pb.clone(), (mnt4_tx.clone(), mnt6_tx.clone())));
//         }
//     }

//     drop(mnt4_tx);
//     drop(mnt6_tx);

//     pb.set_length((parameters_space.len() * SAMPLES * 2) as u64);

//     let handler = thread::spawn(move || {
//         parameters_space.into_par_iter().for_each( |(num_limbs, num_group_limbs, mut rng, pb, (mnt4_tx, mnt6_tx))| {
//             let ate_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_ATE_LOOP_COUNT);
//             let w0_bits_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH);
//             let w1_bits_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_FINAL_EXP_W1_BIT_LENGTH);
//             let mut max_diff_mnt4 = 0i64;
//             let mut max_diff_mnt6 = 0i64;

//             let mut worst_case_report_mnt4 = None;
//             let mut worst_case_report_mnt6 = None;

//             for _ in 0..NUM_BIT_LENGTHS {
//                 {
//                 let ate_bits = ate_rng.sample(&mut rng);
//                 let ate_hamming = Uniform::new_inclusive(1, ate_bits);

//                 let w0_bits = w0_bits_rng.sample(&mut rng);
//                 let w0_hamming = Uniform::new_inclusive(1, w0_bits);
                
//                 let w1_bits = w1_bits_rng.sample(&mut rng);
//                 let w1_hamming = Uniform::new_inclusive(1, w1_bits);

//                 for _ in 0..NUM_HAMMING_PER_BIT_LENGTH {
//                     let ate_hamming = ate_hamming.sample(&mut rng);
//                     let w0_hamming = w0_hamming.sample(&mut rng);
//                     let w1_hamming = w1_hamming.sample(&mut rng);

//                         let curve = gen_params::random_mnt4_params(num_limbs, num_group_limbs, &mut rng);
//                         for num_pairs in pairs.iter() {
//                             let reports = mnt4::process_for_curve_and_bit_sizes(
//                                 curve.clone(), 
//                                 ate_bits, 
//                                 ate_hamming,
//                                 w0_bits,
//                                 w0_hamming,
//                                 w1_bits,
//                                 w1_hamming,                        
//                                 *num_pairs);
//                             for (r, data) in reports.into_iter() {
//                                 let gas_estimated = crate::gas_meter::GasMeter::meter(&data);
//                                 let elapsed = r.run_microseconds;
//                                 if gas_estimated.is_ok() {
//                                     let running_gas = elapsed * GAS_FACTOR_FROM_MICROS;
//                                     let difference = (gas_estimated.unwrap() as i64) - (running_gas as i64);
                    
//                                     if max_diff_mnt4.abs() < difference.abs() {
//                                         max_diff_mnt4 = difference;
//                                         worst_case_report_mnt4 = Some(r.clone());
//                                     }

//                                     mnt4_tx.send(difference).unwrap();
//                                 } else {
//                                     println!("MNT4 gas estimation error {:?}", gas_estimated.err().unwrap());
//                                 }
//                             }
//                         }
//                     }
//                     pb.inc(1);
//                 }
//                 {
//                     let ate_bits = ate_rng.sample(&mut rng);
//                     let ate_hamming = Uniform::new_inclusive(1, ate_bits);
    
//                     let w0_bits = w0_bits_rng.sample(&mut rng);
//                     let w0_hamming = Uniform::new_inclusive(1, w0_bits);
    
//                     let w1_bits = w1_bits_rng.sample(&mut rng);
//                     let w1_hamming = Uniform::new_inclusive(1, w1_bits);

//                     for _ in 0..NUM_HAMMING_PER_BIT_LENGTH {
//                         let ate_hamming = ate_hamming.sample(&mut rng);
//                         let w0_hamming = w0_hamming.sample(&mut rng);
//                         let w1_hamming = w1_hamming.sample(&mut rng);

//                         let curve = gen_params::random_mnt6_params(num_limbs, num_group_limbs, &mut rng);
//                         for num_pairs in pairs.iter() {
//                             let reports = mnt6::process_for_curve_and_bit_sizes(
//                                 curve.clone(), 
//                                 ate_bits, 
//                                 ate_hamming,
//                                 w0_bits,
//                                 w0_hamming,
//                                 w1_bits,
//                                 w1_hamming,                        
//                                 *num_pairs);
//                             for (r, data) in reports.into_iter() {
//                                 let gas_estimated = crate::gas_meter::GasMeter::meter(&data);
//                                 let elapsed = r.run_microseconds;
//                                 if gas_estimated.is_ok() {
//                                     let running_gas = elapsed * GAS_FACTOR_FROM_MICROS;
//                                     let difference = (gas_estimated.unwrap() as i64) - (running_gas as i64);
                    
//                                     if max_diff_mnt6.abs() < difference.abs() {
//                                         max_diff_mnt6 = difference;
//                                         worst_case_report_mnt6 = Some(r.clone());
//                                     }

//                                     mnt6_tx.send(difference).unwrap();
//                                 } else {
//                                     println!("MNT6 gas estimation error {:?}", gas_estimated.err().unwrap());
//                                 }
//                             }
//                         }

//                         pb.inc(1);
//                     }
//                 }

//                 println!("Max MNT4 difference is {} gas for report:\n{:?}", max_diff_mnt4, worst_case_report_mnt4);
//                 println!("Max MNT6 difference is {} gas for report:\n{:?}", max_diff_mnt6, worst_case_report_mnt6);
                
//             }
//         });
//     });


//     let mut max_difference_mnt4 = 0i64;
//     let mut max_difference_mnt6 = 0i64;

//     loop {
//         let mut all_empty = false;
//         let mut all_disconnected = false;
//         {
//             let subres = mnt4_rx.try_recv();
//             match subres {
//                 Ok(subres) => {
//                     if max_difference_mnt4.abs() < subres.abs() {
//                         max_difference_mnt4 = subres;
//                     }
//                 },
//                 Err(TryRecvError::Empty) => {
//                     all_empty = true;
//                 },
//                 Err(TryRecvError::Disconnected) => {
//                     all_disconnected = true;
//                 }
//             }
//         }
//         {
//             let subres = mnt6_rx.try_recv();
//             match subres {
//                 Ok(subres) => {
//                     if max_difference_mnt6.abs() < subres.abs() {
//                         max_difference_mnt6 = subres;
//                     }
//                 },
//                 Err(TryRecvError::Empty) => {
//                     all_empty = all_empty & true;
//                 },
//                 Err(TryRecvError::Disconnected) => {
//                     all_disconnected = all_disconnected & true;
//                 }
//             }
//         }

//         if all_empty {
//             std::thread::sleep(std::time::Duration::from_millis(2000u64));
//         }

//         if all_disconnected {
//             pb.println("Joining threads");
//             handler.join().unwrap();
//             break;
//         }
//     }

//     // pb.println("Joining threads");
//     // handler.join().unwrap();

//     println!("Max difference for MNT4 = {} gas", max_difference_mnt4);
//     println!("Max difference for MNT6 = {} gas", max_difference_mnt6);

//     pb.finish_with_message("Done");
// }
