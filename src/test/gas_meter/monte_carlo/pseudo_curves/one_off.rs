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
fn measure_one_off_costs_monte_carlo() {
    assert!(std::option_env!("GAS_METERING").is_some());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    const SAMPLES: usize = 1_000;

    let mut bls12_writer = bls12::Bls12ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bls12/one_off_{}.csv", SAMPLES));
    let mut bn_writer = bn::BnReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bn/one_off_{}.csv", SAMPLES));
    let mut mnt4_writer = mnt4::Mnt4ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt4/one_off_{}.csv", SAMPLES));
    let mut mnt6_writer = mnt6::Mnt6ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt6/one_off_{}.csv", SAMPLES));

    // use pbr::ProgressBar;

    // let mut pb = ProgressBar::new(SAMPLES as u64);
    let pairs = vec![0];

    for num_limbs in 4..=16 {
        for num_group_limbs in 1..=16 {
            {
                let x_bits = 1;
                let x_hamming = 1;
                let curve = gen_params::random_bls12_params(num_limbs, num_group_limbs, &mut rng);
                for num_pairs in pairs.iter() {
                    let reports = bls12::process_for_curve_and_bit_sizes(curve.clone(), x_bits, x_hamming, *num_pairs);
                    for r in reports.into_iter() {
                        bls12_writer.write_report(r);
                    }
                }    
            }
            {
                let u_bits = 1;
                let u_hamming = 1;
                let curve = gen_params::random_bn_params(num_limbs, num_group_limbs, &mut rng);
                for num_pairs in pairs.iter() {
                    let reports = bn::process_for_curve_and_bit_sizes(curve.clone(), u_bits, u_hamming, *num_pairs);
                    for r in reports.into_iter() {
                        bn_writer.write_report(r);
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
                    for r in reports.into_iter() {
                        mnt4_writer.write_report(r);
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
                    for r in reports.into_iter() {
                        mnt6_writer.write_report(r);
                    }
                }
            }
            
        }
    }

    // pb.finish_print("done");
}

#[test]
#[ignore]
fn parallel_measure_one_off_pairing_costs() {
    assert!(std::option_env!("GAS_METERING").is_some());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    // const SAMPLES: usize = 1_000;
    const SAMPLES: usize = 2_000;

    use std::thread;

    use std::sync::mpsc::{channel, TryRecvError};

    use rayon::prelude::*;

    // let mut bls12_writer = bls12::Bls12ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bls12/one_off_parallel_{}.csv", SAMPLES));
    // let mut bn_writer = bn::BnReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bn/one_off_parallel_{}.csv", SAMPLES));
    let mut mnt4_writer = mnt4::Mnt4ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt4/one_off_parallel_{}.csv", SAMPLES));
    let mut mnt6_writer = mnt6::Mnt6ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt6/one_off_parallel_{}.csv", SAMPLES));

    // let (bls_tx, bls_rx) = channel();
    // let (bn_tx, bn_rx) = channel();
    let (mnt4_tx, mnt4_rx) = channel();
    let (mnt6_tx, mnt6_rx) = channel();

    use indicatif::{ProgressBar, ProgressStyle};

    let pb = ProgressBar::new(1u64);

    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
        .progress_chars("##-"));

    let pairs: [usize; 1] = [0];

    let mut parameters_space = vec![];
    for num_limbs in NUM_LIMBS_MIN..=NUM_LIMBS_MAX {
        for num_group_limbs in NUM_GROUP_LIMBS_MIN..=NUM_GROUP_LIMBS_MAX {
            // parameters_space.push((num_limbs, num_group_limbs, rng.clone(), pb.clone(), (bls_tx.clone(), bn_tx.clone(), mnt4_tx.clone(), mnt6_tx.clone())));
            parameters_space.push((num_limbs, num_group_limbs, rng.clone(), pb.clone(), (mnt4_tx.clone(), mnt6_tx.clone())));
        }
    }

    drop(mnt4_tx);
    drop(mnt6_tx);

    pb.set_length((parameters_space.len() * SAMPLES) as u64);

    let handler = thread::spawn(move || {
        // let results: Vec<_> = parameters_space.into_par_iter().map(|(num_limbs, num_group_limbs, rng, pb, tx)| {
        // parameters_space.into_par_iter().for_each( |(num_limbs, num_group_limbs, mut rng, pb, (bls_tx, bn_tx, mnt4_tx, mnt6_tx))| {
        parameters_space.into_par_iter().for_each( |(num_limbs, num_group_limbs, mut rng, pb, (mnt4_tx, mnt6_tx))| {
            for _ in 0..SAMPLES {
                // {
                //     let x_bits = 1;
                //     let x_hamming = 1;
                //     let curve = gen_params::random_bls12_params(num_limbs, num_group_limbs, &mut rng);
                //     for num_pairs in pairs.iter() {
                //         let reports = bls12::process_for_curve_and_bit_sizes(curve.clone(), x_bits, x_hamming, *num_pairs);
                //         for r in reports.into_iter() {
                //             bls_tx.send(r).unwrap();
                //         }
                //     }    
                // }
                // {
                //     let u_bits = 1;
                //     let u_hamming = 1;
                //     let curve = gen_params::random_bn_params(num_limbs, num_group_limbs, &mut rng);
                //     for num_pairs in pairs.iter() {
                //         let reports = bn::process_for_curve_and_bit_sizes(curve.clone(), u_bits, u_hamming, *num_pairs);
                //         for r in reports.into_iter() {
                //             bn_tx.send(r).unwrap();
                //         }
                //     }
                // }
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
                        for r in reports.into_iter() {
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
                        for r in reports.into_iter() {
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
        // {
        //     let subres = bls_rx.try_recv();
        //     match subres {
        //         Ok(subres) => {
        //             bls12_writer.write_report(subres);
        //         },
        //         Err(TryRecvError::Empty) => {
        //             all_empty = true;
        //         },
        //         Err(TryRecvError::Disconnected) => {
        //             all_disconnected = true;
        //         }
        //     }
        // }
        // {
        //     let subres = bn_rx.try_recv();
        //     match subres {
        //         Ok(subres) => {
        //             bn_writer.write_report(subres);
        //         },
        //         Err(TryRecvError::Empty) => {
        //             all_empty = all_empty & true;
        //         },
        //         Err(TryRecvError::Disconnected) => {
        //             all_disconnected = all_disconnected & true;
        //         }
        //     }
        // }
        {
            let subres = mnt4_rx.try_recv();
            match subres {
                Ok(subres) => {
                    mnt4_writer.write_report(subres);
                },
                Err(TryRecvError::Empty) => {
                    all_empty = true;
                    // all_empty = all_empty & true;
                },
                Err(TryRecvError::Disconnected) => {
                    all_disconnected = true;
                    // all_disconnected = all_disconnected & true;
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
fn parallel_measure_miller_loop_pairing_costs() {
    assert!(std::option_env!("GAS_METERING").is_some());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    // const SAMPLES: usize = 1_000;
    const SAMPLES: usize = 2_000;

    use std::thread;

    use std::sync::mpsc::{channel, TryRecvError};

    use rayon::prelude::*;

    let mut mnt4_writer = mnt4::Mnt4ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt4/miller_loop_parallel_{}.csv", SAMPLES));
    let mut mnt6_writer = mnt6::Mnt6ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt6/miller_loop_parallel_{}.csv", SAMPLES));

    let (mnt4_tx, mnt4_rx) = channel();
    let (mnt6_tx, mnt6_rx) = channel();

    use indicatif::{ProgressBar, ProgressStyle};

    let pb = ProgressBar::new(1u64);

    pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
        .progress_chars("##-"));

    let pairs: [usize; 6] = [2, 4, 6, 8, 12, 16];

    let mut parameters_space = vec![];
    for num_limbs in NUM_LIMBS_MIN..=NUM_LIMBS_MAX {
        for num_group_limbs in NUM_GROUP_LIMBS_MIN..=NUM_GROUP_LIMBS_MAX {
            parameters_space.push((num_limbs, num_group_limbs, rng.clone(), pb.clone(), (mnt4_tx.clone(), mnt6_tx.clone())));
        }
    }

    drop(mnt4_tx);
    drop(mnt6_tx);

    pb.set_length((parameters_space.len() * SAMPLES) as u64);

    let handler = thread::spawn(move || {
        parameters_space.into_par_iter().for_each( |(num_limbs, num_group_limbs, mut rng, pb, (mnt4_tx, mnt6_tx))| {
            for _ in 0..SAMPLES {
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
                        for r in reports.into_iter() {
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
                        for r in reports.into_iter() {
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
fn run_single_one_off_test() {
    assert!(std::option_env!("GAS_METERING").is_some());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    const SAMPLES: u64 = 3;

    let pairs = vec![0];

    for num_limbs in 4..=16 {
        for num_group_limbs in 1..=16 {
            {
                let mut total = 0;
                for _ in 0..SAMPLES {
                    let x_bits = 1;
                    let x_hamming = 1;
                    let curve = gen_params::random_bls12_params(num_limbs, num_group_limbs, &mut rng);
                    for num_pairs in pairs.iter() {
                        let reports = bls12::process_for_curve_and_bit_sizes(curve.clone(), x_bits, x_hamming, *num_pairs);
                        for r in reports.into_iter() {
                            total += r.run_microseconds;
                        }
                    }    
                }
                println!("For {} limbs and {} group units run time is {} microseconds", num_limbs, num_group_limbs, total / SAMPLES);
            }
        }
    }

    // pb.finish_print("done");
}


#[test]
fn run_flamer_one_off_bls12() {
    assert!(std::option_env!("GAS_METERING").is_some());
    let input_data = hex::decode("07018076038b8e27878715b7808bb806020707b3254155972e56f81693bc9c25351423d9e4a07f3cdb4b3244da93caaef3acc8ca54d634d37fac5c6ccf3d88427823276888552176b9381ee65718f8df5145447303da4d23f4f6926c326b2c17514a42a68712df9fb7e4fa6d12e90fee38bb9f8ca3dbd453badc0bbfbf654c8582db8d000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000040a5b4f7fe29437a147a961b48ea268671c461b7960cdbdf8ce9546eee7aabc6cfaca779a25b9f8c7fde0273835b9960a00472b8ff574ea6644ee300d676e16ad2eed31b93c281c293d0cd6f7227ea98ed58e656150d94f71906a4a40ce602342abfbbf83b54a00c32cc069f62f05218c26474499681c1276047f296157ae5a6802b062a0f245febdba0fb1811d53b22b4246866eb1318f56a162be6bccc53bf19973f40c12c92b703d7b7848ef8b4466d40823aad3943a312b57432b91ff68be1ec40da9845626890213bd5ba6a1950042ce4c74efbdd6e9d44f1718ea9326e0a4c51f97bcdda93904ae26991b471e9ea942e2b5b8ed26055da11c58bc7b50029275dd19a92ab736627368720908d8cefbd23bfdf11946da6d2a66936ffbb8d94af8035073c4c710bf0c620f9a85476aacc47292a0c78fb25698636583d6058628539a509d6ebed76f43869f695b92873c93d8af7ec96dc6ac30569da0b59f193d1f525d9a532f6d5d2c02f6f56ccb5d482c8f43a077379f5a5c5096b7b62edb7347e3e473726388339492a59b8a6171bbb5c07c88c24c50af7c6b199cebf8d6817c8a0142c4745160f4cf4aa481867b03d12ed5e3c4406ea77332259aa558f753b6e5e140fabead057fa3ed8c4078c130b40c75b538dcfa7a83be55f904e3c0bd3775295f3da0b36d58a5b580b9645c2a75cd2dedc4afd947c4f62a886803d18731740165d75bb73fc4ff360af3e3d00621078267418e1215798be5ba0e3951273402d08b65f7a2a1c7f460b7a2557747e53e775572e4c1fa8f61b9e84e924137a877844be782e3beae417f35eb1e05a20da5676683ec7143da5abbbf0cb608495a4236da79f6e8808dec411badaa592c8e0f4d1940f474ae29e992cb48045ec0101010000");
    let input_data = input_data.unwrap();
    let _ = API::run(&input_data).unwrap();
}

