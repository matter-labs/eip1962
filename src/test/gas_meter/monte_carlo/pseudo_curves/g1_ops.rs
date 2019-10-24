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

mod gen_params;

#[test]
#[ignore]
fn run_g1_addition_pseudo_curves_monte_carlo() {
    assert!(std::option_env!("GAS_METERING").is_some());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    const SAMPLES: usize = 100_000;
    // const SAMPLES: usize = 1_000_000;

    let mut writer = bls12::Bls12ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bls12/monte_carlo_f_exp_{}.csv", SAMPLES));
    let mut bn_writer = bn::BnReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bn/monte_carlo_f_exp_{}.csv", SAMPLES));
    let mut mnt4_writer = mnt4::Mnt4ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt4/monte_carlo_f_exp_{}.csv", SAMPLES));
    let mut mnt6_writer = mnt6::Mnt6ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt6/monte_carlo_f_exp_{}.csv", SAMPLES));

    let x_bits_rng = Uniform::new_inclusive(1, MAX_BLS12_X_BIT_LENGTH);
    let u_bits_rng = Uniform::new_inclusive(1, MAX_BN_U_BIT_LENGTH);
    let ate_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_ATE_LOOP_COUNT);
    let w0_bits_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH);
    let w1_bits_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_FINAL_EXP_W1_BIT_LENGTH);
    let limbs_rng = Uniform::new_inclusive(4, 16);
    let group_limbs_rng = Uniform::new_inclusive(1, 16);

    use pbr::ProgressBar;

    let mut pb = ProgressBar::new(SAMPLES as u64);
    let pairs = vec![2, 4, 6];

    let mut samples_processed = 0;

    while samples_processed < SAMPLES {
        let curve_type = curve_rng.sample(&mut rng);
        let mut got_results = false;
        match curve_type {
            BLS12 => {
                let x_bits = x_bits_rng.sample(&mut rng);
                let x_hamming = Uniform::new_inclusive(1, x_bits);
                let x_hamming = x_hamming.sample(&mut rng);
                let num_limbs = limbs_rng.sample(&mut rng);
                let num_group_limbs = group_limbs_rng.sample(&mut rng);
                let curve = gen_params::random_bls12_params(num_limbs, num_group_limbs, &mut rng);
                for num_pairs in pairs.iter() {
                    let reports = bls12::process_for_curve_and_bit_sizes(curve.clone(), x_bits, x_hamming, *num_pairs);
                    for r in reports.into_iter() {
                        got_results = true;
                        bls12_writer.write_report(r);
                    }
                }                
            },
            BN => {
                let u_bits = u_bits_rng.sample(&mut rng);
                let u_hamming = Uniform::new_inclusive(1, u_bits);
                let u_hamming = u_hamming.sample(&mut rng);
                let num_limbs = limbs_rng.sample(&mut rng);
                let num_group_limbs = group_limbs_rng.sample(&mut rng);
                let curve = gen_params::random_bn_params(num_limbs, num_group_limbs, &mut rng);
                for num_pairs in pairs.iter() {
                    let reports = bn::process_for_curve_and_bit_sizes(curve.clone(), u_bits, u_hamming, *num_pairs);
                    for r in reports.into_iter() {
                        got_results = true;
                        bn_writer.write_report(r);
                    }
                }
            },
            MNT4 => {
                let ate_bits = ate_rng.sample(&mut rng);
                let ate_hamming = Uniform::new_inclusive(1, ate_bits);
                let ate_hamming = ate_hamming.sample(&mut rng);

                let w0_bits = w0_bits_rng.sample(&mut rng);
                let w0_hamming = Uniform::new_inclusive(1, w0_bits);
                let w0_hamming = w0_hamming.sample(&mut rng);

                let w1_bits = w1_bits_rng.sample(&mut rng);
                let w1_hamming = Uniform::new_inclusive(1, w1_bits);
                let w1_hamming = w1_hamming.sample(&mut rng);

                let num_limbs = limbs_rng.sample(&mut rng);
                let num_group_limbs = group_limbs_rng.sample(&mut rng);
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
                        got_results = true;
                        mnt4_writer.write_report(r);
                    }
                }
            },
            MNT6 => {
                let ate_bits = ate_rng.sample(&mut rng);
                let ate_hamming = Uniform::new_inclusive(1, ate_bits);
                let ate_hamming = ate_hamming.sample(&mut rng);

                let w0_bits = w0_bits_rng.sample(&mut rng);
                let w0_hamming = Uniform::new_inclusive(1, w0_bits);
                let w0_hamming = w0_hamming.sample(&mut rng);

                let w1_bits = w1_bits_rng.sample(&mut rng);
                let w1_hamming = Uniform::new_inclusive(1, w1_bits);
                let w1_hamming = w1_hamming.sample(&mut rng);

                let num_limbs = limbs_rng.sample(&mut rng);
                let num_group_limbs = group_limbs_rng.sample(&mut rng);
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
                        got_results = true;
                        mnt6_writer.write_report(r);
                    }
                }
            },
            _ => {
                unimplemented!();
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
#[ignore]
fn run_bls12_bn_pseudo_curves_monte_carlo() {
    assert!(std::option_env!("GAS_METERING").is_some());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    const SAMPLES: usize = 100_000;
    // const SAMPLES: usize = 1_000_000;

    let curve_rng = Uniform::new_inclusive(BLS12, BN);

    let mut bls12_writer = bls12::Bls12ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bls12/monte_carlo_f_exp_{}.csv", SAMPLES));
    let mut bn_writer = bn::BnReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/bn/monte_carlo_f_exp_{}.csv", SAMPLES));

    let x_bits_rng = Uniform::new_inclusive(1, MAX_BLS12_X_BIT_LENGTH);
    let u_bits_rng = Uniform::new_inclusive(1, MAX_BN_U_BIT_LENGTH);
    let limbs_rng = Uniform::new_inclusive(4, 16);
    let group_limbs_rng = Uniform::new_inclusive(1, 16);

    use pbr::ProgressBar;

    let mut pb = ProgressBar::new(SAMPLES as u64);
    let pairs = vec![2, 4, 6];

    let mut samples_processed = 0;

    while samples_processed < SAMPLES {
        let curve_type = curve_rng.sample(&mut rng);
        let mut got_results = false;
        match curve_type {
            BLS12 => {
                let x_bits = x_bits_rng.sample(&mut rng);
                let x_hamming = Uniform::new_inclusive(1, x_bits);
                let x_hamming = x_hamming.sample(&mut rng);
                let num_limbs = limbs_rng.sample(&mut rng);
                let num_group_limbs = group_limbs_rng.sample(&mut rng);
                let curve = gen_params::random_bls12_params(num_limbs, num_group_limbs, &mut rng);
                for num_pairs in pairs.iter() {
                    let reports = bls12::process_for_curve_and_bit_sizes(curve.clone(), x_bits, x_hamming, *num_pairs);
                    for r in reports.into_iter() {
                        got_results = true;
                        bls12_writer.write_report(r);
                    }
                }                
            },
            BN => {
                let u_bits = u_bits_rng.sample(&mut rng);
                let u_hamming = Uniform::new_inclusive(1, u_bits);
                let u_hamming = u_hamming.sample(&mut rng);
                let num_limbs = limbs_rng.sample(&mut rng);
                let num_group_limbs = group_limbs_rng.sample(&mut rng);
                let curve = gen_params::random_bn_params(num_limbs, num_group_limbs, &mut rng);
                for num_pairs in pairs.iter() {
                    let reports = bn::process_for_curve_and_bit_sizes(curve.clone(), u_bits, u_hamming, *num_pairs);
                    for r in reports.into_iter() {
                        got_results = true;
                        bn_writer.write_report(r);
                    }
                }
            },
            _ => {
                unreachable!();
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
#[ignore]
fn run_mnt_pseudo_curves_monte_carlo() {
    assert!(std::option_env!("GAS_METERING").is_some());

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    const SAMPLES: usize = 100_000;
    // const SAMPLES: usize = 1_000_000;

    let curve_rng = Uniform::new_inclusive(MNT4, MNT6);

    let mut mnt4_writer = mnt4::Mnt4ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt4/monte_carlo_f_exp_{}.csv", SAMPLES));
    let mut mnt6_writer = mnt6::Mnt6ReportWriter::new_for_path(format!("src/test/gas_meter/pseudo_curves/mnt6/monte_carlo_f_exp_{}.csv", SAMPLES));

    let ate_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_ATE_LOOP_COUNT);
    let w0_bits_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH);
    let w1_bits_rng = Uniform::new_inclusive(1, MAX_ATE_PAIRING_FINAL_EXP_W1_BIT_LENGTH);
    let limbs_rng = Uniform::new_inclusive(4, 16);
    let group_limbs_rng = Uniform::new_inclusive(1, 16);

    use pbr::ProgressBar;

    let mut pb = ProgressBar::new(SAMPLES as u64);
    // unfortunately for MNT we need more sample points
    let pairs = vec![2, 4, 8, 16, 24, 32];

    let mut samples_processed = 0;

    while samples_processed < SAMPLES {
        let curve_type = curve_rng.sample(&mut rng);
        let mut got_results = false;
        match curve_type {
            MNT4 => {
                let ate_bits = ate_rng.sample(&mut rng);
                let ate_hamming = Uniform::new_inclusive(1, ate_bits);
                let ate_hamming = ate_hamming.sample(&mut rng);

                let w0_bits = w0_bits_rng.sample(&mut rng);
                let w0_hamming = Uniform::new_inclusive(1, w0_bits);
                let w0_hamming = w0_hamming.sample(&mut rng);

                let w1_bits = w1_bits_rng.sample(&mut rng);
                let w1_hamming = Uniform::new_inclusive(1, w1_bits);
                let w1_hamming = w1_hamming.sample(&mut rng);

                let num_limbs = limbs_rng.sample(&mut rng);
                let num_group_limbs = group_limbs_rng.sample(&mut rng);
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
                        got_results = true;
                        mnt4_writer.write_report(r);
                    }
                }
            },
            MNT6 => {
                let ate_bits = ate_rng.sample(&mut rng);
                let ate_hamming = Uniform::new_inclusive(1, ate_bits);
                let ate_hamming = ate_hamming.sample(&mut rng);

                let w0_bits = w0_bits_rng.sample(&mut rng);
                let w0_hamming = Uniform::new_inclusive(1, w0_bits);
                let w0_hamming = w0_hamming.sample(&mut rng);

                let w1_bits = w1_bits_rng.sample(&mut rng);
                let w1_hamming = Uniform::new_inclusive(1, w1_bits);
                let w1_hamming = w1_hamming.sample(&mut rng);

                let num_limbs = limbs_rng.sample(&mut rng);
                let num_group_limbs = group_limbs_rng.sample(&mut rng);
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
                        got_results = true;
                        mnt6_writer.write_report(r);
                    }
                }
            },
            _ => {
                unreachable!();
            }
        }

        if got_results {
            samples_processed += 1;
            pb.inc();
        }
    }


    pb.finish_print("done");
}