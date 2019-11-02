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
    assert!(std::option_env!("GAS_METERING").is_some());

    use crate::public_interface::decode_utils::*;
    use crate::test::g1_ops::mnt6 as g1_mnt6;

    use crate::test::g2_ops::mnt6 as g2_mnt6;
    use crate::test::gas_meter::arithmetic_ops::*;

    use rand::{SeedableRng};
    use rand_xorshift::XorShiftRng;

    let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    // low number of samples is ok cause we test many things redundantly
    const SAMPLES: u64 = 10;

    let mut multiexp_len = vec![2, 4, 8, 16, 32, 64, 128];
    multiexp_len.reverse();

    for num_limbs in 4..=16 {
        for num_group_limbs in 1..=16 {
    // for num_limbs in 16..=16 {
        // for num_group_limbs in 16..=16 {
            let (mut curve_ext3, g1_worst_case_pair, g2_worst_case_pair) = gen_params::random_mul_params_a_non_zero_ext3(num_limbs, num_group_limbs, 128, &mut rng);
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
                // println!("{}", hex::encode(&input_data));

                let mut total = 0;

                for _ in 0..SAMPLES {
                    let now = Instant::now();
                    run(&input_data);
                    let elapsed = now.elapsed();

                    total += elapsed.as_micros();
                }

                (total as u64) / SAMPLES
            };

            println!("G2 Ext3: For {} limbs and {} group units curve construction taken {} microseconds", num_limbs, num_group_limbs, addition_timing_g2);

        }
    }
}


#[test]
// #[ignore]
fn flame_one_off_field_and_curve_constuction_g2_ext3() {
    assert!(std::option_env!("GAS_METERING").is_some());
    let input_data = hex::decode("048076038b8e27878715b7808bb806020707b3254155972e56f81693bc9c25351423d9e4a07f3cdb4b3244da93caaef3acc8ca54d634d37fac5c6ccf3d88427823276888552176b9381ee65718f8df5145447303da4d23f4f6926c326b2c17514a42a68712df9fb7e4fa6d12e90fee38bb9f8ca3dbd453badc0bbfbf654c8582db8d030b625f820d60c8ebdcfe7ba4a4400a02d4353c7bb6fa5cd31f73aaa8ce09f5f89111fd6232de40566e6d89ad68ba62fdcc4d997743bc7827df477cc162142529a7cb6f94cde717766f147655a79cf1546d5729beeae6c9771051f9924cf3354a585753bece089586ddb8bcded59ff63b64e63b1f354ec6318062870a325c008516a1f3a89f336dc61ef52bcf3198d6f271ee03dbaa8ec2579ad47fd79e1312419fbf42fcd50f47cc4d91e7a41f5d9f712b70ff23c8eed1d570bc96fc66002ab8716bec086cbc51cc462a0ad2f8d6da3b6bf9040c47c7626721363a02528b30660e76d92746ea863fb931a62f3c15aca261762c4c4835373b6b6392b5baa306740c8946451114048da4a6210a2541d25b5f78353c4297e8df81aab0380d6e5fd0b5c72f42ca22a6d4df0b97724b0351c79820f13296f06492c93c7da27aa7ba6628d8148326c16e245b2a5524a9db661d89a883fcc78dffdf7d6cacd462502883b3722e9c2b61a43ff18182e7078b23a156ec7976d5ddde467c4d8c5d2b4565a75b718cc9cd4273a24d7b0040d258e6d3d5e58dfa56b0bd28f9c023ea69ad5613ce592581ad865f5597ffc238fb429b9a4c177c1db472772a3888d8aac94da5ea8fa0e5af3b94dc2a5f8d2dcffc30b352c20bb24964ad0feb34d332565968ee8f1f6b2847a84e7c749ba3d1f6072e8872356e36159a5f498315379109aac5e3d767c1590548590da3ab5bac420665ccdc01d55517ba095bded9600afe4bce5141bde615f8dc3605fe6d8da453c432489cfe0f3fcb93735b1c905ae36b3a9be0b6e73293852d0e934235ebd5333e9f5e81f881469d1b167e4782d16edbbe02bcaed01f5c07ee1f25b0a92a6a1d48728d537645c8b46ee68a2c34bef960519fde58537f4020c13082bfac3420c02a861dcd045fc9d0cf8906c64b6f93070e45c62fe886e5a9d633e847df6d28d713c68a5892c91f06d446b07d62705d12b1efc3df18051673fd4f058bd4e6f9e8caecb297c6a85ec81e830ec5934126fbb2afa30023ad5abb0dd0f1eb95ff8642e66af2d8e25fbab6fca67bcdb3d01da74e9f647107cb89c5076fb6a79595a840d6368647d21cdccc199c06d8b1be06dc611340aa77b05686f0125c6a3cda46832e0a372977d2d11a4f8890905b37b225b92f45b42c11932081375534bfd6e980fc194c1d0871eea50bb39a3889c608f64da2f68b0ae4465d4206a9b050a22500d38cff3eb259f35478de90ac98de8a816e43d14380ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff4280362d15410b961183c5561d1d6ff1e2a7c26806400da49d82aadea0eba45701ab242a33d3e7e4581bd36ec9045a1ea0d82911fdef8dcb0f24723d2413c0f3243973ac50d9250807d6dba48ba871575801fe620ebb37446136263edef4c9e4deddaf76cef7b98f664d2380b836e72b375fae230bdbe9636de23bf904597beb45f0f22c54c9da81c060e3f9213be678fab41f803569c050ca1900e9f0883c4c09c78e3d15b6c015db132658b9513c00384b91cb8bbc38b85d2f86587ab560441770590cbf4873f456294df3e965b559f5c4f479f2cf22225da33c404c1155e4d0dc89b72685f1e99ba2a4a6d2315ef2a2810d12c552a9a6f50c161cac5d9a9b4c70198e53496e140106ffa3270d173cf981c979427b8b1fbcb91f442dfae5c0816d16f50a4903b1cfd468087d73b3e76f708fa14d366f0a6355ef4418eebb3868ed6de23672c3912069b9afe1999548363b1ec561b3c64fd023522bcaedeccf8ab81ffb3341388167bece3b6986e49377a739953a5d4bd612fc3898b2f6582f31ee027766cf2cfb7d64aba43cfbae4cfee7bd2cf04a77b12f58004ec83d4c7588a60db6735a26bdd8acabf1586d1de595aeff3efb9e85b17f29e41d36f09d13d67fd5d071513a4d7826a9713b74498d8b8119f6e0f7376958869551fdc8bb97581694adfd7cca7b9a929086955e25cfffb70ec2d7160bea4e4d0de82c77390d3b883fba058c300870d7d58a1151934ce876ab0a1c837786f5aa3da27a7cf0743a58bd25f250e81ff83bcd50378e29a7606e9a6841cca835f723830c664c30e78ed7625c87d2db2f10d0523a40f8645f858af015562d4ee49fb3304821e4ee301be4dba1eb496a34b2a80bf95e9483bea94406c0536845a6110aaf027e4f2ca0687171c3fb01f0c9a86f66ed1a11605fe39b4489664e81c4f860a5f3abd3afac55ff11a11371a82bbf2715fbd12937461dd9228477e801216f19456af57289b59e82ee8464fe8a9461fec95ef29f60a8dec12903c1e7e363f76d69879980eb2ee4313f5fd1bd6f0d12ca8300a1f4582f83a65390e58998b58e3f206afab7d1540c46693ae9fad9643cec55fda77c04b70463977f80ba7eb01f24307a2e28e2cfd109da63655c9db2f08633d3ccb6ee16d1abcab579efa4f20d98f562d356ceb9ef0e69e9f40219af7c143fbdfb497195229e71e1b9d18edb088187ece349b433343fda5f54a43882c31607e8b01b9cd9ecf2f754f7dae2b244dc117b4b5855ef71eae4c9cb2aba83032f67f084857a02d85023ece5d39dbe02ec4b77ca11cba1bb30985756c3ca075114c92f231575d6befafe4084517f1166a47376867bd108285b4dbf34ce55544c00eebf7c51603cc3d921e9bd8eff13ee00dcdcb84109e2fb730105809f64ea522983d6bbb62f7e2e8cbf702685e9be10e2ef71f8187672643424fe335e160cc3c2e702bfab94c0dd20cc1fd10295a1fa3d589dbdb4f4f8b6a9408625b0ca8fcbfb21d34eec2d8e24e9a30d2d3b32d7a37d110b13afbfeaa99147f8df36eabec4b8963e324f46296a5948edd2d6fb02d11e208918ef04923b77283d0a7bb9e17a27e66851792fdd605cc0a339028b8985390fd024374c765b434349b384fdf0e85e63f74b2706b09be30ec8530e463b09739b4675c8ad1003b4ae2f55bfa341e423466396051e49e6b6bb9367a2ff6c4ee623a8c8ea0255549a541b72ac4a8ea97a5998b17acb81b29eff66763fd37a992fc7de8809bd07c9899e54f8e49e6e5ae548e9ef66dee98188386114562b358f6f20ea8fb8e3dd101e5ffee2dcaa870e6e5881a38aac1278b90514d63a99c57514f387ba6f7cf194c68bc8d91ac8c489ee87dbfc4b94c93c8bbd5fc04c27db8b02303f3a659054ffe0c7e78bf16706968ffcb92694926e13418ab6e015588fa2bb85d270364c28b3682accc3939283b870357cf83683350baf73aa0d3d68bda82a0f6ae7e51746630d9fe37d42d078156a92e62f1bd8601b4811f8b35e0b858ab48e788394ce672e1369df13541805188c8faea95a651b5303fcb88574396789ae6118d10c70a7e3b7a98ec1ac353da87b3515692b69197b2910cfcad9a1994580cf278b780e3b148f7e595988aedb54b3c943ddc761492e53ac848725770b7fe8ad864c45254d00");
    let input_data = input_data.unwrap();
    run(&input_data);
}


fn run(bytes: &[u8]) {
    let _ = API::run(&bytes);
}