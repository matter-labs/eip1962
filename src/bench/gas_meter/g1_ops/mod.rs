extern crate test as rust_test;
use self::rust_test::{Bencher, black_box};

#[bench]
fn bench_metering_from_vectors_g1_mul(b: &mut Bencher) {
    use crate::test::parsers::*;
    use crate::test::g1_ops::bls12::*;
    use crate::public_interface::API;
    use crate::public_interface::constants::*;
    use crate::public_interface::decode_utils::*;
    use crate::gas_meter::GasMeter;
    const GAS_PER_MICROSECOND: f64 = 15f64;
    const SAMPLES: u64 = 1000u64;
    use std::time::Instant;
    let curves = read_dir_and_grab_curves("src/test/test_vectors/bls12/");
    assert!(curves.len() != 0);
    for (curve, _) in curves.into_iter() {
        let (calldata, modulus_len, group_len) = assemble_single_curve_params(curve.clone());
        let modulus = curve.q.clone();
        let order = curve.r.clone();
        let num_limbs = num_limbs_for_modulus(&modulus).expect("must work");
        let num_words_for_order = num_units_for_group_order(&order).expect("must work");
        for pair in curve.g1_mul_vectors.into_iter() {
            // let mut subbencher = b.clone();
            let (points_data, _expected_result) = assemble_single_point_scalar_pair(pair, modulus_len, group_len);

            let mut input_data = vec![OPERATION_G1_MUL];
            input_data.extend(calldata.clone());
            input_data.extend(points_data);

            let input_data = black_box(input_data);

            let proposed_gas = GasMeter::meter(&input_data[..]).expect("Must meter some gas") as f64;
            let now = Instant::now();
            for _ in 0..SAMPLES {
                API::run(&input_data[..]).expect("api call must work");
            }
            let elapsed = now.elapsed().as_micros();
            // b.iter(|| {
            //     API::run(&input_data[..]).unwrap();
            // });
            let ms_per_iter = (elapsed  as f64) / (SAMPLES as f64);

            let gas_per_iter = GAS_PER_MICROSECOND * ms_per_iter;
            println!("Modulus limbs = {}, order words = {}", num_limbs, num_words_for_order);
            println!("Gas from meter = {}", proposed_gas);
            println!("Gas from timing = {}", gas_per_iter);
        }
    }
}

#[bench]
fn bench_metering_from_vectors_g1_add(b: &mut Bencher) {
    use crate::test::parsers::*;
    use crate::test::g1_ops::bls12::*;
    use crate::public_interface::API;
    use crate::public_interface::constants::*;
    use crate::public_interface::decode_utils::*;
    use crate::gas_meter::GasMeter;
    const GAS_PER_MICROSECOND: f64 = 15f64;
    const SAMPLES: u64 = 1000u64;
    use std::time::Instant;
    let curves = read_dir_and_grab_curves("src/test/test_vectors/bls12/");
    assert!(curves.len() != 0);
    for (curve, _) in curves.into_iter() {
        let (calldata, modulus_len, group_len) = assemble_single_curve_params(curve.clone());
        let modulus = curve.q.clone();
        let order = curve.r.clone();
        let num_limbs = num_limbs_for_modulus(&modulus).expect("must work");
        let num_words_for_order = num_units_for_group_order(&order).expect("must work");
        for pair in curve.g1_mul_vectors.into_iter() {
            // let mut subbencher = b.clone();
            let (points_data, _expected_result) = assemble_single_points_addition_pair(pair, modulus_len, group_len);

            let mut input_data = vec![OPERATION_G1_ADD];
            input_data.extend(calldata.clone());
            input_data.extend(points_data);

            let input_data = black_box(input_data);

            let proposed_gas = GasMeter::meter(&input_data[..]).expect("Must meter some gas") as f64;
            let now = Instant::now();
            for _ in 0..SAMPLES {
                API::run(&input_data[..]).expect("api call must work");
            }
            let elapsed = now.elapsed().as_micros();
            // b.iter(|| {
            //     API::run(&input_data[..]).unwrap();
            // });
            let ms_per_iter = (elapsed  as f64) / (SAMPLES as f64);

            let gas_per_iter = GAS_PER_MICROSECOND * ms_per_iter;
            println!("Modulus limbs = {}, order words = {}", num_limbs, num_words_for_order);
            println!("Gas from meter = {}", proposed_gas);
            println!("Gas from timing = {}", gas_per_iter);
        }
    }
}
