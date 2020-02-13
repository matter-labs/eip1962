use crate::test::gas_meter::bls12;
use crate::test::gas_meter::bn;

use crate::public_interface::API;
use crate::public_interface::constants::*;
use crate::public_interface::sane_limits::*;

use crate::test::parsers::*;

use crate::public_interface::constants::*;

use rand::{Rng, thread_rng};
use rand::distributions::Distribution;
use rand::distributions::Uniform;

use num_bigint::BigUint;
use num_traits::Num;
use num_integer::Integer;

fn random_field_element<R: Rng>(modulus: &BigUint, rng: &mut R) -> BigUint {
    let mut tmp = modulus.to_bytes_be();
    rng.try_fill_bytes(&mut tmp).unwrap();
    let mut res = BigUint::from_bytes_be(&tmp);
    res = res % modulus;

    res
}

fn random_biguint<R: Rng>(bytes: usize, rng: &mut R) -> BigUint {
    let mut tmp = vec![0u8; bytes];
    rng.try_fill_bytes(&mut tmp).unwrap();
    let res = BigUint::from_bytes_be(&tmp);

    res
}

pub(crate) fn random_bls12_params<R: Rng>(limbs: usize, group_size_limbs: usize, rng: &mut R) -> JsonBls12PairingCurveParameters {
    let one = BigUint::from(1u64);

    let mut modulus_bytes = vec![0u8; limbs*8];
    let mut group_bytes = vec![0u8; group_size_limbs*8];

    rng.try_fill_bytes(&mut modulus_bytes).unwrap();
    rng.try_fill_bytes(&mut group_bytes).unwrap();

    let mut modulus = BigUint::from_bytes_be(&modulus_bytes);
    if modulus.bits() == limbs*64 {
        modulus >>= 1;
    }
    if modulus.is_even() {
        modulus -= &one;
    }

    let mut group_size = BigUint::from_bytes_be(&group_bytes);
    if group_size.is_even() {
        group_size -= &one;
    }

    let params = JsonBls12PairingCurveParameters {
        non_residue: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        is_d_type: rng.gen_bool(0.5),
        quadratic_non_residue_0: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        quadratic_non_residue_1: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        x: (BigUint::from(0u64), false),
        q: modulus.clone(),
        r: group_size,
        a: BigUint::from(0u64),
        b: random_field_element(&modulus, rng),
        a_twist_0: random_field_element(&modulus, rng),
        a_twist_1: random_field_element(&modulus, rng),
        b_twist_0: random_field_element(&modulus, rng),
        b_twist_1: random_field_element(&modulus, rng),
        g1_x: random_field_element(&modulus, rng),
        g1_y: random_field_element(&modulus, rng),
        g2_x_0: random_field_element(&modulus, rng),
        g2_x_1: random_field_element(&modulus, rng),
        g2_y_0: random_field_element(&modulus, rng),
        g2_y_1: random_field_element(&modulus, rng),
        g1_mul_vectors: vec![],
        g2_mul_vectors: vec![],
    };

    params
}


pub(crate) fn random_bn_params<R: Rng>(limbs: usize, group_size_limbs: usize, rng: &mut R) -> JsonBnPairingCurveParameters {
    let one = BigUint::from(1u64);

    let mut modulus_bytes = vec![0u8; limbs*8];
    let mut group_bytes = vec![0u8; group_size_limbs*8];

    rng.try_fill_bytes(&mut modulus_bytes).unwrap();
    rng.try_fill_bytes(&mut group_bytes).unwrap();

    let mut modulus = BigUint::from_bytes_be(&modulus_bytes);
    if modulus.bits() == limbs*64 {
        modulus >>= 1;
    }
    if modulus.is_even() {
        modulus -= &one;
    }

    let mut group_size = BigUint::from_bytes_be(&group_bytes);
    if group_size.is_even() {
        group_size -= &one;
    }

    let params = JsonBnPairingCurveParameters {
        non_residue: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        is_d_type: rng.gen_bool(0.5),
        quadratic_non_residue_0: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        quadratic_non_residue_1: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        x: (BigUint::from(0u64), false),
        q: modulus.clone(),
        r: group_size,
        a: (BigUint::from(0u64), false),
        b: (random_field_element(&modulus, rng), false),
        a_twist_0: random_field_element(&modulus, rng),
        a_twist_1: random_field_element(&modulus, rng),
        b_twist_0: random_field_element(&modulus, rng),
        b_twist_1: random_field_element(&modulus, rng),
        g1_x: random_field_element(&modulus, rng),
        g1_y: random_field_element(&modulus, rng),
        g2_x_0: random_field_element(&modulus, rng),
        g2_x_1: random_field_element(&modulus, rng),
        g2_y_0: random_field_element(&modulus, rng),
        g2_y_1: random_field_element(&modulus, rng),
        g1_mul_vectors: vec![],
        g2_mul_vectors: vec![],
    };

    params
}

pub(crate) fn random_mnt4_params<R: Rng>(limbs: usize, group_size_limbs: usize, rng: &mut R) -> JsonMnt4PairingCurveParameters {
    let one = BigUint::from(1u64);

    let mut modulus_bytes = vec![0u8; limbs*8];
    let mut group_bytes = vec![0u8; group_size_limbs*8];

    rng.try_fill_bytes(&mut modulus_bytes).unwrap();
    rng.try_fill_bytes(&mut group_bytes).unwrap();

    let mut modulus = BigUint::from_bytes_be(&modulus_bytes);
    if modulus.bits() == limbs*64 {
        modulus >>= 1;
    }
    if modulus.is_even() {
        modulus -= &one;
    }

    let mut group_size = BigUint::from_bytes_be(&group_bytes);
    if group_size.is_even() {
        group_size -= &one;
    }

    let params = JsonMnt4PairingCurveParameters {
        non_residue: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        x: (BigUint::from(0u64), false),
        exp_w0: (BigUint::from(0u64), false),
        exp_w1: BigUint::from(0u64),
        q: modulus.clone(),
        r: group_size,
        a: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        b: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        a_twist_0: random_field_element(&modulus, rng),
        a_twist_1: random_field_element(&modulus, rng),
        b_twist_0: random_field_element(&modulus, rng),
        b_twist_1: random_field_element(&modulus, rng),
        g1_x: random_field_element(&modulus, rng),
        g1_y: random_field_element(&modulus, rng),
        g2_x_0: random_field_element(&modulus, rng),
        g2_x_1: random_field_element(&modulus, rng),
        g2_y_0: random_field_element(&modulus, rng),
        g2_y_1: random_field_element(&modulus, rng),
        g1_mul_vectors: vec![],
        g2_mul_vectors: vec![],
    };

    params
}


pub(crate) fn random_mnt6_params<R: Rng>(limbs: usize, group_size_limbs: usize, rng: &mut R) -> JsonMnt6PairingCurveParameters {
    let one = BigUint::from(1u64);

    let mut modulus_bytes = vec![0u8; limbs*8];
    let mut group_bytes = vec![0u8; group_size_limbs*8];

    rng.try_fill_bytes(&mut modulus_bytes).unwrap();
    rng.try_fill_bytes(&mut group_bytes).unwrap();

    let mut modulus = BigUint::from_bytes_be(&modulus_bytes);
    if modulus.bits() == limbs*64 {
        modulus >>= 1;
    }
    if modulus.is_even() {
        modulus -= &one;
    }

    let mut group_size = BigUint::from_bytes_be(&group_bytes);
    if group_size.is_even() {
        group_size -= &one;
    }

    let params = JsonMnt6PairingCurveParameters {
        non_residue: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        x: (BigUint::from(0u64), false),
        exp_w0: (BigUint::from(0u64), false),
        exp_w1: BigUint::from(0u64),
        q: modulus.clone(),
        r: group_size,
        a: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        b: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        a_twist_0: random_field_element(&modulus, rng),
        a_twist_1: random_field_element(&modulus, rng),
        a_twist_2: random_field_element(&modulus, rng),
        b_twist_0: random_field_element(&modulus, rng),
        b_twist_1: random_field_element(&modulus, rng),
        b_twist_2: random_field_element(&modulus, rng),
        g1_x: random_field_element(&modulus, rng),
        g1_y: random_field_element(&modulus, rng),
        g2_x_0: random_field_element(&modulus, rng),
        g2_x_1: random_field_element(&modulus, rng),
        g2_x_2: random_field_element(&modulus, rng),
        g2_y_0: random_field_element(&modulus, rng),
        g2_y_1: random_field_element(&modulus, rng),
        g2_y_2: random_field_element(&modulus, rng),
        g1_mul_vectors: vec![],
        g2_mul_vectors: vec![],
    };

    params
}

pub(crate) fn random_mul_params_a_non_zero_ext3<R: Rng>(
    limbs: usize, 
    group_size_limbs: usize, 
    num_mul_points: usize, 
    rng: &mut R
) -> (JsonMnt6PairingCurveParameters, JsonG1PointScalarMultiplicationPair, JsonG2Ext3PointScalarMultiplicationPair) {
    let one = BigUint::from(1u64);

    let mut modulus_bytes = vec![0u8; limbs*8];

    rng.try_fill_bytes(&mut modulus_bytes).unwrap();

    let mut modulus = BigUint::from_bytes_be(&modulus_bytes);
    if modulus.bits() == limbs*64 {
        modulus >>= 1;
    }
    if modulus.is_even() {
        modulus -= &one;
    }

    // group size is `group_size_limbs` of u64::MAX
    let mut group_size = BigUint::from(1u64) << (group_size_limbs * 64);
    group_size -= &one;

    // it's the worst case
    let mul_scalar = group_size.clone();

    let g1_worst_case_mul_pair = JsonG1PointScalarMultiplicationPair {
        scalar: mul_scalar.clone(),
        base_x: random_field_element(&modulus, rng),
        base_y: random_field_element(&modulus, rng),
        result_x: BigUint::from(0u64),
        result_y: BigUint::from(0u64),
    };

    let g2_worst_case_mul_pair = JsonG2Ext3PointScalarMultiplicationPair {
        scalar: mul_scalar,
        base_x_0: random_field_element(&modulus, rng),
        base_x_1: random_field_element(&modulus, rng),
        base_x_2: random_field_element(&modulus, rng),
        base_y_0: random_field_element(&modulus, rng),
        base_y_1: random_field_element(&modulus, rng),
        base_y_2: random_field_element(&modulus, rng),
        result_x_0: BigUint::from(0u64),
        result_x_1: BigUint::from(0u64),
        result_x_2: BigUint::from(0u64),
        result_y_0: BigUint::from(0u64),
        result_y_1: BigUint::from(0u64),
        result_y_2: BigUint::from(0u64),
    };

    let mut g1_multiexp_data = vec![];
    let mut g2_multiexp_data = vec![];

    for _ in 0..num_mul_points {
        let g1 = JsonG1PointScalarMultiplicationPair {
            scalar: random_field_element(&group_size, rng),
            base_x: random_field_element(&modulus, rng),
            base_y: random_field_element(&modulus, rng),
            result_x: BigUint::from(0u64),
            result_y: BigUint::from(0u64),
        };

        let g2 = JsonG2Ext3PointScalarMultiplicationPair {
            scalar: random_field_element(&group_size, rng),
            base_x_0: random_field_element(&modulus, rng),
            base_x_1: random_field_element(&modulus, rng),
            base_x_2: random_field_element(&modulus, rng),
            base_y_0: random_field_element(&modulus, rng),
            base_y_1: random_field_element(&modulus, rng),
            base_y_2: random_field_element(&modulus, rng),
            result_x_0: BigUint::from(0u64),
            result_x_1: BigUint::from(0u64),
            result_x_2: BigUint::from(0u64),
            result_y_0: BigUint::from(0u64),
            result_y_1: BigUint::from(0u64),
            result_y_2: BigUint::from(0u64),
        };

        g1_multiexp_data.push(g1);
        g2_multiexp_data.push(g2);
    }

    let params = JsonMnt6PairingCurveParameters {
        non_residue: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        x: (BigUint::from(0u64), false),
        exp_w0: (BigUint::from(0u64), false),
        exp_w1: BigUint::from(0u64),
        q: modulus.clone(),
        r: group_size,
        a: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        b: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        a_twist_0: random_field_element(&modulus, rng),
        a_twist_1: random_field_element(&modulus, rng),
        a_twist_2: random_field_element(&modulus, rng),
        b_twist_0: random_field_element(&modulus, rng),
        b_twist_1: random_field_element(&modulus, rng),
        b_twist_2: random_field_element(&modulus, rng),
        g1_x: random_field_element(&modulus, rng),
        g1_y: random_field_element(&modulus, rng),
        g2_x_0: random_field_element(&modulus, rng),
        g2_x_1: random_field_element(&modulus, rng),
        g2_x_2: random_field_element(&modulus, rng),
        g2_y_0: random_field_element(&modulus, rng),
        g2_y_1: random_field_element(&modulus, rng),
        g2_y_2: random_field_element(&modulus, rng),
        g1_mul_vectors: g1_multiexp_data,
        g2_mul_vectors: g2_multiexp_data,
    };

    (params, g1_worst_case_mul_pair, g2_worst_case_mul_pair)
}


pub(crate) fn random_mul_params_a_is_zero_ext3<R: Rng>(
    limbs: usize, 
    group_size_limbs: usize, 
    num_mul_points: usize, 
    rng: &mut R
) -> (JsonMnt6PairingCurveParameters, JsonG1PointScalarMultiplicationPair, JsonG2Ext3PointScalarMultiplicationPair) {
    let one = BigUint::from(1u64);

    let mut modulus_bytes = vec![0u8; limbs*8];

    rng.try_fill_bytes(&mut modulus_bytes).unwrap();

    let mut modulus = BigUint::from_bytes_be(&modulus_bytes);
    if modulus.bits() == limbs*64 {
        modulus >>= 1;
    }
    if modulus.is_even() {
        modulus -= &one;
    }

    let mut group_size = BigUint::from(1u64) << (group_size_limbs * 64);
    group_size -= &one;

    // it's the worst case
    let mul_scalar = group_size.clone();

    let g1_worst_case_mul_pair = JsonG1PointScalarMultiplicationPair {
        scalar: mul_scalar.clone(),
        base_x: random_field_element(&modulus, rng),
        base_y: random_field_element(&modulus, rng),
        result_x: BigUint::from(0u64),
        result_y: BigUint::from(0u64),
    };

    let g2_worst_case_mul_pair = JsonG2Ext3PointScalarMultiplicationPair {
        scalar: mul_scalar,
        base_x_0: random_field_element(&modulus, rng),
        base_x_1: random_field_element(&modulus, rng),
        base_x_2: random_field_element(&modulus, rng),
        base_y_0: random_field_element(&modulus, rng),
        base_y_1: random_field_element(&modulus, rng),
        base_y_2: random_field_element(&modulus, rng),
        result_x_0: BigUint::from(0u64),
        result_x_1: BigUint::from(0u64),
        result_x_2: BigUint::from(0u64),
        result_y_0: BigUint::from(0u64),
        result_y_1: BigUint::from(0u64),
        result_y_2: BigUint::from(0u64),
    };

    let mut g1_multiexp_data = vec![];
    let mut g2_multiexp_data = vec![];

    for _ in 0..num_mul_points {
        let g1 = JsonG1PointScalarMultiplicationPair {
            scalar: random_field_element(&group_size, rng),
            base_x: random_field_element(&modulus, rng),
            base_y: random_field_element(&modulus, rng),
            result_x: BigUint::from(0u64),
            result_y: BigUint::from(0u64),
        };

        let g2 = JsonG2Ext3PointScalarMultiplicationPair {
            scalar: random_field_element(&group_size, rng),
            base_x_0: random_field_element(&modulus, rng),
            base_x_1: random_field_element(&modulus, rng),
            base_x_2: random_field_element(&modulus, rng),
            base_y_0: random_field_element(&modulus, rng),
            base_y_1: random_field_element(&modulus, rng),
            base_y_2: random_field_element(&modulus, rng),
            result_x_0: BigUint::from(0u64),
            result_x_1: BigUint::from(0u64),
            result_x_2: BigUint::from(0u64),
            result_y_0: BigUint::from(0u64),
            result_y_1: BigUint::from(0u64),
            result_y_2: BigUint::from(0u64),
        };

        g1_multiexp_data.push(g1);
        g2_multiexp_data.push(g2);
    }

    let params = JsonMnt6PairingCurveParameters {
        non_residue: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        x: (BigUint::from(0u64), false),
        exp_w0: (BigUint::from(0u64), false),
        exp_w1: BigUint::from(0u64),
        q: modulus.clone(),
        r: group_size,
        a: (BigUint::from(0u64), false),
        b: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        a_twist_0: BigUint::from(0u64),
        a_twist_1: BigUint::from(0u64),
        a_twist_2: BigUint::from(0u64),
        b_twist_0: random_field_element(&modulus, rng),
        b_twist_1: random_field_element(&modulus, rng),
        b_twist_2: random_field_element(&modulus, rng),
        g1_x: random_field_element(&modulus, rng),
        g1_y: random_field_element(&modulus, rng),
        g2_x_0: random_field_element(&modulus, rng),
        g2_x_1: random_field_element(&modulus, rng),
        g2_x_2: random_field_element(&modulus, rng),
        g2_y_0: random_field_element(&modulus, rng),
        g2_y_1: random_field_element(&modulus, rng),
        g2_y_2: random_field_element(&modulus, rng),
        g1_mul_vectors: g1_multiexp_data,
        g2_mul_vectors: g2_multiexp_data,
    };

    (params, g1_worst_case_mul_pair, g2_worst_case_mul_pair)
}


pub(crate) fn random_mul_params_a_non_zero_ext2<R: Rng>(
    limbs: usize, 
    group_size_limbs: usize, 
    num_mul_points: usize, 
    rng: &mut R
) -> (JsonMnt4PairingCurveParameters, JsonG1PointScalarMultiplicationPair, JsonG2PointScalarMultiplicationPair) {
    let one = BigUint::from(1u64);

    let mut modulus_bytes = vec![0u8; limbs*8];

    rng.try_fill_bytes(&mut modulus_bytes).unwrap();

    let mut modulus = BigUint::from_bytes_be(&modulus_bytes);
    if modulus.bits() == limbs*64 {
        modulus >>= 1;
    }
    if modulus.is_even() {
        modulus -= &one;
    }

    let mut group_size = BigUint::from(1u64) << (group_size_limbs * 64);
    group_size -= &one;

    // it's the worst case
    let mul_scalar = group_size.clone();

    let g1_worst_case_mul_pair = JsonG1PointScalarMultiplicationPair {
        scalar: mul_scalar.clone(),
        base_x: random_field_element(&modulus, rng),
        base_y: random_field_element(&modulus, rng),
        result_x: BigUint::from(0u64),
        result_y: BigUint::from(0u64),
    };

    let g2_worst_case_mul_pair = JsonG2PointScalarMultiplicationPair {
        scalar: mul_scalar,
        base_x_0: random_field_element(&modulus, rng),
        base_x_1: random_field_element(&modulus, rng),
        base_y_0: random_field_element(&modulus, rng),
        base_y_1: random_field_element(&modulus, rng),
        result_x_0: BigUint::from(0u64),
        result_x_1: BigUint::from(0u64),
        result_y_0: BigUint::from(0u64),
        result_y_1: BigUint::from(0u64),
    };

    let mut g1_multiexp_data = vec![];
    let mut g2_multiexp_data = vec![];

    for _ in 0..num_mul_points {
        let g1 = JsonG1PointScalarMultiplicationPair {
            scalar: random_field_element(&group_size, rng),
            base_x: random_field_element(&modulus, rng),
            base_y: random_field_element(&modulus, rng),
            result_x: BigUint::from(0u64),
            result_y: BigUint::from(0u64),
        };

        let g2 = JsonG2PointScalarMultiplicationPair {
            scalar: random_field_element(&group_size, rng),
            base_x_0: random_field_element(&modulus, rng),
            base_x_1: random_field_element(&modulus, rng),
            base_y_0: random_field_element(&modulus, rng),
            base_y_1: random_field_element(&modulus, rng),
            result_x_0: BigUint::from(0u64),
            result_x_1: BigUint::from(0u64),
            result_y_0: BigUint::from(0u64),
            result_y_1: BigUint::from(0u64),
        };

        g1_multiexp_data.push(g1);
        g2_multiexp_data.push(g2);
    }

    let params = JsonMnt4PairingCurveParameters {
        non_residue: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        x: (BigUint::from(0u64), false),
        exp_w0: (BigUint::from(0u64), false),
        exp_w1: BigUint::from(0u64),
        q: modulus.clone(),
        r: group_size,
        a: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        b: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        a_twist_0: random_field_element(&modulus, rng),
        a_twist_1: random_field_element(&modulus, rng),
        b_twist_0: random_field_element(&modulus, rng),
        b_twist_1: random_field_element(&modulus, rng),
        g1_x: random_field_element(&modulus, rng),
        g1_y: random_field_element(&modulus, rng),
        g2_x_0: random_field_element(&modulus, rng),
        g2_x_1: random_field_element(&modulus, rng),
        g2_y_0: random_field_element(&modulus, rng),
        g2_y_1: random_field_element(&modulus, rng),
        g1_mul_vectors: g1_multiexp_data,
        g2_mul_vectors: g2_multiexp_data,
    };

    (params, g1_worst_case_mul_pair, g2_worst_case_mul_pair)
}


pub(crate) fn random_mul_params_a_is_zero_ext2<R: Rng>(
    limbs: usize, 
    group_size_limbs: usize, 
    num_mul_points: usize, 
    rng: &mut R
) -> (JsonMnt4PairingCurveParameters, JsonG1PointScalarMultiplicationPair, JsonG2PointScalarMultiplicationPair) {
    let one = BigUint::from(1u64);

    let mut modulus_bytes = vec![0u8; limbs*8];

    rng.try_fill_bytes(&mut modulus_bytes).unwrap();

    let mut modulus = BigUint::from_bytes_be(&modulus_bytes);
    if modulus.bits() == limbs*64 {
        modulus >>= 1;
    }
    if modulus.is_even() {
        modulus -= &one;
    }

    let mut group_size = BigUint::from(1u64) << (group_size_limbs * 64);
    group_size -= &one;

    // it's almost the worst case
    let mul_scalar = group_size.clone();

    let g1_worst_case_mul_pair = JsonG1PointScalarMultiplicationPair {
        scalar: mul_scalar.clone(),
        base_x: random_field_element(&modulus, rng),
        base_y: random_field_element(&modulus, rng),
        result_x: BigUint::from(0u64),
        result_y: BigUint::from(0u64),
    };

    let g2_worst_case_mul_pair = JsonG2PointScalarMultiplicationPair {
        scalar: mul_scalar,
        base_x_0: random_field_element(&modulus, rng),
        base_x_1: random_field_element(&modulus, rng),
        base_y_0: random_field_element(&modulus, rng),
        base_y_1: random_field_element(&modulus, rng),
        result_x_0: BigUint::from(0u64),
        result_x_1: BigUint::from(0u64),
        result_y_0: BigUint::from(0u64),
        result_y_1: BigUint::from(0u64),
    };

    let mut g1_multiexp_data = vec![];
    let mut g2_multiexp_data = vec![];

    for _ in 0..num_mul_points {
        let g1 = JsonG1PointScalarMultiplicationPair {
            scalar: random_field_element(&group_size, rng),
            base_x: random_field_element(&modulus, rng),
            base_y: random_field_element(&modulus, rng),
            result_x: BigUint::from(0u64),
            result_y: BigUint::from(0u64),
        };

        let g2 = JsonG2PointScalarMultiplicationPair {
            scalar: random_field_element(&group_size, rng),
            base_x_0: random_field_element(&modulus, rng),
            base_x_1: random_field_element(&modulus, rng),
            base_y_0: random_field_element(&modulus, rng),
            base_y_1: random_field_element(&modulus, rng),
            result_x_0: BigUint::from(0u64),
            result_x_1: BigUint::from(0u64),
            result_y_0: BigUint::from(0u64),
            result_y_1: BigUint::from(0u64),
        };

        g1_multiexp_data.push(g1);
        g2_multiexp_data.push(g2);
    }

    let params = JsonMnt4PairingCurveParameters {
        non_residue: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        x: (BigUint::from(0u64), false),
        exp_w0: (BigUint::from(0u64), false),
        exp_w1: BigUint::from(0u64),
        q: modulus.clone(),
        r: group_size,
        a: (BigUint::from(0u64), false),
        b: (random_field_element(&modulus, rng), rng.gen_bool(0.5)),
        a_twist_0: BigUint::from(0u64),
        a_twist_1: BigUint::from(0u64),
        b_twist_0: random_field_element(&modulus, rng),
        b_twist_1: random_field_element(&modulus, rng),
        g1_x: random_field_element(&modulus, rng),
        g1_y: random_field_element(&modulus, rng),
        g2_x_0: random_field_element(&modulus, rng),
        g2_x_1: random_field_element(&modulus, rng),
        g2_y_0: random_field_element(&modulus, rng),
        g2_y_1: random_field_element(&modulus, rng),
        g1_mul_vectors: g1_multiexp_data,
        g2_mul_vectors: g2_multiexp_data,
    };

    (params, g1_worst_case_mul_pair, g2_worst_case_mul_pair)
}