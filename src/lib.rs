#![cfg_attr(feature = "unstable", feature(test))]

extern crate byteorder;
extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;
extern crate repr_derive;

mod arithmetics;
mod traits;
mod representation;
mod field;
mod fp;
mod weierstrass;
mod mont_inverse;
mod multiexp;
mod extension_towers;
mod pairings;
mod sliding_window_exp;
mod errors;
mod gas_meter;

pub mod public_interface;

#[cfg(test)]
mod test;

// pub use api::{API, PrecompileAPI};

#[cfg(all(feature = "unstable", test))]
mod bench;

#[cfg(test)]
mod tests {
    extern crate hex;
    extern crate rand;
    extern crate rand_xorshift;

    use num_bigint::BigUint;
    use num_traits::Num;
    use crate::field::*;
    use crate::fp::Fp;
    use crate::weierstrass::curve::*;
    use crate::traits::FieldElement;
    use crate::multiexp::{peppinger};
    use crate::weierstrass::Group;
    use crate::traits::ZeroAndOne;
    use crate::weierstrass::{CurveParameters, CurveOverFpParameters};

    const MULTIEXP_NUM_POINTS: usize = 100;

    #[test]
    fn test_multiplication_bn254() {
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let group_order = BigUint::from_str_radix("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let group_order = biguint_to_u64_vec(group_order);
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let params = CurveOverFpParameters::new(&field);

        let curve = WeierstrassCurve::new(
            group_order, 
            a_coeff, 
            b_coeff,
            &params
        );

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        // scalar is group order
        let scalar = [0x43e1f593f0000001,
                    0x2833e84879b97091,
                    0xb85045b68181585d,
                    0x30644e72e131a029];

        let res = point.mul(&scalar);

        assert!(res.is_zero());
    }

    // #[test]
    // fn test_ben_coster_bn254() {
    //     use crate::representation::ElementRepr;
    //     use rand::{RngCore, SeedableRng};
    //     use rand_xorshift::XorShiftRng;

    //     let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
    //     let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    //     let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
    //     let one = Fp::one(&field);
    //     let a_coeff = Fp::zero(&field);
    //     let mut b_coeff = one.clone();
    //     b_coeff.double();
    //     b_coeff.add_assign(&one);

    //     let curve = WeierstrassCurve::new(
    //         &group, 
    //         a_coeff, 
    //         b_coeff);

    //     let mut two = one.clone();
    //     two.double();

    //     let point = CurvePoint::point_from_xy(
    //         &curve, 
    //         one, 
    //         two);

    //     let pairs: Vec<_> = (0..MULTIEXP_NUM_POINTS).map(|_| {
    //         let mut scalar = U256Repr::default();
    //         let mut bytes = vec![0u8; 32];
    //         rng.fill_bytes(&mut bytes[1..]);
    //         scalar.read_be(& bytes[..]).unwrap();

    //         (point.clone(), scalar)
    //     }).collect();


    //     let naive_res = {
    //         let mut pairs: Vec<_> = pairs.iter().map(|el| el.0.mul(el.1)).collect();
    //         let mut acc = pairs.pop().unwrap();
    //         while let Some(p) = pairs.pop() {
    //             acc.add_assign(&p);
    //         }

    //         acc.into_xy()
    //     };

    //     let ben_coster_res = ben_coster(pairs).into_xy();

    //     assert!(ben_coster_res.0 == naive_res.0);
    //     assert!(ben_coster_res.1 == naive_res.1);
    // }

    #[test]
    fn test_peppinger_bn254() {
        use crate::representation::ElementRepr;
        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let order = BigUint::from_str_radix("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let group_order = biguint_to_u64_vec(order.clone());
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let params = CurveOverFpParameters::new(&field);

        let curve = WeierstrassCurve::new(
            group_order, 
            a_coeff, 
            b_coeff,
            &params
        );

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        let pairs: Vec<_> = (0..MULTIEXP_NUM_POINTS).map(|_| {
            let mut bytes = vec![0u8; 32];
            rng.fill_bytes(&mut bytes[..]);
            let scalar = BigUint::from_bytes_be(&bytes);
            let scalar = scalar % &order;
            let scalar = biguint_to_u64_vec(scalar);

            (point.clone(), scalar)
        }).collect();


        let naive_res = {
            let mut pairs: Vec<_> = pairs.iter().map(|el| el.0.mul(&el.1)).collect();
            let mut acc = pairs.pop().unwrap();
            while let Some(p) = pairs.pop() {
                acc.add_assign(&p);
            }

            acc.into_xy()
        };

        let ben_coster_res = peppinger(pairs).into_xy();

        assert!(ben_coster_res.0 == naive_res.0);
        assert!(ben_coster_res.1 == naive_res.1);
    }

    #[test]
    fn test_wnaf_decomposition() {
        use crate::representation::ElementRepr;
        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;
        use crate::representation::IntoWnaf;

        let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        let mut scalar = U256Repr::default();
        let mut bytes = vec![0u8; 32];
        bytes[31] = 175u8;
        // rng.fill_bytes(&mut bytes[1..]);
        scalar.read_be(& bytes[..]).unwrap();

        println!("{:#b}", 175u8);
        let wnaf = scalar.wnaf(3);

        println!("wnaf = {:?}", wnaf);
    }

    #[test]
    fn test_wnaf_mul_bn254() {
        use crate::representation::ElementRepr;
        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let group_order = BigUint::from_str_radix("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let group_order = biguint_to_u64_vec(group_order);
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let params = CurveOverFpParameters::new(&field);

        let curve = WeierstrassCurve::new(
            group_order, 
            a_coeff, 
            b_coeff,
            &params
        );

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        let mut scalar = U256Repr::default();
        let mut bytes = vec![0u8; 32];
        // bytes[31] = 2u8;
        rng.fill_bytes(&mut bytes[1..]);
        scalar.read_be(& bytes[..]).unwrap();

        let res_double_and_add  = point.clone().mul(scalar).into_xy();
        let wnaf_res = point.wnaf_mul(scalar).into_xy();

        assert!(res_double_and_add.0 == wnaf_res.0);
        assert!(res_double_and_add.1 == wnaf_res.1);
    }
}