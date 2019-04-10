#![feature(test)]

extern crate byteorder;
extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;
extern crate hex;
extern crate rand;
extern crate rand_xorshift;
extern crate repr_derive;

mod arithmetics;
mod traits;
mod representation;
mod field;
mod fp;
mod weierstrass;
mod mont_inverse;
mod multiexp;
mod api;
mod extension_towers;


#[cfg(test)]
mod test;

pub use api::{API, PrecompileAPI};

extern crate test as rust_test;

#[cfg(test)]
mod tests {
    use crate::field::*;
    use crate::fp::Fp;
    use crate::weierstrass::curve::*;
    use crate::traits::FieldElement;
    use rust_test::Bencher;
    use crate::multiexp::ben_coster;

    const MULTIEXP_NUM_POINTS: usize = 100;

    #[bench]
    fn bench_doubling_bn254(b: &mut Bencher) {
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff);

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        b.iter(|| point.clone().double());
    }

    #[bench]
    fn bench_addition_bn254(b: &mut Bencher) {
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff);

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        let mut another = point.clone();
        another.double();

        b.iter(|| point.clone().add_assign(&another));
    }

    #[bench]
    fn bench_multiplication_bn254(b: &mut Bencher) {
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff);

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        // scalar is order - 1
        let scalar = [0x43e1f593f0000000,
                    0x2833e84879b97091,
                    0xb85045b68181585d,
                    0x30644e72e131a029];
        
        b.iter(|| point.mul(&scalar));
    }

    #[bench]
    fn bench_multiplication_bn254_into_affine(b: &mut Bencher) {
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff);

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        // scalar is order - 1
        let scalar = [0x43e1f593f0000000,
                    0x2833e84879b97091,
                    0xb85045b68181585d,
                    0x30644e72e131a029];
        
        b.iter(|| point.mul(&scalar).into_xy());
    }

    #[bench]
    fn bench_multiplication_bn254_into_affine_wnaf(b: &mut Bencher) {
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff);

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        // scalar is order - 1
        let scalar = U256Repr([0x43e1f593f0000000,
                    0x2833e84879b97091,
                    0xb85045b68181585d,
                    0x30644e72e131a029]);
        
        b.iter(|| point.wnaf_mul_impl(scalar).into_xy());
    }

    #[bench]
    fn bench_field_inverse(b: &mut Bencher) {
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let mut be_repr = vec![0u8; 32];
        be_repr[31] = 7u8;
        let element = Fp::from_be_bytes(&field, &be_repr[..], false).unwrap();
        
        b.iter(|| element.inverse().unwrap());
    }

    #[bench]
    fn bench_field_mont_inverse(b: &mut Bencher) {
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let mut be_repr = vec![0u8; 32];
        be_repr[31] = 7u8;
        let element = Fp::from_be_bytes(&field, &be_repr[..], false).unwrap();
        
        b.iter(|| element.mont_inverse().unwrap());
    }

    #[test]
    fn test_multiplication_bn254() {
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff);

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

    #[bench]
    fn bench_ben_coster_bn254(b: &mut Bencher) {
        use crate::representation::ElementRepr;
        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff);

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        let pairs: Vec<_> = (0..MULTIEXP_NUM_POINTS).map(|_| {
            let mut scalar = U256Repr::default();
            let mut bytes = vec![0u8; 32];
            rng.fill_bytes(&mut bytes[1..]);
            scalar.read_be(& bytes[..]).unwrap();

            (point.clone(), scalar)
        }).collect();

        b.iter(move || ben_coster(pairs.clone()));
    }

    #[bench]
    fn bench_naive_multiexp_bn254(b: &mut Bencher) {
        use crate::representation::ElementRepr;
        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff);

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        let pairs: Vec<_> = (0..MULTIEXP_NUM_POINTS).map(|_| {
            let mut scalar = U256Repr::default();
            let mut bytes = vec![0u8; 32];
            rng.fill_bytes(&mut bytes[1..]);
            scalar.read_be(& bytes[..]).unwrap();

            (point.clone(), scalar)
        }).collect();


        b.iter(move || {
            let mut pairs: Vec<_> = pairs.iter().map(|el| el.0.mul(el.1)).collect();
            let mut acc = pairs.pop().unwrap();
            while let Some(p) = pairs.pop() {
                acc.add_assign(&p);
            }
        });
    }

    #[bench]
    fn bench_wnaf_multiexp_bn254(b: &mut Bencher) {
        use crate::representation::ElementRepr;
        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff);

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        let pairs: Vec<_> = (0..MULTIEXP_NUM_POINTS).map(|_| {
            let mut scalar = U256Repr::default();
            let mut bytes = vec![0u8; 32];
            rng.fill_bytes(&mut bytes[1..]);
            scalar.read_be(& bytes[..]).unwrap();

            (point.clone(), scalar)
        }).collect();


        b.iter(move || {
            let mut pairs: Vec<_> = pairs.iter().map(|el| el.0.wnaf_mul_impl(el.1)).collect();
            let mut acc = pairs.pop().unwrap();
            while let Some(p) = pairs.pop() {
                acc.add_assign(&p);
            }
        });
    }

    #[test]
    fn test_ben_coster_bn254() {
        use crate::representation::ElementRepr;
        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff);

        let mut two = one.clone();
        two.double();

        let point = CurvePoint::point_from_xy(
            &curve, 
            one, 
            two);

        let pairs: Vec<_> = (0..MULTIEXP_NUM_POINTS).map(|_| {
            let mut scalar = U256Repr::default();
            let mut bytes = vec![0u8; 32];
            rng.fill_bytes(&mut bytes[1..]);
            scalar.read_be(& bytes[..]).unwrap();

            (point.clone(), scalar)
        }).collect();


        let naive_res = {
            let mut pairs: Vec<_> = pairs.iter().map(|el| el.0.mul(el.1)).collect();
            let mut acc = pairs.pop().unwrap();
            while let Some(p) = pairs.pop() {
                acc.add_assign(&p);
            }

            acc.into_xy()
        };

        let ben_coster_res = ben_coster(pairs).into_xy();

        assert!(ben_coster_res.0 == naive_res.0);
        assert!(ben_coster_res.1 == naive_res.1);
    }

    #[test]
    fn test_wnaf_decomposition() {
        use crate::representation::ElementRepr;
        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

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
        let one = Fp::one(&field);
        let a_coeff = Fp::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff);

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
        let wnaf_res = point.wnaf_mul_impl(scalar).into_xy();

        assert!(res_double_and_add.0 == wnaf_res.0);
        assert!(res_double_and_add.1 == wnaf_res.1);
    }
}