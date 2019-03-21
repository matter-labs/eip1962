#![feature(test)]

extern crate byteorder;
extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;
extern crate hex;
extern crate rand;
extern crate rand_xorshift;

#[macro_use]
extern crate repr_derive;

mod arithmetics;
mod traits;
mod representation;
mod field;
mod weierstrass;
mod mont_inverse;
mod multiexp;

pub use representation::ElementRepr;

extern crate test;

pub fn add_two(a: i32) -> i32 {
    a + 2
}

#[cfg(test)]
mod tests {
    use crate::field::*;
    use crate::weierstrass::*;
    use crate::traits::FieldElement;
    use test::Bencher;
    use crate::multiexp::ben_coster;

    const MULTIEXP_NUM_POINTS: usize = 100;

    #[bench]
    fn bench_doubling_bn254(b: &mut Bencher) {
        let field = new_field("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = PrimeFieldElement::one(&field);
        let a_coeff = PrimeFieldElement::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff, 
            CurveType::Generic);

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
        let field = new_field("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = PrimeFieldElement::one(&field);
        let a_coeff = PrimeFieldElement::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff, 
            CurveType::Generic);

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
        let field = new_field("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = PrimeFieldElement::one(&field);
        let a_coeff = PrimeFieldElement::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff, 
            CurveType::Generic);

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
        let field = new_field("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = PrimeFieldElement::one(&field);
        let a_coeff = PrimeFieldElement::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff, 
            CurveType::Generic);

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
    fn bench_field_inverse(b: &mut Bencher) {
        let field = new_field("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let mut be_repr = vec![0u8; 32];
        be_repr[31] = 7u8;
        let element = PrimeFieldElement::from_be_bytes(&field, &be_repr[..]).unwrap();
        
        b.iter(|| element.inverse().unwrap());
    }

    #[bench]
    fn bench_field_mont_inverse(b: &mut Bencher) {
        let field = new_field("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let mut be_repr = vec![0u8; 32];
        be_repr[31] = 7u8;
        let element = PrimeFieldElement::from_be_bytes(&field, &be_repr[..]).unwrap();
        
        b.iter(|| element.mont_inverse().unwrap());
    }

    #[test]
    fn test_multiplication_bn254() {
        let field = new_field("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = PrimeFieldElement::one(&field);
        let a_coeff = PrimeFieldElement::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff, 
            CurveType::Generic);

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
        let field = new_field("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = PrimeFieldElement::one(&field);
        let a_coeff = PrimeFieldElement::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff, 
            CurveType::Generic);

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
        let field = new_field("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = PrimeFieldElement::one(&field);
        let a_coeff = PrimeFieldElement::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff, 
            CurveType::Generic);

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

    #[test]
    fn test_ben_coster_bn254() {
        use crate::representation::ElementRepr;
        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
        let field = new_field("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let group = new_field("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let one = PrimeFieldElement::one(&field);
        let a_coeff = PrimeFieldElement::zero(&field);
        let mut b_coeff = one.clone();
        b_coeff.double();
        b_coeff.add_assign(&one);

        let curve = WeierstrassCurve::new(
            &group, 
            a_coeff, 
            b_coeff, 
            CurveType::Generic);

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
}