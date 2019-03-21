#![feature(test)]

extern crate byteorder;
extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;

#[macro_use]
extern crate repr_derive;

mod arithmetics;
mod traits;
mod representation;
mod field;
mod weierstrass;
mod mont_inverse;

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
}