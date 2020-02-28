use crate::representation::*;
use crate::field::*;
use crate::fp::*;
use crate::weierstrass::*;
use crate::weierstrass::curve::*;
use crate::integers::MaxGroupSizeUint;
use crate::traits::*;

// define the test processor

pub(crate) struct ArithmeticProcessor<
    'a, 
    'b: 'a, 
    FE: FieldElement + ZeroAndOne + 'a,
    CP: CurveParameters<BaseFieldElement = FE> + 'a
> {
    curve: &'b WeierstrassCurve<'a, CP>,
    generator: &'b CurvePoint<'a, CP>,
    group_order: &'b [u64],
}

impl<
    'a, 
    'b: 'a, 
    FE: FieldElement + ZeroAndOne + 'a,
    CP: CurveParameters<BaseFieldElement = FE> + 'a
> ArithmeticProcessor<'a, 'b, FE, CP> {
    fn a_plus_a_equal_to_2a(&self) {
        let mut a_plus_a = self.generator.clone();
        let other_a = self.generator.clone();
        a_plus_a.add_assign(&other_a);

        let mut two_a = self.generator.clone();
        two_a.double();

        assert_eq!(a_plus_a.into_xy(), two_a.into_xy());
    }

    fn a_minus_a_equal_zero(&self) {
        let mut a_minus_a = self.generator.clone();
        let mut minus_a = self.generator.clone();
        minus_a.negate();

        a_minus_a.add_assign(&minus_a);

        assert!(a_minus_a.is_zero());
    }

    fn two_a_is_equal_to_two_a(&self) {
        let mut two_a = self.generator.clone();
        two_a.double();

        let other_two_a = self.generator.mul(&[2u64]);

        assert_eq!(other_two_a.into_xy(), two_a.into_xy());
    }

    fn a_plus_b_equal_to_b_plus_a(&self) {
        let mut b = self.generator.clone();
        b.double();

        let a = self.generator.clone();

        let mut a_plus_b = a.clone();
        a_plus_b.add_assign(&b);

        let mut b_plus_a = b.clone();
        b_plus_a.add_assign(&a);

        assert_eq!(a_plus_b.into_xy(), b_plus_a.into_xy());
    }

    fn a_mul_by_zero_is_zero(&self) {
        let a = self.generator.mul(&[0u64]);

        assert!(a.is_zero());
    }

    fn a_mul_by_group_order_is_zero(&self) {
        let a = self.generator.mul(&self.group_order);

        assert!(a.is_zero());
    }

    fn a_mul_by_sclar_wraps_over_group_order(&self) {
        let scalar = MaxGroupSizeUint::from(&[12345][..]);
        let group_order = MaxGroupSizeUint::from(&self.group_order[..]);
        let scalar_plus_group_order = scalar + group_order;
        let a = self.generator.mul(&scalar.as_ref());
        let b = self.generator.mul(&scalar_plus_group_order.as_ref());

        assert_eq!(a.into_xy(), b.into_xy());
    }

    pub fn test(&self) {
        self.a_minus_a_equal_zero();
        self.a_plus_a_equal_to_2a();
        self.two_a_is_equal_to_two_a();
        self.a_plus_b_equal_to_b_plus_a();
        self.a_mul_by_zero_is_zero();
        self.a_mul_by_group_order_is_zero();
        self.a_mul_by_sclar_wraps_over_group_order();
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::engines::bls12_381::*;

    #[test]
    fn test_bls12_381_g1() {
        let tester = ArithmeticProcessor::<_, _> {
            curve: &BLS12_381_PAIRING_ENGINE.curve,
            generator: &BLS12_381_G1_GENERATOR,
            group_order: &BLS12_381_SUBGROUP_ORDER,
        };

        tester.test();
    }

    #[test]
    fn test_bls12_381_g2() {
        let tester = ArithmeticProcessor::<_, _> {
            curve: &BLS12_381_PAIRING_ENGINE.curve_twist,
            generator: &BLS12_381_G2_GENERATOR,
            group_order: &BLS12_381_SUBGROUP_ORDER,
        };

        tester.test();
    }
}