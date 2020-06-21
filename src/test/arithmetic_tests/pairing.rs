use crate::weierstrass::*;
use crate::weierstrass::curve::*;
use crate::integers::MaxGroupSizeUint;
use crate::traits::*;
use crate::pairings::PairingEngine;
use rand::Rng;

// define the test processor
pub(crate) struct PairingProcessor<
    'b, 
    'a: 'b, 
    FE: FieldElement + ZeroAndOne + 'a,
    FEE: FieldElement + ZeroAndOne + 'a,
    GT: FieldElement + ZeroAndOne + 'a,
    CP: CurveParameters<BaseFieldElement = FE> + 'a,
    CPT: CurveParameters<BaseFieldElement = FEE> + 'a,
    E: PairingEngine<G1 = CurvePoint<'a, CP>, G2 = CurvePoint<'a, CPT>, PairingResult = GT> + 'a
> {
    curve_g1: &'b WeierstrassCurve<'a, CP>,
    curve_g2: &'b WeierstrassCurve<'a, CPT>,
    generator_g1: &'b CurvePoint<'a, CP>,
    generator_g2: &'b CurvePoint<'a, CPT>,
    group_order: &'b [u64],
    gt_one: &'b GT,
    engine: &'b E
}

impl<
    'b, 
    'a: 'b, 
    FE: FieldElement + ZeroAndOne + 'a,
    FEE: FieldElement + ZeroAndOne + 'a,
    GT: FieldElement + ZeroAndOne + 'a,
    CP: CurveParameters<BaseFieldElement = FE> + 'a,
    CPT: CurveParameters<BaseFieldElement = FEE> + 'a,
    E: PairingEngine<G1 = CurvePoint<'a, CP>, G2 = CurvePoint<'a, CPT>, PairingResult = GT> + 'a
> PairingProcessor<'b, 'a, FE, FEE, GT, CP, CPT, E> {
    fn bilinearity<R: Rng>(&self, rng: &mut R, num_attempts: usize, scalar_words: usize) {
        for _ in 0..num_attempts {
            let mut limbs = vec![];
            for _ in 0..scalar_words {
                let limb: u64 = rng.gen();
                limbs.push(limb);
            }
            let scalar = MaxGroupSizeUint::from(&limbs[..]);
            let mut g1_by_scalar = self.generator_g1.mul(&scalar.as_ref());
            g1_by_scalar.normalize();
            let g1 = self.generator_g1.clone();

            let mut g2_by_scalar = self.generator_g2.mul(&scalar.as_ref());
            g2_by_scalar.normalize();
            let g2 = self.generator_g2.clone();

            let result_ab = self.engine.pair(&[g1_by_scalar], &[g2]).unwrap();
            let result_ba = self.engine.pair(&[g1], &[g2_by_scalar]).unwrap();

            assert_eq!(result_ab, result_ba);
        }
    }

    fn degeneracy(&self) {
        let zero_g1 = CurvePoint::zero(self.generator_g1.curve);
        let zero_g2 = CurvePoint::zero(self.generator_g2.curve);

        let scalar = MaxGroupSizeUint::from(&[12345][..]);

        let mut g1 = self.generator_g1.mul(&scalar.as_ref());
        g1.normalize();
        let mut g2 = self.generator_g2.mul(&scalar.as_ref());
        g2.normalize();

        let result_ab = self.engine.pair(&[zero_g1], &[g2]).unwrap();
        let result_ba = self.engine.pair(&[g1], &[zero_g2]).unwrap();

        let one = self.gt_one;
        assert_eq!(result_ab, result_ba);
        assert_eq!(&result_ab, one);
    }

    pub fn test<R: Rng>(&self, rng: &mut R, num_attempts: usize, scalar_words: usize) {
        self.bilinearity(rng, num_attempts, scalar_words);
        self.degeneracy();
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use crate::engines::bls12_381::*;
    use crate::engines::bls12_377::*;

    use crate::extension_towers::fp12_as_2_over3_over_2::Fp12;

    #[test]
    fn test_bls12_381_pairing() {
        let gt_one = Fp12::one(&BLS12_381_EXTENSION_12_FIELD);
        let tester = PairingProcessor::<_, _, _, _, _, _> {
            curve_g1: &BLS12_381_G1_CURVE,
            curve_g2: &BLS12_381_G2_CURVE,
            generator_g1: &BLS12_381_G1_GENERATOR,
            generator_g2: &BLS12_381_G2_GENERATOR,
            group_order: &BLS12_381_SUBGROUP_ORDER,
            gt_one: &gt_one,
            engine: &BLS12_381_PAIRING_ENGINE
        };

        tester.test(&mut make_rng(), 1000, 4);
    }

    #[test]
    fn test_bls12_377_pairing() {
        let gt_one = Fp12::one(&BLS12_377_EXTENSION_12_FIELD);
        let tester = PairingProcessor::<_, _, _, _, _, _> {
            curve_g1: &BLS12_377_G1_CURVE,
            curve_g2: &BLS12_377_G2_CURVE,
            generator_g1: &BLS12_377_G1_GENERATOR,
            generator_g2: &BLS12_377_G2_GENERATOR,
            group_order: &BLS12_377_SUBGROUP_ORDER,
            gt_one: &gt_one,
            engine: &BLS12_377_PAIRING_ENGINE
        };

        tester.test(&mut make_rng(), 1000, 4);
    }

    fn make_rng() -> rand_xorshift::XorShiftRng {
        use rand::{SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        rng
    }
}