
use crate::field::SizedPrimeField;
use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, BitIterator};
use crate::weierstrass::Group;
use crate::weierstrass::curve::CurvePoint;
use crate::weierstrass::twist::TwistPoint;
use crate::extension_towers::fp2::{Fp2, Extension2};
use crate::extension_towers::fp6_as_3_over_2::{Fp6, Extension3Over2};
use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};
use num_bigint::BigUint;
use num_traits::FromPrimitive;
use num_integer::Integer;
use num_traits::Zero;

pub mod bls12;

pub trait PairingEngine: Sized {
    type PairingResult: FieldElement;
    type G1: Group;
    type G2: Group;

    // fn prepare_twist_point<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>>
    //     (&self, twist_point: &'a TwistPoint<'a, FE, F, GE, G>) -> PreparedTwistPoint<'a, FE, F>;

    fn pair<'b> (&self, points: &'b [Self::G1], twists: &'b [Self::G2]) -> Option<Self::PairingResult>;
    
    // fn pair<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>>
        // (&self, point: &'a CurvePoint<'a, FE, F, GE, G>, twist_point: &'a TwistPoint<'a, FE, F, GE, G>) -> Self::PairingResult;
}

pub fn frobenius_calculator_fp2<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>(
        // modulus: BigUint,
        // base_field: &'a F,
        extension: &Extension2<'a, FE, F>,
    ) -> Result<[Fp<'a, FE, F>; 2], ()> {
        // use crate::field::biguint_to_u64_vec;
        let non_residue = extension.non_residue.clone();

        // NONRESIDUE**(((q^0) - 1) / 2)
        let f_0 = Fp::one(extension.field);

        let mut two_inv = Fp::one(extension.field);
        two_inv.double();
        let two_inv = two_inv.inverse().ok_or(())?;

        // NONRESIDUE**(((q^1) - 1) / 2)
        let mut power = Fp::one(extension.field);
        power.negate();
        power.mul_assign(&two_inv);
        let f_1 = non_residue.pow(&power.into_repr());

        // let one = BigUint::from_u64(1).unwrap();
        // let two = BigUint::from_u64(2).unwrap();

        // // Fq2(u + 1)**(((q^1) - 1) / 3)
        // let mut q_power = modulus.clone();
        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&two);
        // debug_assert!(rem.is_zero());
        // let f_1 = non_residue.pow(&biguint_to_u64_vec(power));

        Ok([f_0, f_1])
}

pub fn frobenius_calculator_fp6<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>(
        modulus: BigUint,
        extension: &Extension3Over2<'a, FE, F>
    ) -> Result<([Fp2<'a, FE, F>; 6], [Fp2<'a, FE, F>; 6]), ()> {
        use crate::field::biguint_to_u64_vec;
        // let mut two = Fp::one(field);
        // two.double();
        // let two_inv = two.inverse().ok_or(())?;

        // let mut three = Fp::one(field);
        // three.add_assign(&two);
        // let three_inv = three.inverse().ok_or(())?;

        // let mut two_le_bytes = vec![];
        // two.into_repr().write_le(&mut two_le_bytes).map_err(|_| Err())?;

        // let mut three_le_bytes = vec![];
        // three.into_repr().write_le(&mut three_le_bytes).map_err(|_| Err())?;

        let one = BigUint::from_u64(1).unwrap();
        let three = BigUint::from_u64(3).unwrap();
        // let two_inv = BigUint::from_bytes_le(&two_le_bytes);
        // let three_inv = BigUint::from_bytes_le(&three_le_bytes);
    
        // Fq2(u + 1)**(((q^0) - 1) / 3)
        let non_residue = extension.non_residue.clone();
        let f_0 = Fp2::one(extension.field);

        // Fq2(u + 1)**(((q^1) - 1) / 3)
        let mut q_power = modulus.clone();
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&three);
        debug_assert!(rem.is_zero());
        let f_1 = non_residue.pow(&biguint_to_u64_vec(power));

        // Fq2(u + 1)**(((q^2) - 1) / 3)
        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&three);
        debug_assert!(rem.is_zero());
        let f_2 = non_residue.pow(&biguint_to_u64_vec(power));

        // Fq2(u + 1)**(((q^3) - 1) / 3)
        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&three);
        debug_assert!(rem.is_zero());
        let f_3 = non_residue.pow(&biguint_to_u64_vec(power));

        // Fq2(u + 1)**(((q^4) - 1) / 3)
        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&three);
        debug_assert!(rem.is_zero());
        let f_4 = non_residue.pow(&biguint_to_u64_vec(power));

        // Fq2(u + 1)**(((q^5) - 1) / 3)
        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&three);
        debug_assert!(rem.is_zero());
        let f_5 = non_residue.pow(&biguint_to_u64_vec(power));



        let f_0_c2 = f_0.clone();

        let mut f_1_c2 = f_1.clone();
        let mut f_2_c2 = f_2.clone();
        let mut f_3_c2 = f_3.clone();
        let mut f_4_c2 = f_4.clone();
        let mut f_5_c2 = f_5.clone();

        f_1_c2.square();
        f_2_c2.square();
        f_3_c2.square();
        f_4_c2.square();
        f_5_c2.square();

        Ok(([f_0, f_1, f_2, f_3, f_4, f_5], [f_0_c2, f_1_c2, f_2_c2, f_3_c2, f_4_c2, f_5_c2]))
}

pub fn frobenius_calculator_fp12<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>(
        modulus: BigUint,
        extension: &Extension2Over3Over2<'a, FE, F>
    ) -> Result<[Fp2<'a, FE, F>; 12], ()> {
        use crate::field::biguint_to_u64_vec;

        let one = BigUint::from_u64(1).unwrap();
        let six = BigUint::from_u64(6).unwrap();
    
        // Fq2(u + 1)**(((q^0) - 1) / 6)
        let non_residue = extension.field.non_residue.clone();
        let f_0 = Fp2::one(extension.field.field);

        // Fq2(u + 1)**(((q^1) - 1) / 6)
        let mut q_power = modulus.clone();
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_1 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_2 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_3 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_4 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_5 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_6 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_7 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_8 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_9 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_10 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_11 = non_residue.pow(&biguint_to_u64_vec(power));

        Ok([f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11])
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use num_traits::FromPrimitive;
    use num_integer::Integer;
    use num_traits::Zero;
    use crate::field::{U384Repr, new_field};
    use crate::fp::Fp;
    use crate::traits::{FieldElement};
    use crate::extension_towers::fp2::{Fp2, Extension2};
    use crate::extension_towers::fp6_as_3_over_2::{Fp6, Extension3Over2};
    use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};
    use num_traits::Num;

    #[test]
    fn test_bls12_381_frob_fp2() {
        let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let extension = Extension2 {
            field: &base_field,
            non_residue: fp_non_residue,
            frobenius_coeffs_c1: [Fp::zero(&base_field), Fp::zero(&base_field)]
        };

        let coeffs = super::frobenius_calculator_fp2(&extension).unwrap();


        // let coeffs = super::frobenius_calculator_fp2(&base_field, &extension).unwrap();

        println!("C_0 = {}", coeffs[0]);
        println!("C_1 = {}", coeffs[1]);


        println!("C_0_1 = {:x}", coeffs[0].repr.0[0]); 
    }

    #[test]
    fn test_bls12_381_frob_fp6() {
        let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let mut extension_2 = Extension2 {
            field: &base_field,
            non_residue: fp_non_residue,
            frobenius_coeffs_c1: [Fp::zero(&base_field), Fp::zero(&base_field)]
        };

        let coeffs = super::frobenius_calculator_fp2(&extension_2).unwrap();
        extension_2.frobenius_coeffs_c1 = coeffs;

        let one = Fp::one(&base_field);

        let mut fp2_non_residue = Fp2::zero(&extension_2); // non-residue is 1 + u
        fp2_non_residue.c0 = one.clone();
        fp2_non_residue.c1 = one.clone();

        let f_c1 = [Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2)];

        let extension_6 = Extension3Over2 {
            non_residue: fp2_non_residue,
            field: &extension_2,
            frobenius_coeffs_c1: f_c1.clone(),
            frobenius_coeffs_c2: f_c1,
        };

        let coeffs = super::frobenius_calculator_fp6(modulus, &extension_6).unwrap();

        println!("C_0_0 = {}", coeffs.0[0]);
        println!("C_0_1 = {}", coeffs.0[1]);

        println!("C_1_0 = {}", coeffs.1[0]);
        println!("C_1_1 = {}", coeffs.1[1]);


        // println!("C_0_1 = {:x}", coeffs[0].repr.0[0]); 
    }

    #[test]
    fn test_bls12_381_frob_fp12() {
        let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let mut extension_2 = Extension2 {
            field: &base_field,
            non_residue: fp_non_residue,
            frobenius_coeffs_c1: [Fp::zero(&base_field), Fp::zero(&base_field)]
        };

        let coeffs = super::frobenius_calculator_fp2(&extension_2).unwrap();
        extension_2.frobenius_coeffs_c1 = coeffs;

        let one = Fp::one(&base_field);

        let mut fp2_non_residue = Fp2::zero(&extension_2);
        fp2_non_residue.c0 = one.clone();
        fp2_non_residue.c1 = one.clone();

        let f_c1 = [Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2)];

        let mut extension_6 = Extension3Over2 {
            non_residue: fp2_non_residue,
            field: &extension_2,
            frobenius_coeffs_c1: f_c1.clone(),
            frobenius_coeffs_c2: f_c1,
        };

        let (coeffs_c1, coeffs_c2) = super::frobenius_calculator_fp6(modulus.clone(), &extension_6).unwrap();

        extension_6.frobenius_coeffs_c1 = coeffs_c1;
        extension_6.frobenius_coeffs_c2 = coeffs_c2;

        let mut fp2_non_residue = Fp2::zero(&extension_2);

         let f_c1 = [Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2)];

        let mut extension_12 = Extension2Over3Over2 {
            non_residue: Fp6::zero(&extension_6),
            field: &extension_6,
            frobenius_coeffs_c1: f_c1,
        };

        let coeffs = super::frobenius_calculator_fp12(modulus, &extension_12).unwrap();

        println!("C_0 = {}", coeffs[0]);
        println!("C_1 = {}", coeffs[1]);
        println!("C_10 = {}", coeffs[10]);
    }


}
