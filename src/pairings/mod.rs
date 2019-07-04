
use crate::field::SizedPrimeField;
use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, ZeroAndOne};
use crate::weierstrass::Group;
use crate::extension_towers::{fp2::Fp2, fp2::Extension2};
use crate::extension_towers::{fp3::Fp3, fp3::Extension3};
use crate::extension_towers::fp4_as_2_over_2;
use crate::extension_towers::fp6_as_2_over_3;
use crate::extension_towers::fp6_as_3_over_2;
use crate::extension_towers::{fp12_as_2_over3_over_2::Fp12, fp12_as_2_over3_over_2::Extension2Over3Over2};
use num_bigint::BigUint;
use num_traits::FromPrimitive;
use num_integer::Integer;
use num_traits::Zero;

pub(crate) mod bls12;
pub(crate) mod bn;
// pub(crate) mod cp;
pub(crate) mod mnt6;
pub(crate) mod mnt4;

#[derive(Eq, PartialEq, Clone, Copy, Debug)]
pub(crate) enum TwistType {
    D,
    M
}

pub(crate) trait PairingEngine: Sized {
    type PairingResult: FieldElement;
    type G1: Group;
    type G2: Group;

    fn pair<'b> (&self, points: &'b [Self::G1], twists: &'b [Self::G2]) -> Option<Self::PairingResult>;
}

pub(crate) fn frobenius_calculator_fp2<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>(
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

pub(crate) fn frobenius_calculator_fp3<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>(
        modulus: BigUint,
        extension: &Extension3<'a, FE, F>
    ) -> Result<([Fp<'a, FE, F>; 3], [Fp<'a, FE, F>; 3]), ()> {
        use crate::field::biguint_to_u64_vec;
        let one = BigUint::from_u64(1).unwrap();
        let three = BigUint::from_u64(3).unwrap();
    
        // NON_RESIDUE**(((q^0) - 1) / 3)
        let non_residue = extension.non_residue.clone();
        let f_0 = Fp::one(extension.field);

        // NON_RESIDUE**(((q^1) - 1) / 3)
        let mut q_power = modulus.clone();
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&three);
        debug_assert!(rem.is_zero());
        let f_1 = non_residue.pow(&biguint_to_u64_vec(power));

        // NON_RESIDUE**(((q^2) - 1) / 3)
        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&three);
        debug_assert!(rem.is_zero());
        let f_2 = non_residue.pow(&biguint_to_u64_vec(power));


        let f_0_c2 = f_0.clone();

        let mut f_1_c2 = f_1.clone();
        let mut f_2_c2 = f_2.clone();

        f_1_c2.square();
        f_2_c2.square();

        Ok(([f_0, f_1, f_2], [f_0_c2, f_1_c2, f_2_c2]))
}

pub(crate) fn frobenius_calculator_fp4_as_2_over_2<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>(
        modulus: BigUint,
        extension: &fp4_as_2_over_2::Extension2Over2<'a, FE, F>
    ) -> Result<[Fp<'a, FE, F>; 4], ()> {
        use crate::field::biguint_to_u64_vec;

        // TODO: change into bitshift

        let one = BigUint::from_u64(1).unwrap();
        // let divisor = BigUint::from_u64(4).unwrap();
    
        // NON_REDISUE**(((q^0) - 1) / 4)
        let non_residue = extension.field.non_residue.clone();
        let f_0 = Fp::one(extension.field.field);

        // NON_REDISUE**(((q^1) - 1) / 4)
        let mut q_power = modulus.clone();
        let power = (q_power.clone() - &one) >> 2;
        // let (power, rem) = power.div_rem(&divisor);
        // debug_assert!(rem.is_zero());
        let f_1 = non_residue.pow(&biguint_to_u64_vec(power));

        // NON_REDISUE**(((q^2) - 1) / 4)
        q_power *= &modulus;
        let power = (q_power.clone() - &one) >> 2;
        // let (power, rem) = power.div_rem(&divisor);
        // debug_assert!(rem.is_zero());
        let f_2 = non_residue.pow(&biguint_to_u64_vec(power));

        // this one is not needed
        // // NON_REDISUE**(((q^3) - 1) / 4)
        // q_power *= &modulus;
        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&divisor);
        // debug_assert!(rem.is_zero());
        // let f_3 = non_residue.pow(&biguint_to_u64_vec(power));

        let f_3 = Fp::zero(extension.field.field);

        Ok([f_0, f_1, f_2, f_3])
}

pub(crate) fn frobenius_calculator_fp6_as_2_over_3<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>(
        modulus: BigUint,
        extension: &fp6_as_2_over_3::Extension2Over3<'a, FE, F>
    ) -> Result<[Fp<'a, FE, F>; 6], ()> {
        use crate::field::biguint_to_u64_vec;

        let one = BigUint::from_u64(1).unwrap();
        let divisor = BigUint::from_u64(6).unwrap();
    
        // NON_REDISUE**(((q^0) - 1) / 6)
        let non_residue = extension.field.non_residue.clone();
        let f_0 = Fp::one(extension.field.field);

        let mut q_power = modulus.clone();
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&divisor);
        debug_assert!(rem.is_zero());
        let f_1 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let f_2 = Fp::zero(extension.field.field);
        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&divisor);
        // debug_assert!(rem.is_zero());
        // let f_2 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&divisor);
        debug_assert!(rem.is_zero());
        let f_3 = non_residue.pow(&biguint_to_u64_vec(power));

        let f_4 = Fp::zero(extension.field.field);
        let f_5 = Fp::zero(extension.field.field);

        // not needed
        // q_power *= &modulus;
        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&divisor);
        // debug_assert!(rem.is_zero());
        // let f_4 = non_residue.pow(&biguint_to_u64_vec(power));

        // q_power *= &modulus;
        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&divisor);
        // debug_assert!(rem.is_zero());
        // let f_5 = non_residue.pow(&biguint_to_u64_vec(power));

        Ok([f_0, f_1, f_2, f_3, f_4, f_5])
}

pub(crate) fn frobenius_calculator_fp6_as_3_over_2<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>(
        modulus: BigUint,
        extension: &fp6_as_3_over_2::Extension3Over2<'a, FE, F>
    ) -> Result<([Fp2<'a, FE, F>; 6], [Fp2<'a, FE, F>; 6]), ()> {
        use crate::field::biguint_to_u64_vec;

        let one = BigUint::from_u64(1).unwrap();
        let three = BigUint::from_u64(3).unwrap();

        // NON_RESIDUE**(((q^0) - 1) / 3)
        let non_residue = extension.non_residue.clone();
        let f_0 = Fp2::one(extension.field);

        let mut q_power = modulus.clone();
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&three);
        debug_assert!(rem.is_zero());
        let f_1 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&three);
        debug_assert!(rem.is_zero());
        let f_2 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&three);
        debug_assert!(rem.is_zero());
        let f_3 = non_residue.pow(&biguint_to_u64_vec(power));

        // q_power *= &modulus;
        let f_4 = Fp2::zero(extension.field);

        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&three);
        // debug_assert!(rem.is_zero());
        // let f_4 = non_residue.pow(&biguint_to_u64_vec(power));

        // q_power *= &modulus;
        let f_5 = Fp2::zero(extension.field);

        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&three);
        // debug_assert!(rem.is_zero());
        // let f_5 = non_residue.pow(&biguint_to_u64_vec(power));

        let f_0_c2 = f_0.clone();

        let mut f_1_c2 = f_1.clone();
        f_1_c2.square();
        let mut f_2_c2 = f_2.clone();
        f_2_c2.square();
        let mut f_3_c2 = f_3.clone();
        f_3_c2.square();

        let f_4_c2 = f_4.clone();
        let f_5_c2 = f_5.clone();

        Ok(([f_0, f_1, f_2, f_3, f_4, f_5], [f_0_c2, f_1_c2, f_2_c2, f_3_c2, f_4_c2, f_5_c2]))
}


use crate::sliding_window_exp::{WindowExpBase};
pub(crate) fn frobenius_calculator_fp6_as_3_over_2_using_sliding_window<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>(
        modulus: BigUint,
        base: &WindowExpBase<Fp2<'a, FE, F>>,
        extension: &fp6_as_3_over_2::Extension3Over2<'a, FE, F>
    ) -> Result<([Fp2<'a, FE, F>; 6], [Fp2<'a, FE, F>; 6]), ()> {
        use crate::field::biguint_to_u64_vec;

        let one = BigUint::from_u64(1).unwrap();
        let three = BigUint::from_u64(3).unwrap();

        // NON_RESIDUE**(((q^0) - 1) / 3)
        // let non_residue = extension.non_residue.clone();
        let f_0 = Fp2::one(extension.field);

        let mut powers = vec![];

        let mut q_power = modulus.clone();

        {
            let power = q_power.clone() - &one;
            let (power, rem) = power.div_rem(&three);
            debug_assert!(rem.is_zero());
            powers.push(biguint_to_u64_vec(power));
        }
        for _ in 1..3 {
            q_power *= &modulus;
            let power = q_power.clone() - &one;
            let (power, rem) = power.div_rem(&three);
            debug_assert!(rem.is_zero());
            powers.push(biguint_to_u64_vec(power));
        }

        let mut result = base.exponentiate(&powers);
        debug_assert!(result.len() == 3);

        let f_3 = result.pop().unwrap();
        let f_2 = result.pop().unwrap();
        let f_1 = result.pop().unwrap();

        let f_4 = Fp2::zero(extension.field);
        let f_5 = Fp2::zero(extension.field);

        let f_0_c2 = f_0.clone();

        let mut f_1_c2 = f_1.clone();
        f_1_c2.square();
        let mut f_2_c2 = f_2.clone();
        f_2_c2.square();
        let mut f_3_c2 = f_3.clone();
        f_3_c2.square();

        let f_4_c2 = f_4.clone();
        let f_5_c2 = f_5.clone();

        Ok(([f_0, f_1, f_2, f_3, f_4, f_5], [f_0_c2, f_1_c2, f_2_c2, f_3_c2, f_4_c2, f_5_c2]))
}

pub(crate) fn frobenius_calculator_fp12_using_sliding_window<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>(
        modulus: BigUint,
        base: &WindowExpBase<Fp2<'a, FE, F>>,
        extension: &Extension2Over3Over2<'a, FE, F>
    ) -> Result<[Fp2<'a, FE, F>; 12], ()> {
        use crate::field::biguint_to_u64_vec;

        let one = BigUint::from_u64(1).unwrap();
        let six = BigUint::from_u64(6).unwrap();
    
        // Fq2(u + 1)**(((q^0) - 1) / 6)
        let f_0 = Fp2::one(extension.field.field);

        // Fq2(u + 1)**(((q^1) - 1) / 6)
        let mut powers = vec![];

        let mut q_power = modulus.clone();

        // 1
        {
            let power = q_power.clone() - &one;
            let (power, rem) = power.div_rem(&six);
            debug_assert!(rem.is_zero());
            powers.push(biguint_to_u64_vec(power));
        }
        // 2 & 3
        for _ in 1..3 {
            q_power *= &modulus;
            let power = q_power.clone() - &one;
            let (power, rem) = power.div_rem(&six);
            debug_assert!(rem.is_zero());
            powers.push(biguint_to_u64_vec(power));
        }
        // 6
        {
            q_power *= q_power.clone();
            let power = q_power.clone() - &one;
            let (power, rem) = power.div_rem(&six);
            debug_assert!(rem.is_zero());
            powers.push(biguint_to_u64_vec(power));
        }

        let mut result = base.exponentiate(&powers);
        debug_assert!(result.len() == 4);

        let f_6 = result.pop().unwrap();
        let f_3 = result.pop().unwrap();
        let f_2 = result.pop().unwrap();
        let f_1 = result.pop().unwrap();

        let f_4 = Fp2::zero(extension.field.field);
        let f_5 = Fp2::zero(extension.field.field);

        let f_7 = Fp2::zero(extension.field.field);
        let f_8 = Fp2::zero(extension.field.field);
        let f_9 = Fp2::zero(extension.field.field);
        let f_10 = Fp2::zero(extension.field.field);
        let f_11 = Fp2::zero(extension.field.field);

        Ok([f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11])
}


pub(crate) fn frobenius_calculator_fp12<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>(
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
        let f_4 = Fp2::zero(extension.field.field);

        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&six);
        // debug_assert!(rem.is_zero());
        // let f_4 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let f_5 = Fp2::zero(extension.field.field);
        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&six);
        // debug_assert!(rem.is_zero());
        // let f_5 = non_residue.pow(&biguint_to_u64_vec(power));

        q_power *= &modulus;
        let power = q_power.clone() - &one;
        let (power, rem) = power.div_rem(&six);
        debug_assert!(rem.is_zero());
        let f_6 = non_residue.pow(&biguint_to_u64_vec(power));

        let f_7 = Fp2::zero(extension.field.field);
        let f_8 = Fp2::zero(extension.field.field);
        let f_9 = Fp2::zero(extension.field.field);
        let f_10 = Fp2::zero(extension.field.field);
        let f_11 = Fp2::zero(extension.field.field);

        // q_power *= &modulus;
        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&six);
        // debug_assert!(rem.is_zero());
        // let f_7 = non_residue.pow(&biguint_to_u64_vec(power));

        // q_power *= &modulus;
        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&six);
        // debug_assert!(rem.is_zero());
        // let f_8 = non_residue.pow(&biguint_to_u64_vec(power));

        // q_power *= &modulus;
        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&six);
        // debug_assert!(rem.is_zero());
        // let f_9 = non_residue.pow(&biguint_to_u64_vec(power));

        // q_power *= &modulus;
        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&six);
        // debug_assert!(rem.is_zero());
        // let f_10 = non_residue.pow(&biguint_to_u64_vec(power));

        // q_power *= &modulus;
        // let power = q_power.clone() - &one;
        // let (power, rem) = power.div_rem(&six);
        // debug_assert!(rem.is_zero());
        // let f_11 = non_residue.pow(&biguint_to_u64_vec(power));

        Ok([f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11])
}

pub(crate) fn into_ternary_wnaf(repr: &[u64]) -> Vec<i64> {
    fn is_zero(repr: &[u64]) -> bool {

        for el in repr.iter() {
            if *el != 0 {
                return false;
            }
        }

        true
    }

    fn is_odd(repr: &[u64]) -> bool {
        if repr.len() == 0 {
            return false;
        }

        repr[0] & 1u64 == 1u64
    }

    fn div2(repr: &mut [u64]) {
        let mut t = 0;
        for i in repr.iter_mut().rev() {
            let t2 = *i << 63;
            *i >>= 1;
            *i |= t;
            t = t2;
        }
    }

    let mut res = vec![];
    let mut e = repr.to_vec();
    while !is_zero(&e) {
        let z: i64;
        if is_odd(&e) {
            z = 2 - (e[0] % 4) as i64;
            if z >= 0 {
                e[0] -= z as u64;
            } else {
                e[0] += -z as u64;
            }
        } else {
            z = 0;
        }
        res.push(z);
        div2(&mut e);
    }

    res
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use crate::field::{U384Repr, U832Repr, new_field};
    use crate::fp::Fp;
    use crate::traits::{FieldElement, ZeroAndOne};
    use crate::extension_towers::fp2::{Fp2, Extension2};
    use crate::extension_towers::fp3::{Fp3, Extension3};
    use crate::extension_towers::fp6_as_2_over_3;
    use crate::extension_towers::fp6_as_3_over_2;
    use crate::extension_towers::fp12_as_2_over3_over_2::{Extension2Over3Over2};
    use num_traits::Num;

    #[test]
    fn test_bls12_381_frob_fp2() {
        let _modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let extension = Extension2::new(fp_non_residue);

        let coeffs = super::frobenius_calculator_fp2(&extension).unwrap();
    }

    #[test]
    fn test_frob_fp3_calculator() {
        let modulus = BigUint::from_str_radix("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
        let base_field = new_field::<U832Repr>("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
        let nonres_repr = U832Repr::from(13);
        let fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

        let extension = Extension3::new(fp_non_residue);

        let (coeffs_0, coeffs_1) = super::frobenius_calculator_fp3(modulus, &extension).unwrap();
    }

    #[test]
    fn test_frob_fp6_as_2_over_3() {
        let modulus = BigUint::from_str_radix("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
        let base_field = new_field::<U832Repr>("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
        let nonres_repr = U832Repr::from(13);
        let fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

        let mut extension_3 = Extension3::new(fp_non_residue.clone());

        let (coeffs_1, coeffs_2) = super::frobenius_calculator_fp3(modulus.clone(), &extension_3).unwrap();
        extension_3.frobenius_coeffs_c1 = coeffs_1;
        extension_3.frobenius_coeffs_c2 = coeffs_2;
        extension_3.frobenius_coeffs_are_calculated = true;

        let mut fp3_non_residue = Fp3::zero(&extension_3); // non-residue is 13 + 0*u + 0*u^2
        fp3_non_residue.c0 = fp_non_residue;

        let extension_6 = fp6_as_2_over_3::Extension2Over3::new(fp3_non_residue);

        let coeffs = super::frobenius_calculator_fp6_as_2_over_3(modulus, &extension_6).unwrap();

    }

    #[test]
    fn test_bls12_381_frob_fp6() {
        let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let mut extension_2 = Extension2::new(fp_non_residue);

        let coeffs = super::frobenius_calculator_fp2(&extension_2).unwrap();
        extension_2.frobenius_coeffs_c1 = coeffs;

        let one = Fp::one(&base_field);

        let mut fp2_non_residue = Fp2::zero(&extension_2); // non-residue is 1 + u
        fp2_non_residue.c0 = one.clone();
        fp2_non_residue.c1 = one.clone();

        let extension_6 = fp6_as_3_over_2::Extension3Over2::new(fp2_non_residue);

        let coeffs = super::frobenius_calculator_fp6_as_3_over_2(modulus, &extension_6).unwrap();
    }

    #[test]
    fn test_bls12_381_frob_fp12() {
        let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let mut extension_2 = Extension2::new(fp_non_residue);

        let coeffs = super::frobenius_calculator_fp2(&extension_2).unwrap();
        extension_2.frobenius_coeffs_c1 = coeffs;

        let one = Fp::one(&base_field);

        let mut fp2_non_residue = Fp2::zero(&extension_2);
        fp2_non_residue.c0 = one.clone();
        fp2_non_residue.c1 = one.clone();

        let mut extension_6 = fp6_as_3_over_2::Extension3Over2::new(fp2_non_residue);

        let (coeffs_c1, coeffs_c2) = super::frobenius_calculator_fp6_as_3_over_2(modulus.clone(), &extension_6).unwrap();

        extension_6.frobenius_coeffs_c1 = coeffs_c1;
        extension_6.frobenius_coeffs_c2 = coeffs_c2;

        let mut fp2_non_residue = Fp2::zero(&extension_2);

        let mut extension_12 = Extension2Over3Over2::new(fp6_as_3_over_2::Fp6::zero(&extension_6));

        let coeffs = super::frobenius_calculator_fp12(modulus, &extension_12).unwrap();
    }
}
