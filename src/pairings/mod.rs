
use crate::traits::{FieldElement};
use crate::weierstrass::Group;

pub(crate) mod bls12;
pub(crate) mod bn;
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

    fn sub_noborrow(repr: &mut [u64], value: u64) {
        let mut borrow = 0;

        repr[0] = crate::arithmetics::sbb(repr[0], value, &mut borrow);

        for a in repr.iter_mut().skip(1) {
            *a = crate::arithmetics::sbb(*a, 0u64, &mut borrow);
        }
    }

    fn add_nocarry(repr: &mut [u64], value: u64) {
        let mut carry = 0;

        repr[0] = crate::arithmetics::adc(repr[0], value, &mut carry);

        for a in repr.iter_mut().skip(1) {
            *a = crate::arithmetics::adc(*a, 0u64, &mut carry);
        }
    }

    if repr.len() == 0 {
        return vec![];
    }

    let mut res = vec![];
    let mut e = repr.to_vec();

    const WINDOW: u64 = 1u64;
    const MIDPOINT: u64 = 1u64 << WINDOW;
    const MIDPOINT_I64: i64 = MIDPOINT as i64;
    const MASK: u64 = 1u64 << (WINDOW + 1u64);

    while !is_zero(&e) {
        let z: i64;
        if is_odd(&e) {
            z = MIDPOINT_I64 - ((e[0] % MASK) as i64);
            if z >= 0 {
                sub_noborrow(&mut e, z as u64);
            } else {
                add_nocarry(&mut e, (-z) as u64);
            }
        } else {
            z = 0i64;
        }
        res.push(z);
        div2(&mut e);
    }

    res
}

#[cfg(test)]
mod test {
    use num_bigint::BigUint;
    use crate::field::{U256Repr, new_field};
    use crate::fp::Fp;
    use crate::traits::{FieldElement, ZeroAndOne};
    use crate::extension_towers::{Fp2Fp4FrobeniusBaseElements, Fp3Fp6FrobeniusBaseElements, Fp6Fp12FrobeniusBaseElements};
    use crate::extension_towers::fp2::{Fp2, Extension2};
    use crate::extension_towers::fp3::{Fp3, Extension3};
    use crate::extension_towers::fp4_as_2_over_2::{Extension2Over2};
    use crate::extension_towers::fp6_as_2_over_3::{Extension2Over3};
    use crate::extension_towers::fp6_as_3_over_2::{Fp6, Extension3Over2};
    use crate::extension_towers::fp12_as_2_over3_over_2::{Extension2Over3Over2};
    use num_traits::Num;
    use crate::test::{biguint_to_u64_vec};
    use crate::sliding_window_exp::WindowExpBase;
    use crate::constants::MaxFieldUint;

    #[test]
    fn test_fp12_faster_frobenius() {
        let modulus = BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let base_field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let modulus_biguint = modulus.clone();

        let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

        let mut extension_2 = Extension2::new(fp_non_residue);
        extension_2.calculate_frobenius_coeffs(&modulus).expect("must work");

        let one = Fp::one(&base_field);

        // non-residue is u+9
        let mut fp2_non_residue = Fp2::zero(&extension_2);
        let fp_9_repr = U256Repr::from(9u64);
        let fp_9 = Fp::from_repr(&base_field, fp_9_repr).unwrap(); 
        fp2_non_residue.c0 = fp_9.clone();
        fp2_non_residue.c1 = one.clone();

        let mut q_squared_minus_one = modulus_biguint.clone();
        q_squared_minus_one *= &modulus_biguint;
        q_squared_minus_one -= BigUint::from(1u64);

        let one_fp2 = Fp2::one(&extension_2);

        let may_be_one = fp2_non_residue.pow(&biguint_to_u64_vec(q_squared_minus_one.clone()));

        assert!(one_fp2 == may_be_one);

        // (q^6 - 1) = (q^2 - 1)*(q^4 + q^2 + 1)
 
        // largest power that we have to calculate for frobenius is 
        // (q^6 - 1)/6 = (q^2 - 1)(q^4 + q^2 + 1)/6
        // we know that 6 | q^2 - 1 (cause we use it for lower power)

        // (q^4 + q^2 + 1) = (q^2 - 1)*(q^2 + 2) + 3
        // remember that Fp2 ^ (q^2 - 1) == 1, so 
        // making a power of (q^6 - 1)/6 is equivalent to making a power of 
        // (q^2 - 1)/6 * 3

        let mut q_squared_minus_one_by_six = modulus_biguint.clone();
        q_squared_minus_one_by_six *= &modulus_biguint;
        q_squared_minus_one_by_six -= BigUint::from(1u64);
        q_squared_minus_one_by_six /= BigUint::from(6u64);

        let fp2_in_t = fp2_non_residue.pow(&biguint_to_u64_vec(q_squared_minus_one_by_six.clone()));

        let may_be_frob6 = fp2_in_t.pow(&[3u64]);

        let mut q_in_six_minus_one_by_six = modulus_biguint.clone();
        q_in_six_minus_one_by_six *= &modulus_biguint;
        q_in_six_minus_one_by_six *= &modulus_biguint;
        q_in_six_minus_one_by_six *= &modulus_biguint;
        q_in_six_minus_one_by_six *= &modulus_biguint;
        q_in_six_minus_one_by_six *= &modulus_biguint;
        q_in_six_minus_one_by_six -= BigUint::from(1u64);
        q_in_six_minus_one_by_six /= BigUint::from(6u64);

        let frob_6 = fp2_non_residue.pow(&biguint_to_u64_vec(q_in_six_minus_one_by_six.clone()));

        assert!(frob_6 == may_be_frob6);

        let exp_base = WindowExpBase::new(&fp2_non_residue, Fp2::one(&extension_2), 8, 7);

        let mut extension_6 = Extension3Over2::new(fp2_non_residue.clone());
        extension_6.calculate_frobenius_coeffs(&modulus, &exp_base).expect("must work");

        let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(&extension_6));
        extension_12.calculate_frobenius_coeffs(&modulus, &exp_base).expect("must work");

        assert_eq!(extension_12.frobenius_coeffs_c1[6], may_be_frob6);

    }

    #[test]
    fn test_fp12_fast_precomp() {
        let modulus = BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let base_field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

        let mut extension_2 = Extension2::new(fp_non_residue);
        extension_2.calculate_frobenius_coeffs(&modulus).expect("must work");

        let one = Fp::one(&base_field);

        // non-residue is u+9
        let mut fp2_non_residue = Fp2::zero(&extension_2);
        let fp_9_repr = U256Repr::from(9u64);
        let fp_9 = Fp::from_repr(&base_field, fp_9_repr).unwrap(); 
        fp2_non_residue.c0 = fp_9.clone();
        fp2_non_residue.c1 = one.clone();

        let exp_base = WindowExpBase::new(&fp2_non_residue, Fp2::one(&extension_2), 8, 7);

        let mut extension_6 = Extension3Over2::new(fp2_non_residue.clone());
        extension_6.calculate_frobenius_coeffs(&modulus, &exp_base).expect("must work");

        let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(&extension_6));
        extension_12.calculate_frobenius_coeffs(&modulus, &exp_base).expect("must work");

        let precomp_base = Fp6Fp12FrobeniusBaseElements::construct(&modulus, &fp2_non_residue).unwrap();

        let mut extension_6_fast = Extension3Over2::new(fp2_non_residue.clone());
        extension_6_fast.calculate_frobenius_coeffs_with_precomp(&precomp_base).expect("must work");

        let mut extension_12_fast = Extension2Over3Over2::new(Fp6::zero(&extension_6));
        extension_12_fast.calculate_frobenius_coeffs_with_precomp(&precomp_base).expect("must work");

        assert!(extension_6_fast.frobenius_coeffs_c1 == extension_6.frobenius_coeffs_c1);
        assert!(extension_6_fast.frobenius_coeffs_c2 == extension_6.frobenius_coeffs_c2);
        assert!(extension_12_fast.frobenius_coeffs_c1 == extension_12.frobenius_coeffs_c1);
    }

    #[test]
    fn test_fp4_fast_precomp() {
        use crate::field::U768Repr;

        let modulus = BigUint::from_str_radix("41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888253786114353726529584385201591605722013126468931404347949840543007986327743462853720628051692141265303114721689601", 10).unwrap();
        let base_field = new_field::<U768Repr>("41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888253786114353726529584385201591605722013126468931404347949840543007986327743462853720628051692141265303114721689601", 10).unwrap();
        let nonres_repr = U768Repr::from(13);
        let fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();
        
        let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

        let mut extension_2 = Extension2::new(fp_non_residue.clone());
        extension_2.calculate_frobenius_coeffs(&modulus).expect("must work");

        let fp2_non_residue = Fp2::zero(&extension_2); // non-residue is 13 + 0*u + 0*u^2
        // fp2_non_residue.c0 = fp_non_residue;

        let mut extension_4 = Extension2Over2::new(fp2_non_residue.clone());
        extension_4.calculate_frobenius_coeffs(&modulus).expect("must work");

        let precomp_base = Fp2Fp4FrobeniusBaseElements::construct(&modulus, &fp_non_residue).unwrap();

        let mut extension_2_fast = Extension2::new(fp_non_residue.clone());
        extension_2_fast.calculate_frobenius_coeffs_with_precomp(&precomp_base).expect("must work");

        let mut extension_4_fast = Extension2Over2::new(fp2_non_residue);
        extension_4_fast.calculate_frobenius_coeffs_with_precomp(&precomp_base).expect("must work");

        assert!(extension_2.frobenius_coeffs_c1 == extension_2_fast.frobenius_coeffs_c1);
        assert!(extension_4.frobenius_coeffs_c1 == extension_4_fast.frobenius_coeffs_c1);
    }

    #[test]
    fn test_fp6_fast_precomp() {
        use crate::field::{U640Repr};
        let modulus = BigUint::from_str_radix("19050022797317891600939264904924934656417895081121634056186244048763811669585984032184028629480644260294123843823582617865870693473572190965725707704312821545976965077621486794922414287", 10).unwrap();
        let base_field = new_field::<U640Repr>("19050022797317891600939264904924934656417895081121634056186244048763811669585984032184028629480644260294123843823582617865870693473572190965725707704312821545976965077621486794922414287", 10).unwrap();
        let nonres_repr = U640Repr::from(3);
        let fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

        let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

        let mut extension_3 = Extension3::new(fp_non_residue.clone());
        extension_3.calculate_frobenius_coeffs(&modulus).expect("must work");

        let mut fp3_non_residue = Fp3::zero(&extension_3);
        fp3_non_residue.c0 = fp_non_residue.clone();

        let mut extension_6 = Extension2Over3::new(fp3_non_residue.clone());
        extension_6.calculate_frobenius_coeffs(&modulus).expect("must work");

        let precomp_base = Fp3Fp6FrobeniusBaseElements::construct(&modulus, &fp_non_residue).unwrap();

        let mut extension_3_fast = Extension3::new(fp_non_residue.clone());
        extension_3_fast.calculate_frobenius_coeffs_with_precomp(&precomp_base).expect("must work");

        let mut extension_6_fast = Extension2Over3::new(fp3_non_residue);
        extension_6_fast.calculate_frobenius_coeffs_with_precomp(&precomp_base).expect("must work");

        assert!(extension_3.frobenius_coeffs_c1 == extension_3_fast.frobenius_coeffs_c1);
        assert!(extension_3.frobenius_coeffs_c2 == extension_3_fast.frobenius_coeffs_c2);

        for (i, (a, b)) in extension_6.frobenius_coeffs_c1.iter().zip(extension_6_fast.frobenius_coeffs_c1.iter()).enumerate() {
            assert_eq!(a, b, "Invalid coeff number {}", i);
        }
        // assert!(extension_6.frobenius_coeffs_c1 == extension_6_fast.frobenius_coeffs_c1);
    }


}