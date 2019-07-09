use crate::field::{SizedPrimeField};
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, BitIterator, FieldExtension, ZeroAndOne};
use super::fp2::{Fp2, Extension2};


// this implementation assumes extension using polynomial v^3 - xi = 0
// multiply_by_non_residue function in the extension field actually depends
// on the nature of the higher extension, but is supplied at runtime
// both BLS12-381 and BN254 uses `v`
pub struct Fp6<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >{
    pub c0: Fp2<'a, E, F>,
    pub c1: Fp2<'a, E, F>,
    pub c2: Fp2<'a, E, F>,
    pub extension_field: &'a Extension3Over2<'a, E, F>
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >std::fmt::Display for Fp6<'a, E, F> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fq6({} + {} * v + {} * v^2)", self.c0, self.c1, self.c2)
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >std::fmt::Debug for Fp6<'a, E, F> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fq6({} + {} * v + {} * v^2)", self.c0, self.c1, self.c2)
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Clone for Fp6<'a, E, F> {
    #[inline(always)]
    fn clone(&self) -> Self {
        Self{
            c0: self.c0.clone(),
            c1: self.c1.clone(),
            c2: self.c2.clone(),
            extension_field: self.extension_field
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > PartialEq for Fp6<'a, E, F> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && 
        self.c1 == other.c1 &&
        self.c2 == other.c2
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Eq for Fp6<'a, E, F> {
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Fp6<'a, E, F> {
    // pub fn zero(extension_field: &'a Extension3Over2<'a, E, F>) -> Self {
    //     let zero = Fp2::zero(extension_field.field);
        
    //     Self {
    //         c0: zero.clone(),
    //         c1: zero.clone(),
    //         c2: zero,
    //         extension_field: extension_field
    //     }
    // }

    // pub fn one(extension_field: &'a Extension3Over2<'a, E, F>) -> Self {
    //     let zero = Fp2::zero(extension_field.field);
    //     let one = Fp2::one(extension_field.field);
        
    //     Self {
    //         c0: one,
    //         c1: zero.clone(),
    //         c2: zero,
    //         extension_field: extension_field
    //     }
    // }

    pub fn mul_by_1(&mut self, c1: &Fp2<'a, E, F>) {
        let mut b_b = self.c1.clone();
        b_b.mul_assign(c1);

        let mut t1 = c1.clone();
        {
            let mut tmp = self.c1.clone();
            tmp.add_assign(&self.c2);

            t1.mul_assign(&tmp);
            t1.sub_assign(&b_b);
            t1.mul_by_nonresidue(self.extension_field);
        }

        let mut t2 = c1.clone();
        {
            let mut tmp = self.c0.clone();
            tmp.add_assign(&self.c1);

            t2.mul_assign(&tmp);
            t2.sub_assign(&b_b);
        }

        self.c0 = t1;
        self.c1 = t2;
        self.c2 = b_b;
    }

    pub fn mul_by_01(&mut self, c0: &Fp2<'a, E, F>, c1: &Fp2<'a, E, F>) {
        let mut a_a = self.c0.clone();
        let mut b_b = self.c1.clone();
        a_a.mul_assign(c0);
        b_b.mul_assign(c1);

        let mut t1 = c1.clone();
        {
            let mut tmp = self.c1.clone();
            tmp.add_assign(&self.c2);

            t1.mul_assign(&tmp);
            t1.sub_assign(&b_b);
            t1.mul_by_nonresidue(self.extension_field);
            t1.add_assign(&a_a);
        }

        let mut t3 = c0.clone();
        {
            let mut tmp = self.c0.clone();
            tmp.add_assign(&self.c2);

            t3.mul_assign(&tmp);
            t3.sub_assign(&a_a);
            t3.add_assign(&b_b);
        }

        let mut t2 = c0.clone();
        t2.add_assign(c1);
        {
            let mut tmp = self.c0.clone();
            tmp.add_assign(&self.c1);

            t2.mul_assign(&tmp);
            t2.sub_assign(&a_a);
            t2.sub_assign(&b_b);
        }

        self.c0 = t1;
        self.c1 = t2;
        self.c2 = t3;
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > ZeroAndOne for Fp6<'a, E, F> {
    type Params = &'a Extension3Over2<'a, E, F>;

    fn zero(extension_field: &'a Extension3Over2<'a, E, F>) -> Self {
        let zero = Fp2::zero(extension_field.field);
        
        Self {
            c0: zero.clone(),
            c1: zero.clone(),
            c2: zero,
            extension_field: extension_field
        }
    }

    fn one(extension_field: &'a Extension3Over2<'a, E, F>) -> Self {
        let zero = Fp2::zero(extension_field.field);
        let one = Fp2::one(extension_field.field);
        
        Self {
            c0: one,
            c1: zero.clone(),
            c2: zero,
            extension_field: extension_field
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldElement for Fp6<'a, E, F> {
    /// Returns true iff this element is zero.
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && 
        self.c1.is_zero() &&
        self.c2.is_zero()
    }

    fn add_assign(&mut self, other: &Self) {
        self.c0.add_assign(&other.c0);
        self.c1.add_assign(&other.c1);
        self.c2.add_assign(&other.c2);
    }

    fn double(&mut self) {
        self.c0.double();
        self.c1.double();
        self.c2.double();
    }

    fn sub_assign(&mut self, other: &Self) {
        self.c0.sub_assign(&other.c0);
        self.c1.sub_assign(&other.c1);
        self.c2.sub_assign(&other.c2);
    }

    fn negate(&mut self) {
        self.c0.negate();
        self.c1.negate();
        self.c2.negate();
    }

    fn inverse(&self) -> Option<Self> {
        let mut c0 = self.c2.clone();
        c0.mul_by_nonresidue(self.extension_field);
        c0.mul_assign(&self.c1);
        c0.negate();
        {
            let mut c0s = self.c0.clone();
            c0s.square();
            c0.add_assign(&c0s);
        }
        let mut c1 = self.c2.clone();
        c1.square();
        c1.mul_by_nonresidue(self.extension_field);
        {
            let mut c01 = self.c0.clone();
            c01.mul_assign(&self.c1);
            c1.sub_assign(&c01);
        }
        let mut c2 = self.c1.clone();
        c2.square();
        {
            let mut c02 = self.c0.clone();
            c02.mul_assign(&self.c2);
            c2.sub_assign(&c02);
        }

        let mut tmp1 = self.c2.clone();
        tmp1.mul_assign(&c1);
        let mut tmp2 = self.c1.clone();
        tmp2.mul_assign(&c2);
        tmp1.add_assign(&tmp2);
        tmp1.mul_by_nonresidue(self.extension_field);
        tmp2 = self.c0.clone();
        tmp2.mul_assign(&c0);
        tmp1.add_assign(&tmp2);

        match tmp1.inverse() {
            Some(t) => {
                let mut tmp = Fp6 {
                    c0: t.clone(),
                    c1: t.clone(),
                    c2: t,
                    extension_field: self.extension_field
                };
                tmp.c0.mul_assign(&c0);
                tmp.c1.mul_assign(&c1);
                tmp.c2.mul_assign(&c2);

                Some(tmp)
            }
            None => None,
        }
    }

    fn mul_assign(&mut self, other: &Self)
    {
        let mut a_a = self.c0.clone();
        let mut b_b = self.c1.clone();
        let mut c_c = self.c2.clone();
        a_a.mul_assign(&other.c0);
        b_b.mul_assign(&other.c1);
        c_c.mul_assign(&other.c2);

        let mut t1 = other.c1.clone();
        t1.add_assign(&other.c2);
        {
            let mut tmp = self.c1.clone();
            tmp.add_assign(&self.c2);

            t1.mul_assign(&tmp);
            t1.sub_assign(&b_b);
            t1.sub_assign(&c_c);
            t1.mul_by_nonresidue(self.extension_field);
            t1.add_assign(&a_a);
        }

        let mut t3 = other.c0.clone();
        t3.add_assign(&other.c2);
        {
            let mut tmp = self.c0.clone();
            tmp.add_assign(&self.c2);

            t3.mul_assign(&tmp);
            t3.sub_assign(&a_a);
            t3.add_assign(&b_b);
            t3.sub_assign(&c_c);
        }

        let mut t2 = other.c0.clone();
        t2.add_assign(&other.c1);
        {
            let mut tmp = self.c0.clone();
            tmp.add_assign(&self.c1);

            t2.mul_assign(&tmp);
            t2.sub_assign(&a_a);
            t2.sub_assign(&b_b);
            c_c.mul_by_nonresidue(self.extension_field);
            t2.add_assign(&c_c);
        }

        self.c0 = t1;
        self.c1 = t2;
        self.c2 = t3;
    }

    fn square(&mut self)
    {
        let mut s0 = self.c0.clone();
        s0.square();
        let mut ab = self.c0.clone();
        ab.mul_assign(&self.c1);
        let mut s1 = ab;
        s1.double();
        let mut s2 = self.c0.clone();
        s2.sub_assign(&self.c1);
        s2.add_assign(&self.c2);
        s2.square();
        let mut bc = self.c1.clone();
        bc.mul_assign(&self.c2);
        let mut s3 = bc.clone();
        s3.double();
        let mut s4 = self.c2.clone();
        s4.square();

        self.c0 = s3.clone();
        self.c0.mul_by_nonresidue(self.extension_field);
        self.c0.add_assign(&s0);

        self.c1 = s4.clone();
        self.c1.mul_by_nonresidue(self.extension_field);
        self.c1.add_assign(&s1);

        self.c2 = s1;
        self.c2.add_assign(&s2);
        self.c2.add_assign(&s3);
        self.c2.sub_assign(&s0);
        self.c2.sub_assign(&s4);
    }

    fn conjugate(&mut self) {
        unreachable!();
    }

    fn pow<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Self::one(&self.extension_field);

        let mut found_one = false;

        for i in BitIterator::new(exp) {
            if found_one {
                res.square();
            } else {
                found_one = i;
            }

            if i {
                res.mul_assign(self);
            }
        }

        res
    }

    fn mul_by_nonresidue<EXT: FieldExtension<Element = Self>>(&mut self, for_extesion: &EXT) {
        for_extesion.multiply_by_non_residue(self);
        // self.extension_field.multiply_by_non_residue(self);
    }

    fn frobenius_map(&mut self, power: usize) {
        debug_assert!(self.extension_field.frobenius_coeffs_are_calculated);
        match power {
            0 | 1 | 2 | 3 | 6 => {

            },
            _ => {
                unreachable!("can not reach power {}", power);
            }
        }
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);
        self.c2.frobenius_map(power);

        self.c1.mul_assign(&self.extension_field.frobenius_coeffs_c1[power % 6]);
        self.c2.mul_assign(&self.extension_field.frobenius_coeffs_c2[power % 6]);
    }
}

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::Zero;
use crate::sliding_window_exp::{WindowExpBase};

// For example, BLS12-381 has non-residue = 1 + u;
pub struct Extension3Over2<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > {
    pub(crate) non_residue: Fp2<'a, E, F>,
    pub(crate) field: &'a Extension2<'a, E, F>,
    pub(crate) frobenius_coeffs_c1: [Fp2<'a, E, F>; 6],
    pub(crate) frobenius_coeffs_c2: [Fp2<'a, E, F>; 6],
    pub(crate) frobenius_coeffs_are_calculated: bool
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Extension3Over2<'a, E, F> {
    pub (crate) fn new(non_residue: Fp2<'a, E, F>) -> Self {
        let extension_2 = &non_residue.extension_field;

        let zeros = [Fp2::zero(extension_2), Fp2::zero(extension_2), Fp2::zero(extension_2),
                    Fp2::zero(extension_2), Fp2::zero(extension_2), Fp2::zero(extension_2)];
        
        Self {
            non_residue: non_residue.clone(),
            field: & non_residue.extension_field,
            frobenius_coeffs_c1: zeros.clone(),
            frobenius_coeffs_c2: zeros,
            frobenius_coeffs_are_calculated: false
        }
    }

    pub(crate) fn calculate_frobenius_coeffs(
        &mut self,
        modulus: BigUint,
        base: &WindowExpBase<Fp2<'a, E, F>>
    ) -> Result<(), ()> {
        use crate::field::biguint_to_u64_vec;
        use crate::constants::ONE_BIGUINT;
        use crate::constants::THREE_BIGUINT;

        // let one = BigUint::from(1u64);
        // let three = BigUint::from(3u64);

        // NON_RESIDUE**(((q^0) - 1) / 3)
        // let non_residue = extension.non_residue.clone();
        let f_0 = Fp2::one(self.field);

        let mut powers = vec![];

        let mut q_power = modulus.clone();

        {
            let power = q_power.clone() - &*ONE_BIGUINT;
            let (power, rem) = power.div_rem(&*THREE_BIGUINT);
            if !rem.is_zero() {
                return Err(());
            }
            // debug_assert!(rem.is_zero());
            powers.push(biguint_to_u64_vec(power));
        }
        for _ in 1..3 {
            q_power *= &modulus;
            let power = q_power.clone() - &*ONE_BIGUINT;
            let (power, rem) = power.div_rem(&*THREE_BIGUINT);
            if !rem.is_zero() {
                return Err(());
            }
            // debug_assert!(rem.is_zero());
            powers.push(biguint_to_u64_vec(power));
        }

        let mut result = base.exponentiate(&powers);
        debug_assert!(result.len() == 3);

        let f_3 = result.pop().expect("has enough elements");
        let f_2 = result.pop().expect("has enough elements");
        let f_1 = result.pop().expect("has enough elements");

        let f_4 = Fp2::zero(self.field);
        let f_5 = Fp2::zero(self.field);

        let f_0_c2 = f_0.clone();

        let mut f_1_c2 = f_1.clone();
        f_1_c2.square();
        let mut f_2_c2 = f_2.clone();
        f_2_c2.square();
        let mut f_3_c2 = f_3.clone();
        f_3_c2.square();

        let f_4_c2 = f_4.clone();
        let f_5_c2 = f_5.clone();

        self.frobenius_coeffs_c1 = [f_0, f_1, f_2, f_3, f_4, f_5];
        self.frobenius_coeffs_c2 = [f_0_c2, f_1_c2, f_2_c2, f_3_c2, f_4_c2, f_5_c2];
        self.frobenius_coeffs_are_calculated = true;

        Ok(())
    }

}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldExtension for Extension3Over2<'a, E, F> {
    const EXTENSION_DEGREE: usize = 3;
    
    type Element = Fp2<'a, E, F>;

    fn multiply_by_non_residue(&self, el: &mut Self::Element) {
        // this is simply a multiplication by non-residue that is Fp element cause everything else 
        // is covered in explicit formulas for multiplications for Fp2
        el.mul_assign(&self.non_residue);
    }

}
