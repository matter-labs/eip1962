use crate::fp::Fp;
use crate::field::{SizedPrimeField};
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, BitIterator, FieldExtension, ZeroAndOne};
use super::fp3::{Fp3, Extension3};

pub struct Fp6<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >{
    pub c0: Fp3<'a, E, F>,
    pub c1: Fp3<'a, E, F>,
    pub extension_field: &'a Extension2Over3<'a, E, F>
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >std::fmt::Display for Fp6<'a, E, F> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Fq6({} + {} * v)", self.c0, self.c1)
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >std::fmt::Debug for Fp6<'a, E, F> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Fq6({} + {} * v)", self.c0, self.c1)
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Clone for Fp6<'a, E, F> {
    #[inline(always)]
    fn clone(&self) -> Self {
        Self{
            c0: self.c0.clone(),
            c1: self.c1.clone(),
            extension_field: self.extension_field
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > PartialEq for Fp6<'a, E, F> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && 
        self.c1 == other.c1
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Eq for Fp6<'a, E, F> {
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Fp6<'a, E, F> {
    pub fn zero(extension_field: &'a Extension2Over3<'a, E, F>) -> Self {
        let zero = Fp3::zero(extension_field.field);
        
        Self {
            c0: zero.clone(),
            c1: zero,
            extension_field: extension_field
        }
    }

    pub fn one(extension_field: &'a Extension2Over3<'a, E, F>) -> Self {
        let zero = Fp3::zero(extension_field.field);
        let one = Fp3::one(extension_field.field);
        
        Self {
            c0: one,
            c1: zero,
            extension_field: extension_field
        }
    }

    pub fn cyclotomic_exp<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Self::one(self.extension_field);
        let mut self_inverse = self.clone();
        self_inverse.conjugate();

        let mut found_nonzero = false;
        use crate::pairings::into_ternary_wnaf;
        let naf = into_ternary_wnaf(exp.as_ref());

        for &value in naf.iter().rev() {
            if found_nonzero {
                res.square();
            }

            if value != 0 {
                found_nonzero = true;

                if value > 0 {
                    res.mul_assign(&self);
                } else {
                    res.mul_assign(&self_inverse);
                }
            }
        }

        res
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > ZeroAndOne for Fp6<'a, E, F> {
    type Params = &'a Extension2Over3<'a, E, F>;

    fn zero(extension_field: &'a Extension2Over3<'a, E, F>) -> Self {
        let zero = Fp3::zero(extension_field.field);
        
        Self {
            c0: zero.clone(),
            c1: zero,
            extension_field: extension_field
        }
    }

    fn one(extension_field: &'a Extension2Over3<'a, E, F>) -> Self {
        let zero = Fp3::zero(extension_field.field);
        let one = Fp3::one(extension_field.field);
        
        Self {
            c0: one,
            c1: zero,
            extension_field: extension_field
        }
    }

}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldElement for Fp6<'a, E, F> {
    /// Returns true iff this element is zero.
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && 
        self.c1.is_zero()
    }

    fn add_assign(&mut self, other: &Self) {
        self.c0.add_assign(&other.c0);
        self.c1.add_assign(&other.c1);
    }

    fn double(&mut self) {
        self.c0.double();
        self.c1.double();
    }

    fn sub_assign(&mut self, other: &Self) {
        self.c0.sub_assign(&other.c0);
        self.c1.sub_assign(&other.c1);
    }

    fn negate(&mut self) {
        self.c0.negate();
        self.c1.negate();
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // From "High-Speed Software Implementation of the Optimal Ate Pairing over
            // Barreto-Naehrig
            // Curves"; Algorithm 8
            let a = self.c0.clone();
            let b = self.c1.clone();

            let mut t1 = b.clone();
            t1.square();
            let mut t0 = a.clone();
            t0.square();

            let mut v0 = t1.clone();
            v0.mul_by_nonresidue(self.extension_field);
            t0.sub_assign(&v0);

            let t2 = t0.inverse();
            if t2.is_none() {
                return None;
            }
            
            let t2 = t2.expect("is not None");

            let mut c0 = a;
            c0.mul_assign(&t2);
            let mut c1 = b;
            c1.mul_assign(&t2);
            c1.negate();

            Some(Self {
                c0, 
                c1,
                extension_field: self.extension_field
            })
        }
    }

    fn mul_assign(&mut self, other: &Self)
    {
        let a0 = self.c0.clone();
        let b0 = self.c1.clone();
        let a1 = other.c0.clone();
        let b1 = other.c1.clone();

        let mut a0a1 = a0.clone();
        a0a1.mul_assign(&a1);
        let mut b0b1 = b0.clone();
        b0b1.mul_assign(&b1);
        let mut t0 = b0b1.clone();
        t0.mul_by_nonresidue(self.extension_field);

        let mut c0 = a0a1.clone();
        c0.add_assign(&t0);
        let mut c1 = a0;
        c1.add_assign(&b0);

        let mut t1 = a1;
        t1.add_assign(&b1);

        c1.mul_assign(&t1);
        c1.sub_assign(&a0a1);
        c1.sub_assign(&b0b1);

        self.c0 = c0;
        self.c1 = c1;
    }

    fn square(&mut self)
    {
        let a = self.c0.clone();
        let b = self.c1.clone();
        let mut ab_add = a.clone();
        ab_add.add_assign(&b);
        let mut ab_mul = a.clone();
        ab_mul.mul_assign(&b);

        let mut t0 = b.clone();
        t0.mul_by_nonresidue(self.extension_field);
        t0.add_assign(&a);

        let mut t1 = ab_mul.clone();
        t1.mul_by_nonresidue(self.extension_field);

        let mut c0 = ab_add;
        c0.mul_assign(&t0);
        c0.sub_assign(&ab_mul);
        c0.sub_assign(&t1);
        
        let mut c1 = ab_mul;
        c1.double();

        self.c0 = c0;
        self.c1 = c1;
    }

    fn conjugate(&mut self) {
        self.c1.negate();
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
    }

    fn frobenius_map(&mut self, power: usize) {
        debug_assert!(self.extension_field.frobenius_coeffs_are_calculated);
        match power {
            1 | 3 => {

            },
            _ => {
                unreachable!("can not reach power {}", power);
            }
        }
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);
        self.c1.mul_by_fp(&self.extension_field.frobenius_coeffs_c1[power % 6]);
    }
}

pub struct Extension2Over3<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > {
    pub(crate) field: &'a Extension3<'a, E, F>,
    pub(crate) non_residue: Fp3<'a, E, F>,
    pub(crate) frobenius_coeffs_c1: [Fp<'a, E, F>; 6],
    pub(crate) frobenius_coeffs_are_calculated: bool
}

use crate::constants::*;

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Extension2Over3<'a, E, F> {
    pub (crate) fn new(non_residue: Fp3<'a, E, F>) -> Self {
        let field = non_residue.extension_field.field;

        let zeros = [Fp::zero(field), Fp::zero(field), Fp::zero(field),
                    Fp::zero(field), Fp::zero(field), Fp::zero(field)];
        
        Self {
            non_residue: non_residue.clone(),
            field: & non_residue.extension_field,
            frobenius_coeffs_c1: zeros,
            frobenius_coeffs_are_calculated: false
        }
    }

    pub(crate) fn calculate_frobenius_coeffs(
        &mut self,
        modulus: &MaxFieldUint,
    ) -> Result<(), ()> {
    
        // NON_REDISUE**(((q^0) - 1) / 6)
        let non_residue = self.field.non_residue.clone();
        let f_0 = Fp::one(self.field.field);

        let modulus = MaxFrobeniusFp6::from(modulus.as_ref());
        let mut q_power = modulus;
        let one = MaxFrobeniusFp6::from(1u64);
        let six = MaxFrobeniusFp6::from(6u64);

        let power = q_power - one;
        let (power, rem) = power.div_mod(six);
        if !rem.is_zero() {
            if !std::option_env!("GAS_METERING").is_some() {
                return Err(());
            }
        }
        let f_1 = non_residue.pow(power.as_ref());

        q_power *= modulus;
        let f_2 = Fp::zero(self.field.field);

        q_power *= modulus;
        let power = q_power - one;
        let (power, rem) = power.div_mod(six);
        if !rem.is_zero() {
            if !std::option_env!("GAS_METERING").is_some() {
                return Err(());
            }
        }
        let f_3 = non_residue.pow(power.as_ref());

        let f_4 = Fp::zero(self.field.field);
        let f_5 = Fp::zero(self.field.field);

        self.frobenius_coeffs_c1 = [f_0, f_1, f_2, f_3, f_4, f_5];
        self.frobenius_coeffs_are_calculated = true;

        Ok(())
    }

}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldExtension for Extension2Over3<'a, E, F> {
    const EXTENSION_DEGREE: usize = 2;
    
    type Element = Fp3<'a, E, F>;

    fn multiply_by_non_residue(&self, el: &mut Self::Element) {
        // IMPORTANT: This only works cause the structure of extension field for Fp6
        // is w^2 - u = 0!
        // take an element in Fp6 as 2 over 3 and mutplity
        // (c0 + c1 * u)*u with u^2 - xi = 0 -> (c1*xi + c0 * u)
        let mut c0 = el.c2.clone();
        el.c2 = el.c1.clone();
        el.c1 = el.c0.clone();
        c0.mul_by_nonresidue(&*el.extension_field);
        el.c0 = c0;
    }
}
