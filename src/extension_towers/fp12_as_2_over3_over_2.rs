use crate::field::{SizedPrimeField};
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, BitIterator, FieldExtension, ZeroAndOne};
use super::fp6_as_3_over_2::{Fp6, Extension3Over2};
use super::fp2::Fp2;
use super::Fp6Fp12FrobeniusBaseElements;

// this implementation assumes extension using polynomial w^2 - v = 0
pub struct Fp12<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >{
    pub c0: Fp6<'a, E, F>,
    pub c1: Fp6<'a, E, F>,
    pub extension_field: &'a Extension2Over3Over2<'a, E, F>
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >std::fmt::Display for Fp12<'a, E, F> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fq12({} + {} * w)", self.c0, self.c1)
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >std::fmt::Debug for Fp12<'a, E, F> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fq12({} + {} * w)", self.c0, self.c1)
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Clone for Fp12<'a, E, F> {
    #[inline(always)]
    fn clone(&self) -> Self {
        Self{
            c0: self.c0.clone(),
            c1: self.c1.clone(),
            extension_field: self.extension_field
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > PartialEq for Fp12<'a, E, F> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && 
        self.c1 == other.c1
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Eq for Fp12<'a, E, F> {
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Fp12<'a, E, F> {
    pub fn mul_by_034(
        &mut self,
        c0: & Fp2<'a, E, F>,
        c3: & Fp2<'a, E, F>,
        c4: & Fp2<'a, E, F>,
    ) {
        let mut a = self.c0.clone();
        a.c0.mul_assign(c0);
        a.c1.mul_assign(c0);
        a.c2.mul_assign(c0);

        let mut b = self.c1.clone();
        b.mul_by_01(&c3, &c4);

        let mut t0 = c0.clone();
        t0.add_assign(c3);

        let mut e = self.c0.clone();
        e.add_assign(&self.c1);
        e.mul_by_01(&t0, &c4);

        self.c1 = e;
        self.c1.sub_assign(&a);
        self.c1.sub_assign(&b);


        let mut t1 = b.clone();
        t1.mul_by_nonresidue(self.extension_field);
        self.c0 = a;
        self.c0.add_assign(&t1);
    }

    pub fn mul_by_014(
        &mut self,
        c0: & Fp2<'a, E, F>,
        c1: & Fp2<'a, E, F>,
        c4: & Fp2<'a, E, F>,
    ) {
        let mut aa = self.c0.clone();
        aa.mul_by_01(c0, c1);
        let mut bb = self.c1.clone();
        bb.mul_by_1(c4);
        let mut o = c1.clone();
        o.add_assign(c4);
        self.c1.add_assign(&self.c0);
        self.c1.mul_by_01(c0, &o);
        self.c1.sub_assign(&aa);
        self.c1.sub_assign(&bb);
        self.c0 = bb;
        self.c0.mul_by_nonresidue(self.extension_field);
        self.c0.add_assign(&aa);
    }

    pub fn cyclotomic_square(&mut self) {
        let z0 = self.c0.c0.clone();
        let z4 = self.c0.c1.clone();
        let z3 = self.c0.c2.clone();
        let z2 = self.c1.c0.clone();
        let z1 = self.c1.c1.clone();
        let z5 = self.c1.c2.clone();

        // t0 + t1*y = (z0 + z1*y)^2 = a^2
        let mut tmp = z0.clone();
        tmp.mul_assign(&z1);

        let mut a0 = z0.clone();
        a0.add_assign(&z1);
        let mut a1 = z1.clone();
        a1.mul_by_nonresidue(self.extension_field.field);
        a1.add_assign(&z0);

        let mut a2 = tmp.clone();
        a2.mul_by_nonresidue(self.extension_field.field);

        let mut t0 = a0;
        t0.mul_assign(&a1);
        t0.sub_assign(&tmp);
        t0.sub_assign(&a2);
        let mut t1 = tmp;
        t1.double();

        // t2 + t3*y = (z2 + z3*y)^2 = b^2
        let mut tmp = z2.clone();
        tmp.mul_assign(&z3);

        let mut a0 = z2.clone();
        a0.add_assign(&z3);
        let mut a1 = z3.clone();
        a1.mul_by_nonresidue(self.extension_field.field);
        a1.add_assign(&z2);

        let mut a2 = tmp.clone();
        a2.mul_by_nonresidue(self.extension_field.field);

        let mut t2 = a0;
        t2.mul_assign(&a1);
        t2.sub_assign(&tmp);
        t2.sub_assign(&a2);

        let mut t3 = tmp;
        t3.double();

        // t4 + t5*y = (z4 + z5*y)^2 = c^2
        let mut tmp = z4.clone();
        tmp.mul_assign(&z5);

        let mut a0 = z4.clone();
        a0.add_assign(&z5);
        let mut a1 = z5.clone();
        a1.mul_by_nonresidue(self.extension_field.field);
        a1.add_assign(&z4);

        let mut a2 = tmp.clone();
        a2.mul_by_nonresidue(self.extension_field.field);

        let mut t4 = a0;
        t4.mul_assign(&a1);
        t4.sub_assign(&tmp);
        t4.sub_assign(&a2);

        let mut t5 = tmp.clone();
        t5.double();

        // for A

        // g0 = 3 * t0 - 2 * z0
        let mut g0 = t0.clone();
        g0.sub_assign(&z0);
        g0.double();
        g0.add_assign(&t0);

        self.c0.c0 = g0;

        // g1 = 3 * t1 + 2 * z1
        let mut g1 = t1.clone();
        g1.add_assign(&z1);
        g1.double();
        g1.add_assign(&t1);
        self.c1.c1 = g1;

        // for B

        // g2 = 3 * (xi * t5) + 2 * z2
        let mut tmp = t5.clone();
        tmp.mul_by_nonresidue(self.extension_field.field);
        let mut g2 = tmp.clone();
        g2.add_assign(&z2);
        g2.double();
        g2.add_assign(&tmp);
        self.c1.c0 = g2;

        // g3 = 3 * t4 - 2 * z3
        let mut g3 = t4.clone();
        g3.sub_assign(&z3);
        g3.double();
        g3.add_assign(&t4);
        self.c0.c2 = g3;

        // for C

        // g4 = 3 * t2 - 2 * z4
        let mut g4 = t2.clone();
        g4.sub_assign(&z4);
        g4.double();
        g4.add_assign(&t2);
        self.c0.c1 = g4;

        // g5 = 3 * t3 + 2 * z5
        let mut g5 = t3.clone();
        g5.add_assign(&z5);
        g5.double();
        g5.add_assign(&t3);
        self.c1.c2 = g5;
    }

    pub fn cyclotomic_exp<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Self::one(&self.extension_field);

        let mut found_one = false;

        for i in BitIterator::new(exp) {
            if found_one {
                res.cyclotomic_square();
            } else {
                found_one = i;
            }

            if i {
                res.mul_assign(self);
            }
        }

        res
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > ZeroAndOne for Fp12<'a, E, F> {
    type Params = &'a Extension2Over3Over2<'a, E, F>;

    fn zero(extension_field: &'a Extension2Over3Over2<'a, E, F>) -> Self {
        let zero = Fp6::zero(extension_field.field);
        
        Self {
            c0: zero.clone(),
            c1: zero,
            extension_field: extension_field
        }
    }

    fn one(extension_field: &'a Extension2Over3Over2<'a, E, F>) -> Self {
        let zero = Fp6::zero(extension_field.field);
        let one = Fp6::one(extension_field.field);
        
        Self {
            c0: one,
            c1: zero,
            extension_field: extension_field
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldElement for Fp12<'a, E, F> {
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
        let mut c0s = self.c0.clone();
        c0s.square();
        let mut c1s = self.c1.clone();
        c1s.square();
        c1s.mul_by_nonresidue(self.extension_field);
        c0s.sub_assign(&c1s);

        c0s.inverse().map(|t| {
            let mut tmp = Fp12 { 
                c0: t.clone(), 
                c1: t.clone(),
                extension_field: self.extension_field
            };
            tmp.c0.mul_assign(&self.c0);
            tmp.c1.mul_assign(&self.c1);
            tmp.c1.negate();

            tmp
        })
    }

    fn mul_assign(&mut self, other: &Self)
    {
        let mut aa = self.c0.clone();
        aa.mul_assign(&other.c0);
        let mut bb = self.c1.clone();
        bb.mul_assign(&other.c1);
        let mut o = other.c0.clone();
        o.add_assign(&other.c1);
        self.c1.add_assign(&self.c0);
        self.c1.mul_assign(&o);
        self.c1.sub_assign(&aa);
        self.c1.sub_assign(&bb);
        self.c0 = bb;
        self.c0.mul_by_nonresidue(self.extension_field);
        self.c0.add_assign(&aa);
    }

    fn square(&mut self)
    {
        let mut ab = self.c0.clone();
        ab.mul_assign(&self.c1);
        let mut c0c1 = self.c0.clone();
        c0c1.add_assign(&self.c1);
        let mut c0 = self.c1.clone();
        c0.mul_by_nonresidue(self.extension_field);
        c0.add_assign(&self.c0);
        c0.mul_assign(&c0c1);
        c0.sub_assign(&ab);
        self.c1 = ab.clone();
        self.c1.add_assign(&ab);
        ab.mul_by_nonresidue(self.extension_field);
        c0.sub_assign(&ab);
        self.c0 = c0;
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

    fn mul_by_nonresidue<EXT: FieldExtension<Element = Self>>(&mut self, _for_extesion: &EXT) {
        unreachable!();
        // for_extesion.multiply_by_non_residue(self);
    }

    fn frobenius_map(&mut self, power: usize) {
        assert!(self.extension_field.frobenius_coeffs_are_calculated);
        match power {
            1 | 2 | 3 | 6 => {

            },
            _ => {
                unreachable!("can not reach power {}", power);
            }
        }
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);

        self.c1
            .c0
            .mul_assign(&self.extension_field.frobenius_coeffs_c1[power % 12]);
        self.c1
            .c1
            .mul_assign(&self.extension_field.frobenius_coeffs_c1[power % 12]);
        self.c1
            .c2
            .mul_assign(&self.extension_field.frobenius_coeffs_c1[power % 12]);
    }
}

pub struct Extension2Over3Over2<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > {
    pub(crate) non_residue: Fp6<'a, E, F>,
    pub(crate) field: &'a Extension3Over2<'a, E, F>,
    pub(crate) frobenius_coeffs_c1: [Fp2<'a, E, F>; 12],
    pub(crate) frobenius_coeffs_are_calculated: bool
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Clone for Extension2Over3Over2<'a, E, F> {
    fn clone(&self) -> Self {
        Self {
            non_residue: self.non_residue.clone(),
            field: self.field,
            frobenius_coeffs_c1: self.frobenius_coeffs_c1.clone(),
            frobenius_coeffs_are_calculated: self.frobenius_coeffs_are_calculated
        }
    }
}

use crate::integers::*;

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Extension2Over3Over2<'a, E, F> {
    pub (crate) fn new(non_residue: Fp6<'a, E, F>) -> Self {
        let extension_2 = &non_residue.extension_field.field;

        Self {
            non_residue: non_residue.clone(),
            field: & non_residue.extension_field,
            frobenius_coeffs_c1: [Fp2::zero(extension_2), Fp2::zero(extension_2), Fp2::zero(extension_2),
                                Fp2::zero(extension_2), Fp2::zero(extension_2), Fp2::zero(extension_2),
                                Fp2::zero(extension_2), Fp2::zero(extension_2), Fp2::zero(extension_2),
                                Fp2::zero(extension_2), Fp2::zero(extension_2), Fp2::zero(extension_2)],
            frobenius_coeffs_are_calculated: false
        }
    }

    #[allow(dead_code)]
    pub(crate) fn calculate_frobenius_coeffs_optimized(
        &mut self,
        modulus: &MaxFieldUint,
    ) -> Result<(), ()> {    
        // NON_RES**(((q^0) - 1) / 6)
        let f_0 = Fp2::one(self.field.field);

        let non_residue = &self.field.non_residue;

        // use a fact that Fp2 ** (q^2 - 1) == 1 and that 6 | q^2 - 1

        // then
        // c1 = Fp2**( (q^1 - 1) / 6) has to be calculated
        // c2 = Fp2**( (q^2 - 1) / 6) has to be calculated
        // c3 = Fp2**( (q^3 - 1) / 6) = Fp2**(( (q^2 - 1) / 6) * q + (q - 1) / 6 ) =
        // = Fp2**( ( (q^2 - 1) / 6) * q ) * Fp2**((q - 1) / 6) =
        // = c2.frobenius(1) * c1
        // c4 is not calculated
        // c5 is not calculated
        // c6 = Fp2**( (q^6 - 1) / 6) === Fp2**( (q^2 - 1)/6 * 3) = c2 ** 3

        let modulus = MaxFieldSquaredUint::from(modulus.as_ref());
        let mut q_power = modulus;
        let one = MaxFieldSquaredUint::from(1u64);
        let six = MaxFieldSquaredUint::from(6u64);

        // 1
        let f_1 = {
            let power = q_power - one;
            let (power, rem) = power.div_mod(six);
            if !rem.is_zero() {
                if !crate::features::in_gas_metering() {
                    return Err(());
                }
            }

            non_residue.pow(power.as_ref())
        };

        // 2
        let f_2 = {
            q_power = q_power.adaptive_multiplication(modulus);
            let power = q_power - one;
            let (power, rem) = power.div_mod(six);
            if !rem.is_zero() {
                if !crate::features::in_gas_metering() {
                    return Err(());
                }
            }

            non_residue.pow(power.as_ref())
        };

        let mut f_3 = f_2.clone();
        f_3.frobenius_map(1);
        f_3.mul_assign(&f_1);

        let f_6 = f_2.pow(&[3u64]);

        let f_4 = Fp2::zero(self.field.field);
        let f_5 = Fp2::zero(self.field.field);

        let f_7 = Fp2::zero(self.field.field);
        let f_8 = Fp2::zero(self.field.field);
        let f_9 = Fp2::zero(self.field.field);
        let f_10 = Fp2::zero(self.field.field);
        let f_11 = Fp2::zero(self.field.field);

        self.frobenius_coeffs_c1 = [f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11];
        self.frobenius_coeffs_are_calculated = true;

        Ok(())
    }

    pub(crate) fn calculate_frobenius_coeffs_with_precomp(
        &mut self,
        precomp: &Fp6Fp12FrobeniusBaseElements<'a, E, F>
    ) -> Result<(), ()> {    
        let f_0 = Fp2::one(self.field.field);
        let f_1 = precomp.non_residue_in_q_minus_one_by_six.clone();
        let f_2 = precomp.non_residue_in_q_squared_minus_one_by_six.clone();

        // use a fact that Fp2 ** (q^2 - 1) == 1 and that 6 | q^2 - 1 and that 6 | q - 1

        // then
        // c1 = Fp2**( (q^1 - 1) / 6) has to be calculated
        // c2 = Fp2**( (q^2 - 1) / 6) has to be calculated
        // c3 = Fp2**( (q^3 - 1) / 6) = Fp2**(( (q^2 - 1) / 6) * q + (q - 1) / 6 ) =
        // = Fp2**( ( (q^2 - 1) / 6) * q ) * Fp2**((q - 1) / 6) =
        // = c2.frobenius(1) * c1
        // c4 = c2.frobenius(0) ** 2
        // c5 = c1 * (c2.frobenius(1) ** 2)
        // c6 = Fp2**( (q^6 - 1) / 6) = Fp2**( (q^2 - 1)/6 * 3) = c2 ** 3
        // c7 = c1 * (c2.frobenius(1) ** 3)
        // c8 = c2.frobenius(0) ** 4
        // c9 = c1 * (c2.frobenius(1) ** 4)
        // c10 = Fp2**( (q^11 - 1) / 6) = c2 ** (q^8 + q^6 + q^4 + q^2 + 1) = c2.frobenius(0) ** 5
        // c11 = Fp2**( (q^11 - 1) / 6) = Fp2**((q - 1) / 6) * c2 ** (q^9 + q^7 + q^5 + q^3 + q) =
        // = c1 * (c2.frobenius(1) ** 5)

        let mut f_3 = f_2.clone();
        f_3.frobenius_map(1);
        f_3.mul_assign(&f_1);

        let f_6 = f_2.pow(&[3u64]);

        let f_4 = Fp2::zero(self.field.field);
        let f_5 = Fp2::zero(self.field.field);

        let f_7 = Fp2::zero(self.field.field);
        let f_8 = Fp2::zero(self.field.field);
        let f_9 = Fp2::zero(self.field.field);
        let f_10 = Fp2::zero(self.field.field);
        let f_11 = Fp2::zero(self.field.field);

        self.frobenius_coeffs_c1 = [f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11];
        self.frobenius_coeffs_are_calculated = true;

        Ok(())
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldExtension for Extension2Over3Over2<'a, E, F> {
    const EXTENSION_DEGREE: usize = 2;
    
    type Element = Fp6<'a, E, F>;

    fn multiply_by_non_residue(&self, el: &mut Self::Element) {
        // IMPORTANT: This only works cause the structure of extension field for Fp12
        // is w^2 - v = 0!
        // take an element in Fp6 that is 3 over 2 and multiply by non-residue
        // (c0 + c1 * v + c2 * v^2)*v with v^3 - xi = 0 -> (c2*xi + c0 * v + c1 * v^2)
        let mut new_c0 = el.c2.clone();
        new_c0.mul_by_nonresidue(&*el.extension_field);
        el.c2 = el.c1.clone();
        el.c1 = el.c0.clone();
        el.c0 = new_c0;
    }

}
