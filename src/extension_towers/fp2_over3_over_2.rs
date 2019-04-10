use crate::field::{PrimeFieldElement, SizedPrimeField};
use crate::representation::ElementRepr;
use super::{FieldExtension, ExtensionFieldElement};
use crate::traits::{FieldElement, BitIterator};
use super::fp3_over_2::{Fp6, Extension3Over2};

// this implementation assumes extension using polynomial w^2 - v = 0
pub struct Fp12<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >{
    pub c0: Fp6<'a, E, F>,
    pub c1: Fp6<'a, E, F>,
    pub extension_field: &'a Extension2Over3Over2<'a, E, F>
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >std::fmt::Display for Fp12<'a, E, F> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fq2({} + {} * w)", self.c0, self.c1)
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
    pub fn zero(extension_field: &'a Extension2Over3Over2<'a, E, F>) -> Self {
        let zero = Fp6::zero(extension_field.field);
        
        Self {
            c0: zero.clone(),
            c1: zero,
            extension_field: extension_field
        }
    }

    pub fn one(extension_field: &'a Extension2Over3Over2<'a, E, F>) -> Self {
        let zero = Fp6::zero(extension_field.field);
        let one = Fp6::one(extension_field.field);
        
        Self {
            c0: one,
            c1: zero,
            extension_field: extension_field
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > ExtensionFieldElement for Fp12<'a, E, F> {
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
        c1s.mul_by_nonresidue();
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
        self.c0.mul_by_nonresidue();
        self.c0.add_assign(&aa);
    }

    fn square(&mut self)
    {
        let mut ab = self.c0.clone();
        ab.mul_assign(&self.c1);
        let mut c0c1 = self.c0.clone();
        c0c1.add_assign(&self.c1);
        let mut c0 = self.c1.clone();
        c0.mul_by_nonresidue();
        c0.add_assign(&self.c0);
        c0.mul_assign(&c0c1);
        c0.sub_assign(&ab);
        self.c1 = ab.clone();
        self.c1.add_assign(&ab);
        ab.mul_by_nonresidue();
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

    fn mul_by_nonresidue(&mut self) {
        unreachable!();
    }
}

pub struct Extension2Over3Over2<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > {
    pub non_residue_c0: Fp6<'a, E, F>,
    pub non_residue_c1: Fp6<'a, E, F>,
    pub field: &'a Extension3Over2<'a, E, F>
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldExtension for Extension2Over3Over2<'a, E, F> {
    const EXTENSION_DEGREE: usize = 2;
    
    type Element = Fp12<'a, E, F>;

    fn multiply_by_non_residue(&self, el: &mut Self::Element) {
        
    }

}
