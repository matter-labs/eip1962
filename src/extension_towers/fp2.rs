use crate::fp::Fp;
use crate::field::{SizedPrimeField};
use crate::representation::ElementRepr;
use super::{FieldExtension};
use crate::traits::{FieldElement, BitIterator};


// this implementation assumes extension using polynomial u^2 + 1 = 0
// multiply_by_non_residue function in the extension field actually depends
// on the nature of the higher extension, but is supplied at runtime
// e.g. BLS12-381 curve multiplies by (u+1), while BN254 uses (u+9)
pub struct Fp2<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >{
    pub c0: Fp<'a, E, F>,
    pub c1: Fp<'a, E, F>,
    pub extension_field: &'a Extension2<'a, E, F>
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >std::fmt::Display for Fp2<'a, E, F> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fq2({} + {} * u)", self.c0, self.c1)
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >std::fmt::Debug for Fp2<'a, E, F> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "Fq2({} + {} * u)", self.c0, self.c1)
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Clone for Fp2<'a, E, F> {
    #[inline(always)]
    fn clone(&self) -> Self {
        Self{
            c0: self.c0.clone(),
            c1: self.c1.clone(),
            extension_field: self.extension_field
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > PartialEq for Fp2<'a, E, F> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && 
        self.c1 == other.c1
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Eq for Fp2<'a, E, F> {
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Fp2<'a, E, F> {
    pub fn zero(extension_field: &'a Extension2<'a, E, F>) -> Self {
        let zero = Fp::zero(extension_field.field);
        
        Self {
            c0: zero.clone(),
            c1: zero,
            extension_field: extension_field
        }
    }

    pub fn one(extension_field: &'a Extension2<'a, E, F>) -> Self {
        let zero = Fp::zero(extension_field.field);
        let one = Fp::one(extension_field.field);
        
        Self {
            c0: one,
            c1: zero,
            extension_field: extension_field
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldElement for Fp2<'a, E, F> {
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
        let mut t1 = self.c1.clone();
        t1.square();
        let mut t0 = self.c0.clone();
        t0.square();
        t0.add_assign(&t1);
        t0.inverse().map(|t| {
            let mut tmp = Self {
                c0: self.c0.clone(),
                c1: self.c1.clone(),
                extension_field: self.extension_field
            };
            tmp.c0.mul_assign(&t);
            tmp.c1.mul_assign(&t);
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
        self.c0 = aa;
        self.c0.sub_assign(&bb);
    }

    fn square(&mut self)
    {
        let mut ab = self.c0.clone();
        ab.mul_assign(&self.c1);
        let mut c0c1 = self.c0.clone();
        c0c1.add_assign(&self.c1);
        let mut c0 = self.c1.clone();
        c0.negate();
        c0.add_assign(&self.c0);
        c0.mul_assign(&c0c1);
        c0.sub_assign(&ab);
        self.c1 = ab.clone();
        self.c1.add_assign(&ab);
        c0.add_assign(&ab);
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
        self.extension_field.multiply_by_non_residue(self);
    }

    fn frobenius_map(&mut self, power: usize) {
        unimplemented!();
    }
}

pub struct Extension2<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > {
    pub non_residue_c0: Fp<'a, E, F>,
    pub non_residue_c1: Fp<'a, E, F>,
    pub field: &'a F
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldExtension for Extension2<'a, E, F> {
    const EXTENSION_DEGREE: usize = 2;
    
    type Element = Fp2<'a, E, F>;

    fn multiply_by_non_residue(&self, el: &mut Self::Element) {
        let mut as_element = el.clone();
        as_element.c0 = self.non_residue_c0.clone();
        as_element.c1 = self.non_residue_c1.clone();
        el.mul_assign(&as_element);
    }

}
