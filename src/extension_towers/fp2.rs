use crate::fp::Fp;
use crate::field::{SizedPrimeField};
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, BitIterator, FieldExtension};


// this implementation assumes extension using polynomial u^2 + m = 0
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

    pub fn mul_by_fp(&mut self, element: &Fp<'a, E, F>) {
        self.c0.mul_assign(&element);
        self.c1.mul_assign(&element);
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
        if self.is_zero() {
            None
        } else {
            // Guide to Pairing-based Cryptography, Algorithm 5.19.
            // v0 = c0.square()
            let mut v0 = self.c0.clone();
            v0.square();
            // v1 = c1.square()
            let mut v1 = self.c1.clone();
            v1.square();
            // v0 = v0 - beta * v1
            let mut v1_by_nonresidue = v1.clone();
            v1_by_nonresidue.mul_by_nonresidue(self.extension_field);
            v0.sub_assign(&v1_by_nonresidue);
            v0.inverse().map(|v1| {
                let mut c0 = self.c0.clone();
                c0.mul_assign(&v1);
                let mut c1 = self.c1.clone();
                c1.mul_assign(&v1);
                c1.negate();

                Self {
                    c0: c0, 
                    c1: c1,
                    extension_field: self.extension_field
                }
            })
        }
    }

    fn mul_assign(&mut self, other: &Self)
    {
        let mut v0 = self.c0.clone();
        v0.mul_assign(&other.c0);
        let mut v1 = self.c1.clone();
        v1.mul_assign(&other.c1);

        self.c1.add_assign(&self.c0);
        let mut t0 = other.c0.clone();
        t0.add_assign(&other.c1);
        self.c1.mul_assign(&t0);
        self.c1.sub_assign(&v0);
        self.c1.sub_assign(&v1);
        self.c0 = v0;
        v1.mul_by_nonresidue(self.extension_field);
        self.c0.add_assign(&v1);
    }

    fn square(&mut self)
    {
        // v0 = c0 - c1
        let mut v0 = self.c0.clone();
        v0.sub_assign(&self.c1);
        // v3 = c0 - beta * c1
        let mut v3 = self.c0.clone();
        let mut t0 = self.c1.clone();
        t0.mul_by_nonresidue(self.extension_field);
        v3.sub_assign(&t0);
        // v2 = c0 * c1
        let mut v2 = self.c0.clone();
        v2.mul_assign(&self.c1);

        // v0 = (v0 * v3) + v2
        v0.mul_assign(&v3);
        v0.add_assign(&v2);

        self.c1 = v2.clone();
        self.c1.double();
        self.c0 = v0;
        v2.mul_by_nonresidue(self.extension_field);
        self.c0.add_assign(&v2);

    }

    fn conjugate(&mut self) {
        unreachable!();
        // self.c1.negate();
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
        self.c1.mul_assign(&self.extension_field.frobenius_coeffs_c1[power % 2]);
    }
}

// For example, BLS12-381 has non-residue = -1;
pub struct Extension2<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > {
    pub field: &'a F,
    pub non_residue: Fp<'a, E, F>,
    pub frobenius_coeffs_c1: [Fp<'a, E, F>; 2],
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldExtension for Extension2<'a, E, F> {
    const EXTENSION_DEGREE: usize = 2;
    
    type Element = Fp<'a, E, F>;

    fn multiply_by_non_residue(&self, el: &mut Self::Element) {
        el.mul_assign(&self.non_residue);
    }
}
