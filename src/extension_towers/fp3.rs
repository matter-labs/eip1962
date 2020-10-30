use crate::fp::Fp;
use crate::field::{SizedPrimeField};
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, BitIterator, FieldExtension, ZeroAndOne};
use super::Fp3Fp6FrobeniusBaseElements;

// this implementation assumes extension using polynomial u^3 + m = 0
pub struct Fp3<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >{
    pub c0: Fp<'a, E, F>,
    pub c1: Fp<'a, E, F>,
    pub c2: Fp<'a, E, F>,
    pub extension_field: &'a Extension3<'a, E, F>
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >std::fmt::Display for Fp3<'a, E, F> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Fq3({} + {} * u + {} * u^2)", self.c0, self.c1, self.c2)
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >std::fmt::Debug for Fp3<'a, E, F> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Fq3({} + {} * u) + {} * u^2", self.c0, self.c1, self.c2)
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Clone for Fp3<'a, E, F> {
    #[inline(always)]
    fn clone(&self) -> Self {
        Self{
            c0: self.c0,
            c1: self.c1,
            c2: self.c2,
            extension_field: self.extension_field
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Copy for Fp3<'a, E, F> {}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > PartialEq for Fp3<'a, E, F> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && 
        self.c1 == other.c1 &&
        self.c2 == other.c2
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Eq for Fp3<'a, E, F> {
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Fp3<'a, E, F> {
    pub fn mul_by_fp(&mut self, element: &Fp<'a, E, F>) {
        self.c0.mul_assign(&element);
        self.c1.mul_assign(&element);
        self.c2.mul_assign(&element);
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > ZeroAndOne for Fp3<'a, E, F> {
    type Params = &'a Extension3<'a, E, F>;

    fn zero(extension_field: &'a Extension3<'a, E, F>) -> Self {
        let zero = Fp::zero(extension_field.field());
        
        Self {
            c0: zero,
            c1: zero,
            c2: zero,
            extension_field: extension_field
        }
    }

    fn one(extension_field: &'a Extension3<'a, E, F>) -> Self {
        let zero = Fp::zero(extension_field.field());
        let one = Fp::one(extension_field.field());
        
        Self {
            c0: one,
            c1: zero,
            c2: zero,
            extension_field: extension_field
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldElement for Fp3<'a, E, F> {
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
        if self.is_zero() {
            None
        } else {
            let mut t0 = self.c0;
            t0.square();
            let mut t1 = self.c1;
            t1.square();
            let mut t2 = self.c2;
            t2.square();
            let mut t3 = self.c0;
            t3.mul_assign(&self.c1);
            let mut t4 = self.c0;
            t4.mul_assign(&self.c2);
            let mut t5 = self.c1;
            t5.mul_assign(&self.c2);
            let mut n5 = t5;
            n5.mul_by_nonresidue(self.extension_field);

            let mut s0 = t0;
            s0.sub_assign(&n5);
            let mut s1 = t2;
            s1.mul_by_nonresidue(self.extension_field);
            s1.sub_assign(&t3);
            let mut s2 = t1;
            s2.sub_assign(&t4); // typo in paper referenced above. should be "-" as per Scott, but is "*"

            let mut a1 = self.c2;
            a1.mul_assign(&s1);
            let mut a2 = self.c1;
            a2.mul_assign(&s2);
            let mut a3 = a1;
            a3.add_assign(&a2);
            a3.mul_by_nonresidue(self.extension_field);
            let mut t6 = self.c0;
            t6.mul_assign(&s0);
            t6.add_assign(&a3);
            let t6 = t6.inverse();
            if t6.is_none() {
                return None;
            }

            let t6 = t6.expect("is not None");

            let mut c0 = t6;
            c0.mul_assign(&s0);
            let mut c1 = t6;
            c1.mul_assign(&s1);
            let mut c2 = t6;
            c2.mul_assign(&s2);


            Some(Self{
                c0, 
                c1, 
                c2,
                extension_field: self.extension_field
            })
        }
    }

    fn mul_assign(&mut self, other: &Self)
    {
        let a = other.c0;
        let b = other.c1;
        let c = other.c2;

        let d = self.c0;
        let e = self.c1;
        let f = self.c2;

        let mut ad = d;
        ad.mul_assign(&a);
        let mut be = e;
        be.mul_assign(&b);
        let mut cf = f;
        cf.mul_assign(&c);

        let mut t0 = b;
        t0.add_assign(&c);

        let mut x = e;
        x.add_assign(&f);
        x.mul_assign(&t0);
        x.sub_assign(&be);
        x.sub_assign(&cf);

        let mut t0 = a;
        t0.add_assign(&b);

        let mut y = d;
        y.add_assign(&e);
        y.mul_assign(&t0);
        y.sub_assign(&ad);
        y.sub_assign(&be);

        let mut t0 = a;
        t0.add_assign(&c);

        let mut z = d;
        z.add_assign(&f);
        z.mul_assign(&t0);
        z.sub_assign(&ad);
        z.add_assign(&be);
        z.sub_assign(&cf);

        let mut t0 = x;
        t0.mul_by_nonresidue(self.extension_field);

        self.c0 = t0;
        self.c0.add_assign(&ad);

        let mut t0 = cf;
        t0.mul_by_nonresidue(self.extension_field);

        self.c1 = t0;
        self.c1.add_assign(&y);

        self.c2 = z;
    }

    fn square(&mut self)
    {
        let a = self.c0;
        let b = self.c1;
        let c = self.c2;

        let mut s0 = a;
        s0.square();
        let mut ab = a;
        ab.mul_assign(&b);
        let mut s1 = ab;
        s1.double();
        let mut s2 = a;
        s2.sub_assign(&b);
        s2.add_assign(&c);
        s2.square();
        let mut bc = b;
        bc.mul_assign(&c);
        let mut s3 = bc;
        s3.double();
        let mut s4 = c;
        s4.square();

        self.c0 = s0;
        let mut t0 = s3;
        t0.mul_by_nonresidue(self.extension_field);
        self.c0.add_assign(&t0);

        self.c1 = s1;
        let mut t1 = s4;
        t1.mul_by_nonresidue(self.extension_field);
        self.c1.add_assign(&t1);

        self.c2 = s1;
        self.c2.add_assign(&s2);
        self.c2.add_assign(&s3);
        self.c2.sub_assign(&s0);
        self.c2.sub_assign(&s4);
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
        assert!(self.extension_field.frobenius_coeffs_are_calculated);
        self.c1.mul_assign(&self.extension_field.frobenius_coeffs_c1[power % 3]);
        self.c2.mul_assign(&self.extension_field.frobenius_coeffs_c2[power % 3]);
    }
}

pub struct Extension3<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > {
    // pub(crate) field: &'a F,
    pub(crate) non_residue: Fp<'a, E, F>,
    pub(crate) frobenius_coeffs_c1: [Fp<'a, E, F>; 3],
    pub(crate) frobenius_coeffs_c2: [Fp<'a, E, F>; 3],
    pub(crate) frobenius_coeffs_are_calculated: bool
}

use crate::integers::*;

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Extension3<'a, E, F> {
    #[inline(always)]
    pub fn field(&self) -> &'a F {
        self.non_residue.field
    }

    pub (crate) fn new(non_residue: Fp<'a, E, F>) -> Self {
        let field = non_residue.field;

        let zeros = [Fp::zero(field), Fp::zero(field), Fp::zero(field)];
        
        Self {
            non_residue,
            // field: & field,
            frobenius_coeffs_c1: zeros,
            frobenius_coeffs_c2: zeros,
            frobenius_coeffs_are_calculated: false
        }
    }

    pub(crate) fn calculate_frobenius_coeffs_optimized(
        &mut self,
        modulus: &MaxFieldUint,
    ) -> Result<(), ()> {
        // NON_RESIDUE**(((q^0) - 1) / 3)
        let non_residue = &self.non_residue;
        let f_0 = Fp::one(self.field());

        // use a fact that Fp ** (q - 1) == 1 and that 3 | q - 1

        // then
        // c1 = Fp**( (q^1 - 1) / 3) has to be calculated
        // c2 = Fp**( (q^2 - 1) / 3) = Fp**( ((q - 1)/3) *(q+1)) = 
        // = c1 * c1.frobenius(1) = c1 ** 2

        let modulus = MaxFieldUint::from(modulus.as_ref());
        let q_power = modulus;
        let one = MaxFieldUint::from(1u64);
        let three = MaxFieldUint::from(3u64);

        // 1
        let f_1 = {
            let power = q_power - one;
            let (power, rem) = power.div_mod(three);
            if !rem.is_zero() {
                if !crate::features::in_gas_metering() {
                    return Err(());
                }
            }

            non_residue.pow(power.as_ref())
        };

        // c1 * c1.frobenius(1) == c1^2
        let mut f_2 = f_1;
        f_2.square();

        let f_0_c2 = f_0;

        let mut f_1_c2 = f_1;
        let mut f_2_c2 = f_2;

        f_1_c2.square();
        f_2_c2.square();

        self.frobenius_coeffs_c1 = [f_0, f_1, f_2];
        self.frobenius_coeffs_c2 = [f_0_c2, f_1_c2, f_2_c2];
        self.frobenius_coeffs_are_calculated = true;

        Ok(())
    }

    pub(crate) fn calculate_frobenius_coeffs_with_precomp(
        &mut self,
        precomp: &Fp3Fp6FrobeniusBaseElements<'a, E, F>
    ) -> Result<(), ()> {    
        let f_0 = Fp::one(self.field());

        // precomputation has it by 6, so square
        let mut f_1 = precomp.non_residue_in_q_minus_one_by_six;
        f_1.square();
        
        // c1 * c1.frobenius(1) == c1^2
        let mut f_2 = f_1;
        // f_2.frobenius_map(1); // we could leave it formally, but it's an identity
        f_2.square();

        let f_0_c2 = f_0;

        let mut f_1_c2 = f_1;
        let mut f_2_c2 = f_2;

        f_1_c2.square();
        f_2_c2.square();

        self.frobenius_coeffs_c1 = [f_0, f_1, f_2];
        self.frobenius_coeffs_c2 = [f_0_c2, f_1_c2, f_2_c2];
        self.frobenius_coeffs_are_calculated = true;

        Ok(())
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldExtension for Extension3<'a, E, F> {
    const EXTENSION_DEGREE: usize = 3;
    
    type Element = Fp<'a, E, F>;

    fn multiply_by_non_residue(&self, el: &mut Self::Element) {
        // this is simply a multiplication by non-residue that is Fp element cause everything else 
        // is covered in explicit formulas for multiplications for Fp3
        el.mul_assign(&self.non_residue);
    }
}
