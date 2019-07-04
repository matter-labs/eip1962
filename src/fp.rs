use crate::representation::{ElementRepr, RepresentationDecodingError};
use crate::traits::FieldElement;
use crate::traits::BitIterator;
use crate::traits::FieldExtension;
use crate::field::SizedPrimeField;
use crate::traits::ZeroAndOne;

pub struct Fp<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > {
    pub(crate) field: &'a F,
    pub(crate) repr: E
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Clone for Fp<'a, E, F> {
    #[inline(always)]
    fn clone(&self) -> Self {
        Self {
            field: &self.field,
            repr: self.repr
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Ord for Fp<'a, E, F> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // use non-montgommery form
        let modulus = self.field.modulus();
        let mont_inv = self.field.mont_inv();
        let this = self.repr.into_normal_repr(&modulus, mont_inv);
        let that = other.repr.into_normal_repr(&modulus, mont_inv);
        for (a, b) in this.as_ref().iter().rev().zip(that.as_ref().iter().rev()) {
            if a < b {
                return std::cmp::Ordering::Less
            } else if a > b {
                return std::cmp::Ordering::Greater
            }
        }

        std::cmp::Ordering::Equal
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > PartialEq for Fp<'a, E, F> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        for (a, b) in self.repr.as_ref().iter().rev().zip(other.repr.as_ref().iter().rev()) {
            if a != b {
                return false;
            }
        }

        true
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Eq for Fp<'a, E, F> {
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > PartialOrd for Fp<'a, E, F> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > std::fmt::Debug for Fp<'a, E, F>
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "0x")?;
        // for i in self.repr.as_ref().iter().rev() {
        for i in self.into_repr().as_ref().iter().rev() {
            write!(f, "{:016x}", *i)?;
        }

        Ok(())
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > std::fmt::Display for Fp<'a, E, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "0x")?;
        // for i in self.repr.as_ref().iter().rev() {
        for i in self.into_repr().as_ref().iter().rev() {
            write!(f, "{:016x}", *i)?;
        }

        Ok(())
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Fp<'a, E, F> {
    // #[inline(always)]
    // pub fn zero(field: &'a F) -> Self {
    //     Self {
    //         field: field,
    //         repr: E::default()
    //     }
    // }

    // #[inline(always)]
    // pub fn one(field: &'a F) -> Self {
    //     Self {
    //         field: field,
    //         repr: field.mont_r()
    //     }
    // }

    pub fn from_repr(field: &'a F, repr: E) -> Result<Self, RepresentationDecodingError> {
        if field.is_valid_repr(repr) {
            let mut r = Self {
                field: field,
                repr: repr
            };

            let r2 = Self {
                field: field,
                repr: field.mont_r2().clone()
            };

            r.mul_assign(&r2);

            Ok(r)
        } else {
            Err(RepresentationDecodingError::NotInField(format!("{}", repr)))
        }
    }

    pub fn into_repr(&self) -> E {
        let modulus = self.field.modulus();
        let mont_inv = self.field.mont_inv();
        self.repr.into_normal_repr(&modulus, mont_inv)
    }

    pub fn from_be_bytes(field: &'a F, bytes: &[u8], allow_padding: bool) -> Result<Self, RepresentationDecodingError> {
        let mut repr = E::default();
        if bytes.len() >= repr.as_ref().len() * 8 {
            repr.read_be(bytes).map_err(|e| RepresentationDecodingError::NotInField(format!("Failed to read big endian bytes, {}", e)))?;
        } else {
            if allow_padding {
                let mut padded = vec![0u8; repr.as_ref().len() * 8 - bytes.len()];
                padded.extend_from_slice(bytes);
                repr.read_be(&padded[..]).map_err(|e| RepresentationDecodingError::NotInField(format!("Failed to read big endian bytes, {}", e)))?;
            } else {
                repr.read_be(&bytes[..]).map_err(|e| RepresentationDecodingError::NotInField(format!("Failed to read big endian bytes without padding, {}", e)))?;
            }
        }
        Self::from_repr(field, repr)
    }

    /// Subtracts the modulus from this element if this element is not in the
    /// field. Only used interally.
    #[inline(always)]
    fn reduce(&mut self) {
        if !self.field.is_valid_repr(self.repr) {
            self.repr.sub_noborrow(&self.field.modulus());
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > ZeroAndOne for Fp<'a, E, F> {
    type Params = &'a F;

    #[inline(always)]
    fn zero(field: &'a F) -> Self {
        Self {
            field: field,
            repr: E::default()
        }
    }

    #[inline(always)]
    fn one(field: &'a F) -> Self {
        Self {
            field: field,
            repr: *field.mont_r()
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldElement for Fp<'a, E, F> {
    /// Returns true iff this element is zero.
    #[inline]
    fn is_zero(&self) -> bool {
        self.repr.is_zero()
    }

    #[inline]
    fn add_assign(&mut self, other: &Self) {
        // This cannot exceed the backing capacity.
        self.repr.add_nocarry(&other.repr);

        // However, it may need to be reduced.
        self.reduce();
    }

    #[inline]
    fn double(&mut self) {
        // This cannot exceed the backing capacity.
        self.repr.mul2();

        // However, it may need to be reduced.
        self.reduce();
    }

    #[inline]
    fn sub_assign(&mut self, other: &Self) {
        // If `other` is larger than `self`, we'll need to add the modulus to self first.
        if other.repr > self.repr {
            self.repr.add_nocarry(&self.field.modulus());
        }

        self.repr.sub_noborrow(&other.repr);
    }

    #[inline]
    fn negate(&mut self) {
        if !self.is_zero() {
            let mut tmp = self.field.modulus().clone();
            tmp.sub_noborrow(&self.repr);
            self.repr = tmp;
        }
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // Guajardo Kumar Paar Pelzl
            // Efficient Software-Implementation of Finite Fields with Applications to Cryptography
            // Algorithm 16 (BEA for Inversion in Fp)

            let one = F::Repr::from(1);

            let modulus = self.field.modulus();
            let mut u = self.repr;
            let mut v = modulus.clone();
            let mut b = Self {
                field: &self.field,
                repr: self.field.mont_r2().clone()
            }; // Avoids unnecessary reduction step.
            let mut c = Self::zero(&self.field);

            while u != one && v != one {
                while u.is_even() {
                    u.div2();

                    if b.repr.is_even() {
                        b.repr.div2();
                    } else {
                        b.repr.add_nocarry(&modulus);
                        b.repr.div2();
                    }
                }

                while v.is_even() {
                    v.div2();

                    if c.repr.is_even() {
                        c.repr.div2();
                    } else {
                        c.repr.add_nocarry(&modulus);
                        c.repr.div2();
                    }
                }

                if v < u {
                    u.sub_noborrow(&v);
                    b.sub_assign(&c);
                } else {
                    v.sub_noborrow(&u);
                    c.sub_assign(&b);
                }
            }

            if u == one {
                Some(b)
            } else {
                Some(c)
            }
        }
    }

    #[inline]
    fn mul_assign(&mut self, other: &Self)
    {
        self.repr.mont_mul_assign(&other.repr, &self.field.modulus(), self.field.mont_inv());
    }

    #[inline]
    fn square(&mut self)
    {
        self.repr.mont_square(&self.field.modulus(), self.field.mont_inv());
    }

    fn pow<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Self::one(&self.field);

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

    fn conjugate(&mut self) {
        unreachable!();
    }

    fn frobenius_map(&mut self, _power: usize) {
        unreachable!();
    }
}