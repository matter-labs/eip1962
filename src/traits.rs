use std::fmt;

/// This trait represents an element of a field.
pub trait FieldElement:
    Sized 
    + Eq 
    + Clone 
    + Send 
    + Sync 
    + fmt::Debug 
    + fmt::Display
{
    /// Returns true iff this element is zero.
    fn is_zero(&self) -> bool;

    /// Squares this element.
    fn square(&mut self);

    /// Doubles this element.
    fn double(&mut self);

    /// Negates this element.
    fn negate(&mut self);

    /// Adds another element to this element.
    fn add_assign(&mut self, other: &Self);

    /// Subtracts another element from this element.
    fn sub_assign(&mut self, other: &Self);

    /// Multiplies another element by this element.
    fn mul_assign(&mut self, other: &Self);

    /// Computes the multiplicative inverse of this element, if nonzero.
    fn inverse(&self) -> Option<Self>;

    /// Exponentiates this element by a number represented with `u64` limbs,
    /// least significant digit first.
    fn pow<S: AsRef<[u64]>>(&self, exp: S) -> Self;

    fn conjugate(&mut self);

    fn mul_by_nonresidue<EXT: FieldExtension<Element = Self>>(&mut self, for_extension: &EXT);

    fn frobenius_map(&mut self, power: usize);
}

pub trait ZeroAndOne {
    type Params;
    fn zero(f: Self::Params) -> Self;
    fn one(f: Self::Params) -> Self;
    // fn get_params(&self) -> Self::Params;
}

pub trait FieldExtension {
    const EXTENSION_DEGREE: usize;

    type Element;

    fn multiply_by_non_residue(&self, el: &mut Self::Element);
}

#[derive(Debug)]
pub struct BitIterator<E> {
    t: E,
    n: usize,
}

impl<E: AsRef<[u64]>> BitIterator<E> {
    pub fn new(t: E) -> Self {
        let n = t.as_ref().len() * 64;

        BitIterator { t, n }
    }
}

impl<E: AsRef<[u64]>> Iterator for BitIterator<E> {
    type Item = bool;

    fn next(&mut self) -> Option<bool> {
        if self.n == 0 {
            None
        } else {
            self.n -= 1;
            let part = self.n / 64;
            let bit = self.n - (64 * part);

            Some(self.t.as_ref()[part] & (1 << bit) > 0)
        }
    }
}


// this bit iterator skips initial zeroes until reaches MSB
#[derive(Debug)]
pub struct MsbBitIterator<E> {
    t: E,
    n: usize,
}

impl<E: AsRef<[u64]>> MsbBitIterator<E> {
    pub fn new(t: E) -> Self {
        let r = t.as_ref();
        let mut n = r.len() * 64;
        let mut found_one = false;
        for limb in (0..r.len()).rev() {
            if found_one {
                break;
            }
            if r[limb] == 0 {
                n -= 64;
                continue;
            }
            // counting from MSB in limb
            for bit in 0..64 {
                if r[limb] & (1 << (63 - bit)) == 0 {
                    n -= 1;
                } else {
                    found_one = true;
                    break;
                }
            }
        }

        MsbBitIterator { t, n }
    }
}

impl<E: AsRef<[u64]>> Iterator for MsbBitIterator<E> {
    type Item = bool;

    fn next(&mut self) -> Option<bool> {
        if self.n == 0 {
            None
        } else {
            self.n -= 1;
            let part = self.n / 64;
            let bit = self.n - (64 * part);

            Some(self.t.as_ref()[part] & (1 << bit) > 0)
        }
    }
}

// this is LSB bit iterator
#[derive(Debug)]
pub struct LsbBitIterator<E> {
    t: E,
    n: usize,
    max: usize
}

impl<E: AsRef<[u64]>> LsbBitIterator<E> {
    pub fn new(t: E) -> Self {
        let max = t.as_ref().len() * 64;
        let n = 0;
        LsbBitIterator { t, n, max}
    }
}

impl<E: AsRef<[u64]>> Iterator for LsbBitIterator<E> {
    type Item = bool;

    fn next(&mut self) -> Option<bool> {
        if self.n == self.max {
            None
        } else {
            let part = self.n / 64;
            let bit = self.n - (64 * part);
            self.n += 1;

            Some(self.t.as_ref()[part] & (1 << bit) > 0)
        }
    }
}

pub use crate::weierstrass::Group;

// /// This trait represents an element of a field that has a square root operation described for it.
// pub trait SqrtFieldElement: FieldElement {
//     /// Returns the Legendre symbol of the field element.
//     fn legendre(&self) -> LegendreSymbol;

//     /// Returns the square root of the field element, if it is
//     /// quadratic residue.
//     fn sqrt(&self) -> Option<Self>;
// }

// #[derive(Debug, PartialEq)]
// pub enum LegendreSymbol {
//     Zero = 0,
//     QuadraticResidue = 1,
//     QuadraticNonResidue = -1,
// }


#[cfg(test)]
mod bit_iter_tests {
    #[test]
    fn test_msb_iter() {
        use super::MsbBitIterator;
        let word: u64 = 0x0103;
        let iter = MsbBitIterator::new([word, 0]);
        let bits: Vec<bool> = iter.collect();
        assert!(bits.len() == 9);
        assert!(bits == vec![true,
                            false, false, false, false,
                            false, false, true, true]);
    }
}
