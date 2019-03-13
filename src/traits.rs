use std::fmt;

/// This trait represents an element of a field.
pub trait FieldElement:
    Sized + Eq + Clone + Send + Sync + fmt::Debug + fmt::Display
{
    // /// Returns the zero element of the field, the additive identity.
    // fn zero() -> Self;

    // /// Returns the one element of the field, the multiplicative identity.
    // fn one() -> Self;

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

    // /// Exponentiates this element by a number represented with `u64` limbs,
    // /// least significant digit first.
    // fn pow<S: AsRef<[u64]>>(&self, exp: S) -> Self;
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

/// This trait represents an element of a field that has a square root operation described for it.
pub trait SqrtFieldElement: FieldElement {
    /// Returns the Legendre symbol of the field element.
    fn legendre(&self) -> LegendreSymbol;

    /// Returns the square root of the field element, if it is
    /// quadratic residue.
    fn sqrt(&self) -> Option<Self>;
}

#[derive(Debug, PartialEq)]
pub enum LegendreSymbol {
    Zero = 0,
    QuadraticResidue = 1,
    QuadraticNonResidue = -1,
}



