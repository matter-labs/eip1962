pub trait FieldExtension {
    const EXTENSION_DEGREE: usize;

    type Element;

    fn multiply_by_non_residue(&self, el: &mut Self::Element);
}

pub mod fp2;
pub mod fp6_as_3_over_2;
pub mod fp12_as_2_over3_over_2;

// /// This trait represents an element of a field.
// pub trait ExtensionFieldElement:
//     Sized 
//     + Eq 
//     + Clone 
//     + Send 
//     + Sync 
//     // + std::fmt::Debug 
//     // + std::fmt::Display
// {
//     /// Returns true iff this element is zero.
//     fn is_zero(&self) -> bool;

//     /// Squares this element.
//     fn square(&mut self);

//     /// Doubles this element.
//     fn double(&mut self);

//     /// Negates this element.
//     fn negate(&mut self);

//     /// Adds another element to this element.
//     fn add_assign(&mut self, other: &Self);

//     /// Subtracts another element from this element.
//     fn sub_assign(&mut self, other: &Self);

//     /// Multiplies another element by this element.
//     fn mul_assign(&mut self, other: &Self);

//     /// Computes the multiplicative inverse of this element, if nonzero.
//     fn inverse(&self) -> Option<Self>;

//     /// Exponentiates this element by a number represented with `u64` limbs,
//     /// least significant digit first.
//     fn pow<S: AsRef<[u64]>>(&self, exp: S) -> Self;

//     /// Exponentiates this element by a number represented with `u64` limbs,
//     /// least significant digit first.
//     fn conjugate(&mut self);

//     fn mul_by_nonresidue(&mut self);
// }