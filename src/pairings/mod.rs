
use crate::field::SizedPrimeField;
use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, BitIterator};
use crate::weierstrass::Group;
use crate::weierstrass::curve::CurvePoint;
use crate::weierstrass::twist::TwistPoint;
use crate::extension_towers::fp2::{Fp2, Extension2};

pub mod bls12;

// pub struct PreparedTwistPoint<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> {
//     is_infinity: bool,
//     coeffs: Vec<(Fp2<'a, FE, F>, Fp2<'a, FE, F>, Fp2<'a, FE, F>)>
// }

// impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> PreparedTwistPoint<'a, FE, F> {
//     pub fn is_zero(&self) -> bool {
//         self.is_infinity
//     }
// }

pub trait PairingEngine: Sized {
    type PairingResult: FieldElement;
    type G1: Group;
    type G2: Group;

    // fn prepare_twist_point<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>>
    //     (&self, twist_point: &'a TwistPoint<'a, FE, F, GE, G>) -> PreparedTwistPoint<'a, FE, F>;

    fn pair<'b> (&self, point: &'b Self::G1, twist: &'b Self::G2) -> Self::PairingResult;
    
    // fn pair<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>>
        // (&self, point: &'a CurvePoint<'a, FE, F, GE, G>, twist_point: &'a TwistPoint<'a, FE, F, GE, G>) -> Self::PairingResult;
}