#[derive(Eq, PartialEq, Clone, Copy, Debug)]
pub enum CurveType {
    Generic,
    AIsMinus3,
    AIsZero,
    BIsZero,
}

use crate::field::SizedPrimeField;
use crate::representation::ElementRepr;
use crate::traits::FieldElement;
use crate::traits::ZeroAndOne;

pub trait CurveParameters {
    type BaseFieldElement: FieldElement + ZeroAndOne;
    fn params(&self) -> <Self::BaseFieldElement as ZeroAndOne>::Params;
}

pub struct CurveOverFpParameters<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> {
    pub(crate) field: &'a F,
}

use crate::fp::Fp;

impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> CurveParameters for &'a CurveOverFpParameters<FE, F> {
    type BaseFieldElement = Fp<'a, FE, F>;
    fn params(&self) -> <Self::BaseFieldElement as ZeroAndOne>::Params {
        self.field
    }
}


pub trait Group: Sized + Clone {
    fn add_assign(&mut self, other: &Self);
    fn add_assign_mixed(&mut self, other: &Self);
    fn sub_assign(&mut self, other: &Self);
    fn negate(&mut self);
    fn double(&mut self);
    fn mul<S: AsRef<[u64]>>(&self, exp: S) -> Self;
    fn wnaf_mul<S: crate::representation::IntoWnaf>(&self, exp: S) -> Self;
    fn is_zero(&self) -> bool;
    fn check_correct_subgroup(&self) -> bool;
}

pub mod curve;
pub mod twist;
pub mod cubic_twist;