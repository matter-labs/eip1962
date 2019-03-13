use super::field::{PrimeFieldElement, SizedPrimeField};
use super::representation::ElementRepr;

#[derive(Eq, PartialEq, Clone, Copy)]
pub enum CurveType {
    Generic,
    AIsMinus3,
    AIsZero,
    BIsZero,
}

pub struct WeierstrassCurve<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> {
    field: &'a F,
    group: &'a G,
    a: PrimeFieldElement<'a, FE, F>,
    b: PrimeFieldElement<'a, FE, F>,
    curve_type: CurveType
}

pub struct GroupElement<'a, R: ElementRepr, G: SizedPrimeField<Repr = R>> {
    group: &'a G,
    repr: R,
}

pub trait Group {
    fn add_assign(&mut self, other: &Self);
    fn sub_assign(&mut self, other: &Self);
    fn negate(&mut self);
    fn mul<S: AsRef<[u64]>>(&self, exp: S) -> Self;
}

impl<'a, R: ElementRepr, G: SizedPrimeField<Repr = R>> GroupElement<'a, R, G> {
    // add_assign(&mut self, other: &Self) 
}