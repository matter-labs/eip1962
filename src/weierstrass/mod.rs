#[derive(Eq, PartialEq, Clone, Copy, Debug)]
pub enum CurveType {
    Generic,
    AIsMinus3,
    AIsZero,
    BIsZero,
}

pub trait Group: Sized {
    fn add_assign(&mut self, other: &Self);
    fn add_assign_mixed(&mut self, other: &Self);
    fn sub_assign(&mut self, other: &Self);
    fn negate(&mut self);
    fn double(&mut self);
    fn mul<S: AsRef<[u64]>>(&self, exp: S) -> Self;
    fn wnaf_mul<S: crate::representation::IntoWnaf>(&self, exp: S) -> Self;
    fn is_zero(&self) -> bool;
}

pub mod curve;
pub mod twist;
pub mod cubic_twist;