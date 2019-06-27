pub mod fp2;
pub mod fp3;
pub mod fp4_as_2_over_2;
pub mod fp6_as_2_over_3;
pub mod fp6_as_3_over_2;
pub mod fp12_as_2_over3_over_2;

use crate::fp::Fp;
use crate::field::{SizedPrimeField};
use crate::traits::FieldElement;
use crate::representation::{ElementRepr, LegendreSymbol};

pub(crate) fn legendre_symbol<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, S: AsRef<[u64]>>
(
    element: & Fp<'a, FE, F>,
    modulus_minus_one_by_2: S
) -> LegendreSymbol {
    if element.is_zero() {
        return LegendreSymbol::Zero;
    }
    let l = element.pow(modulus_minus_one_by_2);
    let one = Fp::one(element.field);
    if l == one {
        return LegendreSymbol::QuadraticResidue;
    } else {
        return LegendreSymbol::QuadraticNonResidue;
    }
}

pub(crate) fn is_non_cube<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, S: AsRef<[u64]>>
(
    element: & Fp<'a, FE, F>,
    modulus_minus_one_by_3: S
) -> bool {
    if element.is_zero() {
        return false;
    }
    let l = element.pow(modulus_minus_one_by_3);
    let one = Fp::one(element.field);
    if l == one {
        return false;
    } else {
        return true;
    }
}

pub(crate) fn is_non_cube_fp2<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, S: AsRef<[u64]>>
(
    element: & self::fp2::Fp2<'a, FE, F>,
    modulus_squared_minus_one_by_3: S
) -> bool {
    if element.is_zero() {
        return false;
    }
    let l = element.pow(modulus_squared_minus_one_by_3);
    let one = self::fp2::Fp2::one(element.extension_field);
    if l == one {
        return false;
    } else {
        return true;
    }
}