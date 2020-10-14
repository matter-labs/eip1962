use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::field::*;
use crate::extension_towers::fp2::{Extension2, Fp2};
use crate::traits::FieldElement;
use crate::traits::ZeroAndOne;

pub(crate) fn modulus_is_one_mod_four<E: ElementRepr, F: SizedPrimeField<Repr = E>>(field: &F) -> bool {
    const MASK: u64 = 3; // last two bits

    let last_limb = field.modulus().as_ref()[0];

    last_limb & MASK == 1
}

pub(crate) fn modulus_is_three_mod_four<E: ElementRepr, F: SizedPrimeField<Repr = E>>(field: &F) -> bool {
    const MASK: u64 = 3; // last two bits

    let last_limb = field.modulus().as_ref()[0];

    last_limb & MASK == 3
}

pub(crate) fn modulus_is_one_mod_sixteen<E: ElementRepr, F: SizedPrimeField<Repr = E>>(field: &F) -> bool {
    const MASK: u64 = 16; // last four bits

    let last_limb = field.modulus().as_ref()[0];

    last_limb & MASK == 1
}

pub(crate) fn modulus_is_one_mod_four_ext2<E: ElementRepr, F: SizedPrimeField<Repr = E>>(field: &Extension2<E, F>) -> bool {
    modulus_is_one_mod_four(field.field)
}

pub(crate) fn modulus_is_three_mod_four_ext2<E: ElementRepr, F: SizedPrimeField<Repr = E>>(field: &Extension2<E, F>) -> bool {
    modulus_is_three_mod_four(field.field)
}

pub(crate) fn modulus_is_one_mod_sixteen_ext2<E: ElementRepr, F: SizedPrimeField<Repr = E>>(field: &Extension2<E, F>) -> bool {
    modulus_is_one_mod_sixteen(field.field)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LegendreSymbol {
    Zero = 0,
    QuadraticResidue = 1,
    QuadraticNonResidue = -1,
}

pub fn legendre_symbol_fp<'a, E: ElementRepr, F: SizedPrimeField<Repr = E>>(element: &Fp<'a, E, F>) -> LegendreSymbol {
    let mut modulus_minus_one_by_two = *element.field.modulus();
    modulus_minus_one_by_two.shr(1);

    let a = element.pow(&modulus_minus_one_by_two.as_ref());

    if a.is_zero() {
        LegendreSymbol::Zero
    } else if a == Fp::one(element.field) {
        LegendreSymbol::QuadraticResidue
    } else {
        LegendreSymbol::QuadraticNonResidue
    }
}

pub fn legendre_symbol_fp2<'a, E: ElementRepr, F: SizedPrimeField<Repr = E>>(element: &Fp2<'a, E, F>) -> LegendreSymbol {
    let a = element.norm();

    legendre_symbol_fp(&a)
}

fn sqrt_for_one_mod_four<'a, E: ElementRepr, F: SizedPrimeField<Repr = E>>(_element: &Fp<'a, E, F>) -> Option<Fp<'a, E, F>> {
    // TODO, or consider to slow
    
    None
}

pub fn sqrt_for_three_mod_four<'a, E: ElementRepr, F: SizedPrimeField<Repr = E>>(element: &Fp<'a, E, F>) -> Option<Fp<'a, E, F>> {
    // this is a simple case: we compute the power 
    // we know that it's 3 mod 4, so just bit shift

    let mut modulus_minus_three_by_four = *element.field.modulus();
    modulus_minus_three_by_four.shr(2);

    let mut a = element.pow(&modulus_minus_three_by_four.as_ref());

    let mut minus_one = Fp::one(element.field);
    minus_one.negate();

    let mut tmp = a.clone();
    tmp.square();
    tmp.mul_assign(&element);

    if tmp == minus_one {
        None
    } else {
        a.mul_assign(&element);

        Some(a)
    }
}

fn sqrt_for_one_mod_sixteen<'a, E: ElementRepr, F: SizedPrimeField<Repr = E>>(element: &Fp<'a, E, F>) -> Option<Fp<'a, E, F>> {
    // Requires to know a generator and roots of unity
    // postpone for now
    // we know that it's 1 mod 16, so we can use simple bit shifts

    match legendre_symbol_fp(&element) {
        LegendreSymbol::Zero => {
            // it's zero, so clone
            Some(element.clone())
        },
        LegendreSymbol::QuadraticNonResidue => {
            None
        },
        LegendreSymbol::QuadraticResidue => {
            None
        }
    }
}

pub fn sqrt<'a, E: ElementRepr, F: SizedPrimeField<Repr = E>>(element: &Fp<'a, E, F>) -> Option<Fp<'a, E, F>> {
    if modulus_is_three_mod_four(element.field) {
        sqrt_for_three_mod_four(&element)
    } else {
        None
    }
}

pub(crate) fn sqrt_for_three_mod_four_ext2<'a, E: ElementRepr, F: SizedPrimeField<Repr = E>>(element: &Fp2<'a, E, F>) -> Option<Fp2<'a, E, F>> {
    // this is a simple case: we compute the power 
    // we know that it's 3 mod 4, so just bit shifts by 1 or 2 bits are ok

    if element.is_zero() {
        Some(element.clone())
    } else {
        let mut modulus_minus_three_by_four = *element.extension_field.field.modulus();
        modulus_minus_three_by_four.shr(2);
        let mut a1 = element.pow(&modulus_minus_three_by_four.as_ref());

        let mut alpha = a1.clone();
        alpha.square();
        alpha.mul_assign(&element);

        let mut a0 = alpha.clone();
        a0.frobenius_map(1);
        a0.mul_assign(&alpha);

        let one_fp2 = Fp2::one(element.extension_field);

        let mut minus_one_fp2 = one_fp2.clone();
        minus_one_fp2.negate();

        if a0 == minus_one_fp2 {
            None
        } else {
            a1.mul_assign(&element);

            if alpha == minus_one_fp2 {
                let mut tmp = Fp2::zero(element.extension_field);
                tmp.c1 = Fp::one(element.extension_field.field);
                a1.mul_assign(&tmp);
            } else {
                let mut modulus_minus_one_by_two = *element.extension_field.field.modulus();
                modulus_minus_one_by_two.shr(1);

                alpha.add_assign(&one_fp2);

                alpha = alpha.pow(&modulus_minus_one_by_two.as_ref());
                a1.mul_assign(&alpha);
            }

            Some(a1)
        }
    }
}

pub fn sqrt_ext2<'a, E: ElementRepr, F: SizedPrimeField<Repr = E>>(element: &Fp2<'a, E, F>) -> Option<Fp2<'a, E, F>> {
    if modulus_is_three_mod_four_ext2(element.extension_field) {
        sqrt_for_three_mod_four_ext2(&element)
    } else {
        None
    }
}