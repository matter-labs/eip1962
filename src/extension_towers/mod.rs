pub mod fp2;
pub mod fp3;
pub mod fp4_as_2_over_2;
pub mod fp6_as_2_over_3;
pub mod fp6_as_3_over_2;
pub mod fp12_as_2_over3_over_2;

use crate::fp::Fp;
use crate::field::{SizedPrimeField};
use crate::traits::FieldElement;
use crate::representation::{ElementRepr};
use crate::traits::ZeroAndOne;
use crate::integers::*;


pub(crate) struct Fp2Fp4FrobeniusBaseElements<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >{
    pub(crate) non_residue_in_q_minus_one_by_four: Fp<'a, E, F>,
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Fp2Fp4FrobeniusBaseElements<'a, E, F> {
    pub(crate) fn construct(modulus: &MaxFieldUint, non_residue: &Fp<'a, E, F>) -> Result<Self, ()> {
        if !is_one_mod_four(&modulus) {
            if !crate::features::in_gas_metering() {
                return Err(());
            }
        }

        let power = *modulus >> 2;

        let result = Fp2Fp4FrobeniusBaseElements::<'a, E, F> {
            non_residue_in_q_minus_one_by_four: non_residue.pow(power.as_ref())
        };

        Ok(result)
    }
}

pub(crate) struct Fp3Fp6FrobeniusBaseElements<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >{
    pub(crate) non_residue_in_q_minus_one_by_six: Fp<'a, E, F>,
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Fp3Fp6FrobeniusBaseElements<'a, E, F> {
    pub(crate) fn construct(modulus: &MaxFieldUint, non_residue: &Fp<'a, E, F>) -> Result<Self, ()> {
        // we perform an explicit division, so we do not check that
        // modulus == 1 mod 6 due to expensive division
    
        let one = MaxFieldUint::from(1u64);
        let six = MaxFieldUint::from(6u64);

        // 1
        let f_1 = {
            let power = *modulus - one;
            let (power, rem) = power.div_mod(six);
            if !rem.is_zero() {
                if !crate::features::in_gas_metering() {
                    return Err(());
                }
            }

            non_residue.pow(power.as_ref())
        };

        let result = Fp3Fp6FrobeniusBaseElements::<'a, E, F> {
            non_residue_in_q_minus_one_by_six: f_1,
        };

        Ok(result)
    }
}

pub(crate) struct Fp6Fp12FrobeniusBaseElements<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >{
    pub(crate) non_residue_in_q_minus_one_by_six: fp2::Fp2<'a, E, F>,
    pub(crate) non_residue_in_q_squared_minus_one_by_six: fp2::Fp2<'a, E, F>,
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Fp6Fp12FrobeniusBaseElements<'a, E, F> {
    pub(crate) fn construct(modulus: &MaxFieldUint, non_residue: &fp2::Fp2<'a, E, F>) -> Result<Self, ()> {
        let modulus = MaxFieldSquaredUint::from(modulus.as_ref());
        
        // we perform an explicit division, so we do not check that
        // modulus == 1 mod 6 due to expensive division

        let mut q_power = modulus;
        let one = MaxFieldSquaredUint::from(1u64);
        let six = MaxFieldSquaredUint::from(6u64);

        // 1
        let f_1 = {
            let power = q_power - one;
            let (power, rem) = power.div_mod(six);
            if !rem.is_zero() {
                if !crate::features::in_gas_metering() {
                    return Err(());
                }
            }

            non_residue.pow(power.as_ref())
        };

        // 2
        let f_2 = {
            q_power = q_power.adaptive_multiplication(modulus);
            let power = q_power - one;
            let (power, rem) = power.div_mod(six);
            if !rem.is_zero() {
                if !crate::features::in_gas_metering() {
                    return Err(());
                }
            }

            non_residue.pow(power.as_ref())
        };

        let result = Fp6Fp12FrobeniusBaseElements::<'a, E, F> {
            non_residue_in_q_minus_one_by_six: f_1,
            non_residue_in_q_squared_minus_one_by_six: f_2
        };

        Ok(result)
    }
}

pub(crate) fn is_non_nth_root<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>
(
    element: & Fp<'a, FE, F>,
    modulus: &MaxFieldUint,
    n: u64
) -> bool {
    if element.is_zero() {
        return false;
    }
    let one = MaxFieldUint::from(1u64);
    let mut power = *modulus;
    power -= one;

    let divisor = MaxFieldUint::from(n);
    let (power, rem) = power.div_mod(divisor);
    if !rem.is_zero() {
        if !crate::features::in_gas_metering() {
            return false;
        }
    }

    let l = if crate::features::in_gas_metering() {
        element.pow(&vec![core::u64::MAX; power.as_ref().len()])
    } else {
        element.pow(power.as_ref())
    };

    let one = Fp::one(element.field);
    if l == one {
        return false;
    } else {
        return true;
    }
}

pub(crate) fn is_non_nth_root_fp2<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>
(
    element: & self::fp2::Fp2<'a, FE, F>,
    modulus: &MaxFieldUint,
    n: u64
) -> bool {
    if element.is_zero() {
        return false;
    }
    let mut power = MaxFieldSquaredUint::from(modulus.as_ref());
    // power *= power;
    power = power.adaptive_multiplication(power);
    power -= MaxFieldSquaredUint::from(1u64);

    let divisor = MaxFieldSquaredUint::from(n);
    let (power, rem) = power.div_mod(divisor);
    if !rem.is_zero() {
        if !crate::features::in_gas_metering() {
            return false;
        }
    }

    let l = if crate::features::in_gas_metering() {
        element.pow(&vec![core::u64::MAX; power.as_ref().len()])
    } else {
        element.pow(power.as_ref())
    };
    
    let one = self::fp2::Fp2::one(element.extension_field);
    if l == one {
        return false;
    } else {
        return true;
    }
}

pub(crate) fn is_one_mod_two
(
    modulus: &MaxFieldUint,
) -> bool {
    modulus.low_u64() & 1u64 == 1u64
}

#[allow(dead_code)]
pub(crate) fn is_one_mod_three
(
    modulus: &MaxFieldUint,
) -> bool {
    let one = MaxFieldUint::from(1u64);
    let three = MaxFieldUint::from(3u64);
    let mut tmp = *modulus;
    tmp -= one;
    let (_, rem) = tmp.div_mod(three);

    rem.is_zero()
}

pub(crate) fn is_one_mod_four
(
    modulus: &MaxFieldUint,
) -> bool {
    const LAST_TWO_BITS_MASK :u64 = 3;
    modulus.low_u64() & LAST_TWO_BITS_MASK == 1u64
}

pub(crate) fn is_one_mod_six
(
    modulus: &MaxFieldUint,
) -> bool {
    if !is_one_mod_two(modulus) {
        return false;
    }
    let div_2 = *modulus >> 1;
    
    is_one_mod_three(&div_2)
}