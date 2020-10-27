use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::field::*;
use crate::extension_towers::fp2::{Extension2, Fp2};
use crate::traits::FieldElement;
use crate::traits::ZeroAndOne;

#[derive(Clone, Copy, Debug)]
pub struct SqrtContext<'a>{
    pub t_repr: &'a [u64],
    pub t_plus_1_over_2_repr: &'a [u64],
    /// 2^s * t = MODULUS - 1 with t odd
    pub root_of_unity_raw_repr: &'a [u64],
    pub s: u64,
}

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
    const MASK: u64 = 15; // last four bits

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
    // TODO, or consider too slow
    
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

fn sqrt_for_one_mod_sixteen<'a, 'b, E: ElementRepr, F: SizedPrimeField<Repr = E>>(element: &Fp<'a, E, F>, ctx: &SqrtContext<'b>) -> Option<Fp<'a, E, F>> {
    // Requires to know a generator and roots of unity
    // postpone for now
    // we know that it's 1 mod 16, so we can use simple bit shifts

    // Tonelli-Shank's algorithm for q mod 16 = 1
    // https://eprint.iacr.org/2012/685.pdf (page 12, algorithm 5)


    match legendre_symbol_fp(&element) {
        LegendreSymbol::Zero => {
            // it's zero, so clone
            Some(element.clone())
        },
        LegendreSymbol::QuadraticNonResidue => {
            None
        },
        LegendreSymbol::QuadraticResidue => {
            let field = element.field;
            let mut repr = E::default();
            let len_to_copy = ctx.root_of_unity_raw_repr.len();
            debug_assert_eq!(len_to_copy, repr.as_ref().len());
            repr.as_mut()[..len_to_copy].copy_from_slice(ctx.root_of_unity_raw_repr);
            let mut c = if let Ok(c) = Fp::from_raw_repr(field, repr) {
                c
            } else {
                return None
            };
            let mut r = element.pow(ctx.t_plus_1_over_2_repr);
            let mut t = element.pow(ctx.t_repr);
            let mut m = ctx.s;

            let one = Fp::one(field);

            while t != one {
                let mut i = 1;
                {
                    let mut t2i = t;
                    t2i.square();
                    loop {
                        if t2i == one {
                            break;
                        }
                        t2i.square();
                        i += 1;
                    }
                }

                for _ in 0..(m - i - 1) {
                    c.square();
                }
                r.mul_assign(&c);
                c.square();
                t.mul_assign(&c);
                m = i;
            }

            Some(r)
        }
    }
}

pub fn sqrt<'a, 'b, E: ElementRepr, F: SizedPrimeField<Repr = E>>(element: &Fp<'a, E, F>, ctx: Option<SqrtContext<'b>>) -> Option<Fp<'a, E, F>> {
    if modulus_is_three_mod_four(element.field) {
        sqrt_for_three_mod_four(&element)
    } else if modulus_is_one_mod_sixteen(element.field) {
        let ctx = if let Some(ctx) = ctx.as_ref() {
            ctx
        } else {
            return None
        };
        sqrt_for_one_mod_sixteen(element, ctx)
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


pub(crate) fn sqrt_for_one_mod_sixteen_ext_2<'a, 'b, E: ElementRepr, F: SizedPrimeField<Repr = E>>(element: &Fp2<'a, E, F>, ctx: &SqrtContext<'b>) -> Option<Fp2<'a, E, F>> {
    // https://eprint.iacr.org/2012/685.pdf (algorithm 8)

    if element.is_zero() {
        Some(element.clone())
    } else {
        let extension_field = element.extension_field;
        let base_field = extension_field.field;

        let mut two = Fp::one(base_field);
        two.double();
        let two_inv = two.inverse()?;

        let alpha = element.norm();
        let gamma = legendre_symbol_fp(&alpha);

        if gamma == LegendreSymbol::QuadraticNonResidue {
            None
        } else {
            let alpha = sqrt_for_one_mod_sixteen(&alpha, ctx)?;
            let mut delta = element.c0;
            delta.add_assign(&alpha);
            delta.mul_assign(&two_inv);

            let gamma = legendre_symbol_fp(&delta);
            if gamma == LegendreSymbol::QuadraticNonResidue {
                delta = element.c0;
                delta.sub_assign(&alpha);
                delta.mul_assign(&two_inv);
            }

            let x0 = sqrt_for_one_mod_sixteen(&delta, ctx)?;
            let x0_inv = x0.inverse()?;

            let mut x1 = element.c1;
            x1.mul_assign(&x0_inv);
            x1.mul_assign(&two_inv);

            let mut result = Fp2::zero(extension_field);
            result.c0 = x0;
            result.c1 = x1;

            Some(result)
        }
    }
}

pub fn sqrt_ext2<'a, 'b, E: ElementRepr, F: SizedPrimeField<Repr = E>>(element: &Fp2<'a, E, F>, ctx: Option<SqrtContext<'b>>) -> Option<Fp2<'a, E, F>> {
    if modulus_is_three_mod_four_ext2(element.extension_field) {
        sqrt_for_three_mod_four_ext2(&element)
    } else if modulus_is_one_mod_sixteen_ext2(element.extension_field) {
        let ctx = if let Some(ctx) = ctx.as_ref() {
            ctx
        } else {
            return None
        };
        sqrt_for_one_mod_sixteen_ext_2(element, ctx)
    } else {
        None
    }
}

#[cfg(test)]
mod test {
    use crate::field::*;
    use crate::fp::*;
    use crate::traits::*;
    use super::*;

    #[test]
    fn test_bls12_377_sqrt_params() {
        let base_field = &crate::engines::bls12_377::BLS12_377_FIELD;
        let extension = &crate::engines::bls12_377::BLS12_377_EXTENSION_2_FIELD;
        let ctx = crate::engines::bls12_377::BLS12_377_FQ_SQRT_CONTEXT;
        let one = crate::engines::bls12_377::BLS12_377_FP_ONE;
        // simple check that root of unity is actually a root of unity
        let mut repr = U384Repr::default();
        repr.as_mut().copy_from_slice(&ctx.root_of_unity_raw_repr);
        let may_be_one = Fp::from_raw_repr(base_field, repr).unwrap();
        let may_be_one = may_be_one.pow(&[1u64 << ctx.s]);

        assert_eq!(one, may_be_one);
        // now let's compute some square root
        let mut el = one;
        for _ in 0..1024 {
            let may_be_fp_sqrt = sqrt(&el, Some(ctx));
            if let Some(sqrt) = may_be_fp_sqrt {
                let mut tmp = sqrt;
                tmp.square();

                assert_eq!(tmp, el);
            }

            el.add_assign(&one);
        }

        use rand::{Rng, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = &mut XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);
        let mut roots_found = 0;
        // now try on random
        for _ in 0..100000 {
            for limb in repr.as_mut().iter_mut() {
                *limb = rng.gen();
            }

            if let Ok(el) = Fp::from_repr(base_field, repr) {
                if let Some(sqrt) = sqrt(&el, Some(ctx)) {
                    let mut tmp = sqrt;
                    tmp.square();
    
                    assert_eq!(tmp, el);
                    roots_found += 1;
                }
            }
        }

        assert!(roots_found > 0);
        roots_found = 0;

        // and on random fp2
        for _ in 0..1000000 {
            for limb in repr.as_mut().iter_mut() {
                *limb = rng.gen();
            }

            if let Ok(c0) = Fp::from_repr(base_field, repr) {
                for limb in repr.as_mut().iter_mut() {
                    *limb = rng.gen();
                }
                if let Ok(c1) = Fp::from_repr(base_field, repr) {
                    let mut el = Fp2::zero(extension);
                    el.c0 = c0;
                    el.c1 = c1;
                    if let Some(sqrt) = sqrt_ext2(&el, Some(ctx)) {
                        let mut tmp = sqrt;
                        tmp.square();
        
                        assert_eq!(tmp, el);
                        roots_found += 1;
                    }
                }
            }
        }

        assert!(roots_found > 0);
    }
}