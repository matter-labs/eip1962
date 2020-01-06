use eth_pairings_repr_derive::*;

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "4"]
struct U256(U256Repr);

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "5"]
struct U320(U320Repr);

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "6"]
struct U384(U384Repr); 

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "7"]
struct U448(U448Repr);

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "8"]
struct U512(U512Repr);

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "9"]
struct U576(U576Repr);

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "10"]
struct U640(U640Repr);

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "11"]
struct U704(U704Repr);

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "12"]
struct U768(U768Repr);

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "13"]
struct U832(U832Repr);

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "14"]
struct U896(U896Repr);

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "15"]
struct U960(U960Repr);

#[derive(ElementRepresentation)]
#[NumberOfLimbs = "16"]
struct U1024(U1024Repr);

/// PrimeField is a structure that it instantiated at the runtime 
/// and holds all the necessary information for further arithmetic
/// operations (mainly precompiled Montgommery constants)

use crate::representation::{ElementRepr};
use crate::constants::*;

pub(crate) fn slice_to_fixed_length_u64_vec<S: AsRef<[u64]>>(v: S, limbs: usize) -> Vec<u64> {
    let as_ref = v.as_ref();
    let mut ret = Vec::with_capacity(limbs);

    let mut num_words = as_ref.len();
    for v in as_ref.iter().rev() {
        if *v == 0 {
            num_words -= 1;
        } else {
            break;
        }
    }

    debug_assert!(num_words <= limbs);

    ret.extend(&as_ref[0..num_words]);

    if ret.len() < limbs {
        ret.resize(limbs, 0u64);
    }

    assert!(ret.len() == limbs);

    ret
}

pub(crate) fn slice_to_u64_vec<S: AsRef<[u64]>>(v: S) -> Vec<u64> {
    let as_ref = v.as_ref();

    let mut num_words = as_ref.len();
    for v in as_ref.iter().rev() {
        if *v == 0 {
            num_words -= 1;
        } else {
            break;
        }
    }

    let mut ret = Vec::with_capacity(num_words);

    ret.extend(&as_ref[0..num_words]);

    ret
}

// use crate::arrayvec::{ArrayVec, Array};

// pub(crate) fn slice_to_fixed_size_array<S: AsRef<[u64]>, A: Array<Item = u64>>(v: S) -> ArrayVec<A> {
//     let as_ref = v.as_ref();

//     let mut num_words = as_ref.len();
//     for v in as_ref.iter().rev() {
//         if *v == 0 {
//             num_words -= 1;
//         } else {
//             break;
//         }
//     }

//     debug_assert!(num_words <= A::CAPACITY);

//     let mut ret = ArrayVec::<A>::new();
//     ret.try_extend_from_slice(&as_ref[0..num_words]).expect("must extend with enough capacity");

//     ret
// }

fn num_words(number: &MaxFieldSquaredUint) -> usize {
    let bits = number.bits();

    (bits + 63) / 64
}

/// This trait represents an element of a field.
pub trait SizedPrimeField: Sized + Send + Sync + std::fmt::Debug
{
    type Repr: ElementRepr;

    fn mont_power(&self) -> u64;
    fn modulus_bits(&self) -> u64;
    fn modulus(&self) -> &Self::Repr;
    fn mont_r(&self) -> &Self::Repr;
    fn mont_r2(&self) -> &Self::Repr;
    fn mont_inv(&self) -> u64;
    fn is_valid_repr(&self, repr: &Self::Repr) -> bool;
}

#[derive(Debug)]
pub struct PrimeField<E: ElementRepr> {
    mont_power: u64,
    modulus_bits: u64,
    modulus: E,
    mont_r: E,
    mont_r2: E,
    mont_inv: u64
}

impl<E: ElementRepr> SizedPrimeField for PrimeField<E> {
    type Repr = E;

    #[inline(always)]
    fn mont_power(&self) -> u64 { self.mont_power }

    #[inline(always)]
    fn modulus_bits(&self) -> u64 { self.modulus_bits }

    #[inline(always)]
    fn modulus(&self) -> &Self::Repr { &self.modulus }

    #[inline(always)]
    fn mont_r(&self) -> &Self::Repr { &self.mont_r }

    #[inline(always)]
    fn mont_r2(&self) -> &Self::Repr { &self.mont_r2 }

    #[inline(always)]
    fn mont_inv(&self) -> u64 { self.mont_inv }

    #[inline(always)]
    fn is_valid_repr(&self, repr: &Self::Repr) -> bool {
        repr < &self.modulus
    }
}

pub(crate) fn calculate_num_limbs(bitlength: usize) -> Result<usize, ()> {
    let mut num_limbs = (bitlength / 64) + 1;
    if num_limbs < 4 {
        num_limbs = 4;
    }

    if num_limbs > 16 {
        return Err(());
    }

    Ok(num_limbs)
}

pub fn field_from_modulus<R: ElementRepr>(modulus: &MaxFieldUint) -> Result<PrimeField<R>, ()> {
    let bitlength = modulus.bits();
    let num_limbs = calculate_num_limbs(bitlength)?;

    let modulus = MaxFieldSquaredUint::from(modulus.as_ref());

    if R::NUM_LIMBS != num_limbs {
        return Err(());
    }

    let r = (MaxFieldSquaredUint::one() << ((num_limbs * 64) as u32)) % modulus;
    if num_words(&r) > R::NUM_LIMBS {
        return Err(());
    }

    let r2 = (r * r) % modulus;
    if num_words(&r2) > R::NUM_LIMBS {
        return Err(());
    }

    let modulus_lowest_limb = modulus.as_ref()[0];

    let mut inv = 1u64;
    for _ in 0..63 {
        inv = inv.wrapping_mul(inv);
        inv = inv.wrapping_mul(modulus_lowest_limb);
    }
    inv = inv.wrapping_neg();

    let mut modulus_repr = R::default();
    let mut r_repr = R::default();
    let mut r2_repr = R::default();

    let modulus_ref = modulus.as_ref();
    let r_ref = r.as_ref();
    let r2_ref = r2.as_ref();

    for (i, ((m_el, r_el), r2_el)) in modulus_repr.as_mut().iter_mut()
                                    .zip(r_repr.as_mut().iter_mut())
                                    .zip(r2_repr.as_mut().iter_mut())
                                    .enumerate() {
        *m_el = modulus_ref[i];
        *r_el = r_ref[i];
        *r2_el = r2_ref[i];
    }

    let concrete = PrimeField {
        mont_power: (num_limbs*64) as u64,
        modulus_bits: bitlength as u64,
        modulus: modulus_repr,
        mont_r: r_repr,
        mont_r2: r2_repr,
        mont_inv: inv,  
    };

    Ok(concrete)
}

#[cfg(test)]
pub(crate) fn new_field<R: ElementRepr>(modulus: &str, radix: usize) -> Result<PrimeField<R>, ()> {
    use num_bigint::BigUint;
    use num_traits::*;
    let modulus = BigUint::from_str_radix(modulus, radix as u32).unwrap();
    let modulus = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());

    field_from_modulus(&modulus)
}