use repr_derive::*;

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

use num_bigint::BigUint;
use num_traits::{One, ToPrimitive, Zero};

use crate::representation::{ElementRepr};

/// Convert BigUint into a vector of 64-bit limbs.
fn biguint_to_fixed_length_u64_vec(mut v: BigUint, limbs: usize) -> Vec<u64> {
    let m = BigUint::one() << 64;
    let mut ret = vec![];

    while v > BigUint::zero() {
        ret.push((&v % &m).to_u64().unwrap());
        v = v >> 64;
    }

    while ret.len() < limbs {
        ret.push(0);
    }

    assert!(ret.len() == limbs);

    ret
}

pub fn biguint_to_u64_vec(mut v: BigUint) -> Vec<u64> {
    let m = BigUint::one() << 64;
    let mut ret = vec![];

    while v > BigUint::zero() {
        ret.push((&v % &m).to_u64().unwrap());
        v = v >> 64;
    }

    ret
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
    fn is_valid_repr(&self, repr: Self::Repr) -> bool;
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
    fn is_valid_repr(&self, repr: Self::Repr) -> bool {
        repr < self.modulus
    }
}

fn calculate_field_dimension(modulus: BigUint) -> Result<((usize, usize), (Vec<u64>, Vec<u64>, Vec<u64>, u64)), ()> {
    let bitlength = modulus.bits();

    let num_limbs = (bitlength / 64) + 1;

    if num_limbs > 16 {
        return Err(());
    }

    // Compute R = 2**(64 * limbs) mod m
    let r = (BigUint::one() << (num_limbs * 64)) % &modulus;
    // Compute R^2 mod m
    let r2 = biguint_to_fixed_length_u64_vec((&r * &r) % &modulus, num_limbs);

    let r = biguint_to_fixed_length_u64_vec(r, num_limbs);
    let modulus = biguint_to_fixed_length_u64_vec(modulus, num_limbs);

    // Compute -m^-1 mod 2**64 by exponentiating by totient(2**64) - 1
    let mut inv = 1u64;
    for _ in 0..63 {
        inv = inv.wrapping_mul(inv);
        inv = inv.wrapping_mul(modulus[0]);
    }
    inv = inv.wrapping_neg();

    Ok(((bitlength, num_limbs), (modulus, r, r2, inv)))
}

pub fn field_from_modulus<R: ElementRepr>(modulus: BigUint) -> Result<PrimeField<R>, ()> {
    let ((bitlength, num_limbs), (modulus, r, r2, inv)) = calculate_field_dimension(modulus)?;
    
    if R::NUM_LIMBS != num_limbs {
        return Err(());
    }

    let mut modulus_repr = R::default();
    let mut r_repr = R::default();
    let mut r2_repr = R::default();
    for (i, ((m_el, r_el), r2_el)) in modulus.into_iter()
                                    .zip(r.into_iter())
                                    .zip(r2.into_iter())
                                    .enumerate() {
        modulus_repr.as_mut()[i] = m_el;
        r_repr.as_mut()[i] = r_el;
        r2_repr.as_mut()[i] = r2_el;
    }

    let concrete = PrimeField {
        mont_power: (num_limbs*4) as u64,
        modulus_bits: bitlength as u64,
        modulus: modulus_repr,
        mont_r: r_repr,
        mont_r2: r2_repr,
        mont_inv: inv,  
    };

    Ok(concrete)
}

pub fn new_field<R: ElementRepr>(modulus: &str, radix: u32) -> Result<PrimeField<R>, ()> {
    use num_traits::Num;
    let modulus = BigUint::from_str_radix(&modulus, radix).unwrap();

    field_from_modulus::<R>(modulus)
}

