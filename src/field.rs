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
use crate::integers::*;

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

fn num_words(number: &MaxFieldSquaredUint) -> usize {
    let bits = number.bits();

    (bits + 63) / 64
}

/// This trait represents an element of a field.
pub trait SizedPrimeField: Sized + Send + Sync + std::fmt::Debug + 'static + Copy + Clone
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
    pub mont_inv: u64,
    pub modulus: E,
    pub mont_r: E,
    pub mont_r2: E,
    pub mont_power: u64,
    pub modulus_bits: u64,
}

impl<E: ElementRepr> Clone for PrimeField<E> {
    fn clone(&self) -> Self {
        Self {
            mont_power: self.mont_power,
            modulus_bits: self.modulus_bits,
            modulus: self.modulus,
            mont_r: self.mont_r,
            mont_r2: self.mont_r2,
            mont_inv: self.mont_inv
        }
    }
}

impl<E: ElementRepr> Copy for PrimeField<E> {}

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
    if modulus.low_u64() & 1 == 0 {
        // modulus is even
        return Err(());
    }

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

    let r2 = r.adaptive_multiplication(r);
    let r2 = r2 % modulus;

    // debug_assert!((r * r) % modulus == r2);

    // let r2 = (r * r) % modulus;
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

#[cfg(test)]
mod test {
    use num_bigint::BigUint;

    #[test]
    fn test_64_subs_vs_division() {
        use crate::integers::*;
        use num_traits::Num;

        let modulus_biguint = BigUint::from_str_radix("90478214930942163331059450080490962782645063145694412829751296064846357295922334563151067188227020604570222293387098730257627925580549709781572499823579422083606813980978533736029952476604411769448188461808725324994973770830757245236797352800919619593992035038055539141605771028907859068447994637683496172821", 10).unwrap();

        assert_eq!(modulus_biguint.bits(), 1024);

        let modulus = MaxFieldSquaredUint::from_big_endian(&modulus_biguint.to_bytes_be());

        let num_limbs = 16;

        let t = MaxFieldSquaredUint::one() << ((num_limbs * 64) as u32);

        const REPEATS: usize = 10000;
        
        let start = std::time::Instant::now();

        let mut i = 0;
        for _ in 0..REPEATS {
            let _ = t % modulus;
            i += 1;
        }

        println!("Time for one division is {} ns", start.elapsed().as_nanos() / (REPEATS as u128));

        assert!(i == REPEATS);
    }

    #[test]
    fn test_field_construction_speed() {
        use crate::integers::*;
        use num_traits::Num;

        let modulus_biguint = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        assert_eq!(modulus_biguint.bits(), 381);

        let modulus = MaxFieldUint::from_big_endian(&modulus_biguint.to_bytes_be());

        const REPEATS: usize = 10000;
        
        let start = std::time::Instant::now();

        let mut i = 0;
        for _ in 0..REPEATS {
            let _ = super::field_from_modulus::<super::U384Repr>(&modulus).unwrap();
            i += 1;
        }

        println!("Time to construct 381 field is {} ns", start.elapsed().as_nanos() / (REPEATS as u128));

        assert!(i == REPEATS);
    }
}