#[macro_use]
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

/// PrimeField is a structure that it instantiated at the runtime 
/// and holds all the necessary information for further arithmetic
/// operations (mainly precompiled Montgommery constants)

use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::{One, ToPrimitive, Zero};

use crate::representation::{ElementRepr, RepresentationDecodingError};
use crate::traits::FieldElement;


/// Convert BigUint into a vector of 64-bit limbs.
fn biguint_to_real_u64_vec(mut v: BigUint, limbs: usize) -> Vec<u64> {
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

fn biguint_num_bits(mut v: BigUint) -> u32 {
    let mut bits = 0;

    while v != BigUint::zero() {
        v = v >> 1;
        bits += 1;
    }

    bits
}

/// This trait represents an element of a field.
pub trait SizedPrimeField: Sized + Send + Sync + std::fmt::Debug
{
    type Repr: ElementRepr;

    fn modulus(&self) -> Self::Repr;
    fn mont_r(&self) -> Self::Repr;
    fn mont_r2(&self) -> Self::Repr;
    fn mont_inv(&self) -> u64;
    fn is_valid_repr(&self, repr: Self::Repr) -> bool;
}

#[derive(Debug)]
pub struct PrimeField<E: ElementRepr> {
    modulus: E,
    mont_r: E,
    mont_r2: E,
    mont_inv: u64
}

// impl<E: ElementRepr> PrimeField<E> {
//     fn zero<'a>(&self) -> PrimeFieldElement<'a, E, Self> {
//         PrimeFieldElement<'a, E, Self> {
//             field: &self,
//             repr: E::default()
//         }
//     }

//     // fn one(&self) -> PrimeFieldElement<'a, E, Self> {
//     //     PrimeFieldElement {
//     //         field: self,
//     //         repr: self.mont_r
//     //     }
//     // }
// }

impl<E: ElementRepr> SizedPrimeField for PrimeField<E> {
    type Repr = E;

    #[inline(always)]
    fn modulus(&self) -> Self::Repr { self.modulus }

    #[inline(always)]
    fn mont_r(&self) -> Self::Repr { self.mont_r }

    #[inline(always)]
    fn mont_r2(&self) -> Self::Repr { self.mont_r2 }

    #[inline(always)]
    fn mont_inv(&self) -> u64 { self.mont_inv }

    #[inline(always)]
    fn is_valid_repr(&self, repr: Self::Repr) -> bool {
        repr < self.modulus
    }
}

pub fn new_field(modulus: &str, radix: u32) -> Option<impl SizedPrimeField> {
    use num_traits::Num;

    let modulus = BigUint::from_str_radix(&modulus, radix).unwrap();
    let bitlength = modulus.bits();

    let mut num_limbs = 0;

    if bitlength < 256 {
        num_limbs = 4;
    } else if bitlength < 320 {
        num_limbs = 5;
    }

    if num_limbs == 0 {
        return None;
    }

    // Compute R = 2**(64 * limbs) mod m
    let r = (BigUint::one() << (num_limbs * 64)) % &modulus;
    // Compute R^2 mod m
    let r2 = biguint_to_real_u64_vec((&r * &r) % &modulus, num_limbs);

    let r = biguint_to_real_u64_vec(r, num_limbs);
    let modulus = biguint_to_real_u64_vec(modulus, num_limbs);

    // Compute -m^-1 mod 2**64 by exponentiating by totient(2**64) - 1
    let mut inv = 1u64;
    for _ in 0..63 {
        inv = inv.wrapping_mul(inv);
        inv = inv.wrapping_mul(modulus[0]);
    }
    inv = inv.wrapping_neg();

    match num_limbs {
        4 => {
            let mut modulus_repr = U256Repr([0,0,0,0]);
            let mut r_repr = U256Repr([0,0,0,0]);
            let mut r2_repr = U256Repr([0,0,0,0]);
            for (i, ((m_el, r_el), r2_el)) in modulus.into_iter()
                                            .zip(r.into_iter())
                                            .zip(r2.into_iter())
                                            .enumerate() {
                modulus_repr.0[i] = m_el;
                r_repr.0[i] = r_el;
                r2_repr.0[i] = r2_el;
            }
            let concrete = PrimeField {
                        modulus: modulus_repr,
                        mont_r: r_repr,
                        mont_r2: r2_repr,
                        mont_inv: inv,  
                    };

            return Some(concrete);
        },
        _ => {},
    }

    None
}

pub struct PrimeFieldElement<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > {
    field: &'a F,
    repr: E
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Clone for PrimeFieldElement<'a, E, F> {
    #[inline(always)]
    fn clone(&self) -> Self {
        Self {
            field: &self.field,
            repr: self.repr
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Ord for PrimeFieldElement<'a, E, F> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        for (a, b) in self.repr.as_ref().iter().rev().zip(other.repr.as_ref().iter().rev()) {
            if a < b {
                return std::cmp::Ordering::Less
            } else if a > b {
                return std::cmp::Ordering::Greater
            }
        }

        std::cmp::Ordering::Equal
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > PartialEq for PrimeFieldElement<'a, E, F> {
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        for (a, b) in self.repr.as_ref().iter().rev().zip(other.repr.as_ref().iter().rev()) {
            if a != b {
                return false;
            }
        }

        true
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > Eq for PrimeFieldElement<'a, E, F> {
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > PartialOrd for PrimeFieldElement<'a, E, F> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > std::fmt::Debug for PrimeFieldElement<'a, E, F>
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "0x")?;
        for i in self.repr.as_ref().iter().rev() {
            write!(f, "{:016x}", *i)?;
        }

        Ok(())
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > std::fmt::Display for PrimeFieldElement<'a, E, F> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "0x")?;
        for i in self.repr.as_ref().iter().rev() {
            write!(f, "{:016x}", *i)?;
        }

        Ok(())
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > PrimeFieldElement<'a, E, F> {
    #[inline(always)]
    fn zero(field: &'a F) -> Self {
        Self {
            field: field,
            repr: E::default()
        }
    }

    #[inline(always)]
    fn one(field: &'a F) -> Self {
        Self {
            field: field,
            repr: field.mont_r()
        }
    }

    fn from_repr(field: &'a F, repr: E) -> Result<Self, RepresentationDecodingError> {
        if field.is_valid_repr(repr) {
            let mut r = Self {
                field: field,
                repr: repr
            };

            let r2 = Self {
                field: field,
                repr: field.mont_r2()
            };

            r.mul_assign(&r2);

            Ok(r)
        } else {
            Err(RepresentationDecodingError::NotInField(format!("{}", repr)))
        }
    }

    pub fn from_be_bytes(field: &'a F, bytes: &[u8]) -> Result<Self, RepresentationDecodingError> {
        let mut repr = E::default();
        repr.read_be(bytes).unwrap();
        Self::from_repr(field, repr)
    }

    /// Subtracts the modulus from this element if this element is not in the
    /// field. Only used interally.
    #[inline(always)]
    fn reduce(&mut self) {
        if !self.field.is_valid_repr(self.repr) {
            self.repr.sub_noborrow(&self.field.modulus());
        }
    }
}

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> > FieldElement for PrimeFieldElement<'a, E, F> {
    /// Returns true iff this element is zero.
    #[inline]
    fn is_zero(&self) -> bool {
        self.repr.is_zero()
    }

    #[inline]
    fn add_assign(&mut self, other: &Self) {
        // This cannot exceed the backing capacity.
        self.repr.add_nocarry(&other.repr);

        // However, it may need to be reduced.
        self.reduce();
    }

    #[inline]
    fn double(&mut self) {
        // This cannot exceed the backing capacity.
        self.repr.mul2();

        // However, it may need to be reduced.
        self.reduce();
    }

    #[inline]
    fn sub_assign(&mut self, other: &Self) {
        // If `other` is larger than `self`, we'll need to add the modulus to self first.
        if other.repr > self.repr {
            self.repr.add_nocarry(&self.field.modulus());
        }

        self.repr.sub_noborrow(&other.repr);
    }

    #[inline]
    fn negate(&mut self) {
        if !self.is_zero() {
            let mut tmp = self.field.modulus();
            tmp.sub_noborrow(&self.repr);
            self.repr = tmp;
        }
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // Guajardo Kumar Paar Pelzl
            // Efficient Software-Implementation of Finite Fields with Applications to Cryptography
            // Algorithm 16 (BEA for Inversion in Fp)

            let one = F::Repr::from(1);

            let modulus = self.field.modulus();
            let mut u = self.repr;
            let mut v = modulus;
            let mut b = Self {
                            field: &self.field,
                            repr: self.field.mont_r2()
                        }; // Avoids unnecessary reduction step.
            let mut c = Self::zero(&self.field);

            while u != one && v != one {
                while u.is_even() {
                    u.div2();

                    if b.repr.is_even() {
                        b.repr.div2();
                    } else {
                        b.repr.add_nocarry(&modulus);
                        b.repr.div2();
                    }
                }

                while v.is_even() {
                    v.div2();

                    if c.repr.is_even() {
                        c.repr.div2();
                    } else {
                        c.repr.add_nocarry(&modulus);
                        c.repr.div2();
                    }
                }

                if v < u {
                    u.sub_noborrow(&v);
                    b.sub_assign(&c);
                } else {
                    v.sub_noborrow(&u);
                    c.sub_assign(&b);
                }
            }

            if u == one {
                Some(b)
            } else {
                Some(c)
            }
        }
    }

    #[inline]
    fn mul_assign(&mut self, other: &Self)
    {
        self.repr.mul_assign(&other.repr, &self.field.modulus(), self.field.mont_inv());
        // self.reduce();
    }

    #[inline]
    fn square(&mut self)
    {
        self.repr.square(&self.field.modulus(), self.field.mont_inv());
        // self.reduce();
    }
}

#[test]
fn test_bn256_field() {
    let field = new_field("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    let modulus = field.modulus();
    println!("{:?}", modulus);

    let one = PrimeFieldElement::one(&field);
    println!("Modulus = {}", one);
    println!("R = {}", field.mont_r());
    println!("R2 = {}", field.mont_r2());
    println!("Inv = {:016x}", field.mont_inv());

    // this is 7 in BE form
    let mut be_repr = vec![0u8; 32];
    be_repr[31] = 7u8;
    let element = PrimeFieldElement::from_be_bytes(&field, &be_repr[..]).unwrap();
    println!("Mont form element = {}", element);
}

#[test]
fn test_babyjubjub_field() {
    let field = new_field("2736030358979909402780800718157159386076813972158567259200215660948447373041", 10).unwrap();
    let modulus = field.modulus();
    println!("{:?}", modulus);

    let one = PrimeFieldElement::one(&field);
    println!("Modulus = {}", one);
    println!("R = {}", field.mont_r());
    println!("R2 = {}", field.mont_r2());
    println!("Inv = {:016x}", field.mont_inv());

    // this is 2 in BE form
    let mut be_repr = vec![0u8; 32];
    be_repr[31] = 2u8;
    let element = PrimeFieldElement::from_be_bytes(&field, &be_repr[..]).unwrap();
    println!("Mont form element = {}", element);
}

