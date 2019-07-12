
use crate::traits::{FieldElement};
use crate::weierstrass::Group;

pub(crate) mod bls12;
pub(crate) mod bn;
pub(crate) mod mnt6;
pub(crate) mod mnt4;

#[derive(Eq, PartialEq, Clone, Copy, Debug)]
pub(crate) enum TwistType {
    D,
    M
}

pub(crate) trait PairingEngine: Sized {
    type PairingResult: FieldElement;
    type G1: Group;
    type G2: Group;

    fn pair<'b> (&self, points: &'b [Self::G1], twists: &'b [Self::G2]) -> Option<Self::PairingResult>;
}

pub(crate) fn into_ternary_wnaf(repr: &[u64]) -> Vec<i64> {
    fn is_zero(repr: &[u64]) -> bool {
        for el in repr.iter() {
            if *el != 0 {
                return false;
            }
        }

        true
    }

    fn is_odd(repr: &[u64]) -> bool {
        if repr.len() == 0 {
            return false;
        }

        repr[0] & 1u64 == 1u64
    }

    fn div2(repr: &mut [u64]) {
        let mut t = 0;
        for i in repr.iter_mut().rev() {
            let t2 = *i << 63;
            *i >>= 1;
            *i |= t;
            t = t2;
        }
    }

    fn sub_noborrow(repr: &mut [u64], value: u64) {
        let mut borrow = 0;

        repr[0] = crate::arithmetics::sbb(repr[0], value, &mut borrow);

        for a in repr.iter_mut().skip(1) {
            *a = crate::arithmetics::sbb(*a, 0u64, &mut borrow);
        }
    }

    fn add_nocarry(repr: &mut [u64], value: u64) {
        let mut carry = 0;

        repr[0] = crate::arithmetics::adc(repr[0], value, &mut carry);

        for a in repr.iter_mut().skip(1) {
            *a = crate::arithmetics::adc(*a, 0u64, &mut carry);
        }
    }

    if repr.len() == 0 {
        return vec![];
    }

    let mut res = vec![];
    let mut e = repr.to_vec();

    const WINDOW: u64 = 1u64;
    const MIDPOINT: u64 = 1u64 << WINDOW;
    const MIDPOINT_I64: i64 = MIDPOINT as i64;
    const MASK: u64 = 1u64 << (WINDOW + 1u64);

    while !is_zero(&e) {
        let z: i64;
        if is_odd(&e) {
            z = MIDPOINT_I64 - ((e[0] % MASK) as i64);
            if z >= 0 {
                sub_noborrow(&mut e, z as u64);
            } else {
                add_nocarry(&mut e, (-z) as u64);
            }
        } else {
            z = 0i64;
        }
        res.push(z);
        div2(&mut e);
    }

    res
}