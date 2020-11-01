use crate::traits::FieldElement;
use crate::representation::IntoWnaf;

#[allow(dead_code)]
pub(crate) fn calculate_wnaf_table<F: FieldElement>(base: &F, window: usize) -> Vec<F> {
    let mut table: Vec<F> = vec![];
    table.reserve(1 << (window - 1));

    let mut acc = *base;

    let mut double = acc;
    double.double();

    // pushed 1*G, 3*G, 5*G, etc

    for _ in 0..(1 << (window - 1)) {
        let to_push = acc;
        table.push(to_push);
        acc.add_assign(&double);
    }

    table
}

#[allow(dead_code)]
pub(crate) struct WnafBase<F: FieldElement> {
    pub bases: Vec<F>,
    pub window_size: usize,
    zero: F
}

#[allow(dead_code)]
impl<F: FieldElement> WnafBase<F> {
    pub fn new(base: &F, zero: F, window: usize, _num_scalars: usize) -> Self {
        let recommended_window_size = window; // TODO
        let recommended_size_accounding_for_scalars = recommended_window_size;

        let bases = calculate_wnaf_table::<F>(base, recommended_size_accounding_for_scalars);

        Self {
            bases: bases,
            window_size: recommended_size_accounding_for_scalars,
            zero: zero
        }
    }

    // pub fn exponentiate<W: IntoWnaf>(&self, scalars: &[W]) -> Vec<F> {
    //     let mut result = vec![];
    //     for s in scalars.iter() {
    //         let wnaf = s.wnaf(self.window_size as u32);
    //         println!("wnaf = {:?}", wnaf);

    //         let mut res = self.zero.clone();
    //         let mut found_nonzero = false;

    //         for w in wnaf.into_iter().rev() {
    //             if found_nonzero {
    //                 res.double();
    //             }
    //             if w != 0 {
    //                 println!("w = {}", w);
    //                 found_nonzero = true;
    //                 if w > 0 {
    //                     let idx = (w >> 1) as usize;
    //                     let base = &(self.bases[idx]);
    //                     res.add_assign(&base);
    //                 } else {
    //                     let idx = ((-w) >> 1) as usize;
    //                     let base = &(self.bases[idx]);
    //                     res.sub_assign(&base);
    //                 }
    //             }
    //         }

    //         result.push(res)
    //     }
    //     result
    // }
} 

impl<'a> IntoWnaf for &'a [u64] {
// impl IntoWnaf for Vec<u64> {
    fn wnaf(&self, window: u32) -> Vec<i64> {

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
    
        if self.len() == 0 {
            return vec![];
        }
    
        let mut res = Vec::with_capacity(self.len() * 64);
        let mut e = self.to_vec();

        let max = (1 << window) as i64;
        let midpoint = (1 << (window - 1)) as i64;
        let modulus_mask = ((1 << window) - 1) as u64;

        while !is_zero(&e) {
            let z: i64;
            if is_odd(&e) {
                let masked_bits = (e[0] & modulus_mask) as i64;
                if masked_bits > midpoint {
                    z = masked_bits - max;
                    add_nocarry(&mut e, (-z) as u64);
                } else {
                    z = masked_bits;
                    sub_noborrow(&mut e, z as u64);
                }
            } else {
                z = 0i64;
            }
            res.push(z);
            div2(&mut e);
        }

        res
    
        // let window: u64 = window as u64;
        // let midpoint: u64 = 1u64 << window;
        // let midpoint_i64: i64 = midpoint as i64;
        // let mask: u64 = (1u64 << (window + 1u64)) - 1;
    
        // while !is_zero(&e) {
        //     let z: i64;
        //     if is_odd(&e) {
        //         z = midpoint_i64 - ((e[0] & mask) as i64);
        //         if z >= 0 {
        //             sub_noborrow(&mut e, z as u64);
        //         } else {
        //             add_nocarry(&mut e, (-z) as u64);
        //         }
        //     } else {
        //         z = 0i64;
        //     }
        //     res.push(z);
        //     div2(&mut e);
        // }
    
        // res
    }
}

#[cfg(test)]
mod tests {
    use crate::representation::IntoWnaf;

    #[test]
    fn test_wnaf_form_calculation() {
        let repr = vec![13u64];
        // b1101
        let res = (&repr[..]).wnaf(3);
        println!("{:?}", res);
    }

    // #[test]
    // fn test_correctness_wnaf_mul() {
    //     use crate::field::{U256Repr, new_field};
    //     use crate::fp::Fp;
    //     use crate::traits::{FieldElement};
    //     use super::WnafBase;
    //     use crate::traits::ZeroAndOne;

    //     let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
    //     let base = Fp::from_repr(&field, U256Repr::from(4)).unwrap();
    //     let zero = Fp::zero(&field);

    //     // let scalar = vec![0x43e1f593f0000000,
    //     //             0x2833e84879b97091,
    //     //             0xb85045b68181585d,
    //     //             0x30644e72e131a029];

    //     let scalar = vec![2u64];

    //     let naive_result = base.pow(&scalar[..]);
    //     let wnaf_base = WnafBase::new(&base, zero, 3, 1);
    //     let mut results = wnaf_base.exponentiate(&vec![&scalar[..]]);
    //     assert!(results.len() == 1);
    //     let wnaf_result = results.pop().unwrap();
    //     assert!(wnaf_result == naive_result);
    // }


}