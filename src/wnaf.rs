use crate::representation::ElementRepr;
use crate::traits::FieldElement;
use crate::representation::IntoWnaf;
use crate::traits::LsbBitIterator;

pub(crate) fn calculate_wnaf_table<F: FieldElement>(base: &F, window: usize) -> Vec<F> {
    let mut table: Vec<F> = vec![];
    table.reserve(1 << (window - 1));

    let mut acc = base.clone();

    let mut double = acc.clone();
    double.double();

    // pushed 1*G, 3*G, 5*G, etc

    for _ in 0..(1 << (window - 1)) {
        let to_push = acc.clone();
        table.push(to_push);
        acc.add_assign(&double);
    }

    table
}

pub(crate) struct WnafBase<F: FieldElement> {
    pub bases: Vec<F>,
    pub window_size: usize,
    zero: F
}

impl<F: FieldElement> WnafBase<F> {
    pub fn new(base: &F, zero: F, window: usize, num_scalars: usize) -> Self {
        let recommended_window_size = window; // TODO
        let recommended_size_accounding_for_scalars = recommended_window_size;

        let bases = calculate_wnaf_table::<F>(base, recommended_size_accounding_for_scalars);

        Self {
            bases: bases,
            window_size: recommended_size_accounding_for_scalars,
            zero: zero
        }
    }

    pub fn exponentiate<W: IntoWnaf>(&self, scalars: &[W]) -> Vec<F> {
        let mut result = vec![];
        for s in scalars.iter() {
            let wnaf = s.wnaf(self.window_size as u32);
            println!("wnaf = {:?}", wnaf);

            let mut res = self.zero.clone();
            let mut found_nonzero = false;

            for w in wnaf.into_iter().rev() {
                if found_nonzero {
                    res.double();
                }
                if w != 0 {
                    println!("w = {}", w);
                    found_nonzero = true;
                    if w > 0 {
                        let idx = (w >> 1) as usize;
                        let base = &(self.bases[idx]);
                        res.add_assign(&base);
                    } else {
                        let idx = ((-w) >> 1) as usize;
                        let base = &(self.bases[idx]);
                        res.sub_assign(&base);
                    }
                }
            }

            result.push(res)
        }
        result
    }
} 

impl IntoWnaf for Vec<u64> {
    fn wnaf(&self, window: u32) -> Vec<i64> {
        let mut result = vec![];
        let mut found_begining = false;
        let mut w = 0u64;
        let mut bit_count = 0u64;
        let middle_point = 1u64 << window;
        let top_point = middle_point << 1;
        let mut iter = LsbBitIterator::new(&self);
        for b in iter.into_iter() {
            if b {
                if found_begining {
                    w |= 1u64 << bit_count;
                    bit_count += 1;
                } else {
                    found_begining = true;
                    w |= 1u64 << bit_count;
                    bit_count += 1;
                }
            } else {
                if found_begining {
                    bit_count += 1;
                } else {
                    result.push(0i64);
                    continue;
                }
            }
            if found_begining && bit_count == ((window + 1) as u64) {
                if w >= middle_point {
                    let r = -((top_point - w) as i64);
                    result.push(r);
                } else {
                    let r = w as i64;
                    result.push(r);
                }
                w = 0u64;
                found_begining = false;
                bit_count = 0u64;
            }
        }

        if w != 0 {
            // this is a last chunk if bit length is not divisible by window size
            if w > middle_point {
                let r = -((top_point - w) as i64);
                result.push(r);
            } else {
                let r = w as i64;
                result.push(r);
            }
        }

        for _ in 0..result.len() {
            if let Some(v) = result.pop() {
                if v == 0 {
                    continue;
                } else {
                    result.push(v);
                    break;
                }
            }
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use crate::representation::IntoWnaf;

    #[test]
    fn test_into_wnaf() {
        let repr = vec![7u64];
        let res = repr.wnaf(2);
        println!("{:?}", res);
    }

    #[test]
    fn test_exp() {
        use crate::field::{U256Repr, new_field};
        use crate::fp::Fp;
        use crate::traits::{FieldElement};
        use super::WnafBase;

        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let base = Fp::from_repr(&field, U256Repr::from(4)).unwrap();
        let zero = Fp::zero(&field);

        // let scalar = vec![0x43e1f593f0000000,
        //             0x2833e84879b97091,
        //             0xb85045b68181585d,
        //             0x30644e72e131a029];

        let scalar = vec![2u64];

        let naive_result = base.pow(&scalar[..]);
        let wnaf_base = WnafBase::new(&base, zero, 3, 1);
        let mut results = wnaf_base.exponentiate(&vec![scalar]);
        assert!(results.len() == 1);
        let wnaf_result = results.pop().unwrap();
        assert!(wnaf_result == naive_result);
    }


}