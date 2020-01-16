use crate::traits::FieldElement;
use crate::traits::LsbBitIterator;

pub(crate) fn calculate_window_table<F: FieldElement>(base: &F, window: usize) -> Vec<F> {
    let mut table: Vec<F> = vec![];
    table.reserve(1 << (window - 1));

    let mut acc = base.clone();
    table.push(acc.clone());
    let mut square = acc.clone();
    square.square();

    // pushed 1*G, 3*G, 5*G, etc (notation abuse, it's actually exp)
    for _ in 1..(1 << (window - 1)) {
        acc.mul_assign(&square);
        table.push(acc.clone());
        
    }

    table
}

pub(crate) struct WindowExpBase<F: FieldElement> {
    pub bases: Vec<F>,
    pub window_size: usize,
    one: F
}

impl<F: FieldElement> WindowExpBase<F> {
    pub fn new(base: &F, one: F, window: usize, _num_scalars: usize) -> Self {
        let recommended_window_size = window; // TODO
        let recommended_size_accounding_for_scalars = recommended_window_size;

        let bases = calculate_window_table::<F>(base, recommended_size_accounding_for_scalars);

        Self {
            bases: bases,
            window_size: recommended_size_accounding_for_scalars,
            one: one
        }
    }

    pub fn exponentiate<W: IntoWindows>(&self, scalars: &[W]) -> Vec<F> {
        let mut result = Vec::with_capacity(scalars.len());
        for s in scalars.iter() {
            let wnaf = s.windows(self.window_size as u32);

            let mut res = self.one.clone();
            let mut found_nonzero = false;

            for w in wnaf.into_iter().rev() {
                if w == 0 && found_nonzero {
                    res.square();
                } else if w != 0 {
                    found_nonzero = true;
                    for _ in 0..self.window_size {
                        res.square();
                    }
                    let idx = (w >> 1) as usize;
                    let base = &(self.bases[idx]);
                    res.mul_assign(&base);
                }
            }

            result.push(res)
        }
        result
    }
} 

pub trait IntoWindows {
    fn windows(&self, window: u32) -> Vec<u64>;
}

impl IntoWindows for Vec<u64> {
    fn windows(&self, window: u32) -> Vec<u64> {
        let mut result = vec![];
        let mut found_begining = false;
        let mut w = 0u64;
        let mut bit_count = 0u64;
        let iter = LsbBitIterator::new(&self);
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
                    result.push(0u64);
                    continue;
                }
            }
            if found_begining && bit_count == (window as u64) {
                result.push(w);
                w = 0u64;
                found_begining = false;
                bit_count = 0u64;
            }
        }

        if w != 0 {
            // this is a last chunk if bit length is not divisible by window size
            result.push(w);
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
    use super::IntoWindows;

    #[test]
    fn test_into_windows() {
        let repr = vec![13u64];
        // b1101
        let res = repr.windows(3);
        println!("{:?}", res);
    }

    #[test]
    fn test_windowed_exp() {
        use crate::field::{U256Repr, new_field};
        use crate::fp::Fp;
        use crate::traits::{FieldElement, ZeroAndOne};
        use super::WindowExpBase;

        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let base = Fp::from_repr(&field, U256Repr::from(2)).unwrap();
        let one = Fp::one(&field);

        let scalar = vec![0x43e1f593f0000000,
                    0x2833e84879b97091,
                    0xb85045b68181585d,
                    0x30644e72e131a029];

        let mut square = base.clone();
        square.square();

        let mut cube = base.clone();
        cube.mul_assign(&square);

        let mut fifth = cube.clone();
        fifth.mul_assign(&square);

        let naive_result = base.pow(&scalar[..]);
        let exp_base = WindowExpBase::new(&base, one, 3, 1);

        assert!(exp_base.bases[0] == base);
        assert!(exp_base.bases[1] == cube);
        assert!(exp_base.bases[2] == fifth);
        let mut results = exp_base.exponentiate(&vec![scalar]);
        assert!(results.len() == 1);
        let w_result = results.pop().unwrap();
        assert!(w_result == naive_result);
    }


}