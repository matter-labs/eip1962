use crate::representation::{ElementRepr};
use crate::traits::FieldElement;
use crate::field::{SizedPrimeField, PrimeFieldElement};

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >PrimeFieldElement<'a, E, F> {
    pub fn mont_inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // The Montgomery Modular Inverse - Revisited

            // Phase 1
            let modulus = self.field.modulus();
            let mut u = modulus;
            let mut v = self.repr;
            let mut r = F::Repr::from(0);
            let mut s = F::Repr::from(1);
            let mut k = 0u64;

            while !v.is_zero() {
                if u.is_even() {
                    u.div2();
                    s.mul2();
                } else if v.is_even() {
                    v.div2();
                    r.mul2();
                } else if u > v {
                    u.sub_noborrow(&v);
                    u.div2();
                    r.add_nocarry(&s);
                    s.mul2();
                } else if v >= u {
                    v.sub_noborrow(&u);
                    v.div2();
                    s.add_nocarry(&r);
                    r.mul2();
                }

                k += 1;
            }

            if r >= modulus {
                r.sub_noborrow(&modulus);
            }

            let mut tmp = modulus;
            tmp.sub_noborrow(&r);

            r = tmp;

            // phase 2

            for _ in 0..(k-self.field.mont_power()) {
                if r.is_even() {
                    r.div2();
                } else {
                    r.add_nocarry(&modulus);
                    r.div2();
                }
            }

            Some(PrimeFieldElement::from_repr(self.field, r).unwrap())
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::FieldElement;
    use crate::field::{PrimeFieldElement};
    use crate::field::U256Repr;
    #[test]
    fn test_mont_inverse() {
        use crate::field::new_field;
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        // this is 7 in BE form
        let mut be_repr = vec![0u8; 32];
        be_repr[31] = 7u8;
        let element = PrimeFieldElement::from_be_bytes(&field, &be_repr[..]).unwrap();
        let inverse = element.inverse().unwrap();
        let mont_inverse = element.mont_inverse().unwrap();
        println!("inv = {}, mont_inv = {}", inverse, mont_inverse);
        assert!(inverse == mont_inverse);
    }
}