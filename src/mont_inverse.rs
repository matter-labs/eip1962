use crate::representation::{ElementRepr};
use crate::traits::FieldElement;
use crate::field::{SizedPrimeField};
use crate::fp::Fp;

impl<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >Fp<'a, E, F> {
    pub(crate) fn mont_inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // The Montgomery Modular Inverse - Revisited

            // Phase 1
            let modulus = *self.field.modulus();
            let mut u = modulus;
            let mut v = self.repr;
            let mut r = F::Repr::from(0);
            let mut s = F::Repr::from(1);
            let mut k = 0u64;

            let mut found = false;
            for _ in 0..E::NUM_LIMBS*128 {
            // while !v.is_zero() {
                if v.is_zero() {
                    found = true;
                    break;
                }
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
            if !found {
                return None;
            }

            if r >= modulus {
                r.sub_noborrow(&modulus);
            }

            let mut tmp = modulus;
            tmp.sub_noborrow(&r);

            r = tmp;

            // r = (a ^ -1) 2*(k - m)

            // phase 2

            let mont_power_param = self.field.mont_power();
            if k > mont_power_param {
                for _ in 0..(k - mont_power_param) {
                    if r.is_even() {
                        r.div2();
                    } else {
                        r.add_nocarry(&modulus);
                        r.div2();
                    }
                }
            } else {
                for _ in 0..(mont_power_param - k) {
                    r.mul2();
                    if r >= modulus {
                        r.sub_noborrow(&modulus);
                    }
                }
            }

            // if k < mont_power_param {
            //     return None;
            // }

            // for _ in 0..(k - mont_power_param) {
            //     if r.is_even() {
            //         r.div2();
            //     } else {
            //         r.add_nocarry(&modulus);
            //         r.div2();
            //     }
            // }

            let el = Fp::from_repr(self.field, r);
            if el.is_err() {
                return None;
            }
            let el = el.expect("guaranteed to exist");

            Some(el)
        }
    }

    pub(crate) fn new_mont_inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // The Montgomery Modular Inverse - Revisited

            // Phase 1
            let modulus = *self.field.modulus();
            let mut u = modulus;
            // v = a * 2^m
            let mut v = self.repr;
            let mut r = F::Repr::from(0);
            let mut s = F::Repr::from(1);
            let mut k = 0u64;

            let mut found = false;
            for _ in 0..E::NUM_LIMBS*128 {
                if v.is_zero() {
                    found = true;
                    break;
                }
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
            if !found {
                return None;
            }

            if r >= modulus {
                r.sub_noborrow(&modulus);
            }

            let mut tmp = modulus;
            tmp.sub_noborrow(&r);

            r = tmp;

            // r = (a ^ -1) 2*(k - m)

            println!("Phase 1: r = {}", r);

            // phase 2

            let mont_power = self.field.mont_power();
            let modulus_bits_ceil = self.field.modulus_bits();
            if modulus_bits_ceil <= k && k <= mont_power + modulus_bits_ceil {
                let mut r_by_r2 = Self {
                    field: self.field,
                    repr: r
                };
    
                let r2 = Self {
                    field: self.field,
                    repr: *self.field.mont_r2()
                };
    
                r_by_r2.mul_assign(&r2);

                r = r_by_r2.repr;

                if k < mont_power {
                    k += mont_power;
                }
            } else {
                return None;
            }

            if k > 2*mont_power {
                return None;
            }

            // now we need montgomery(!) multiplication by 2^(2m - k)
            let mut two_in_two_m_minus_k_repr = F::Repr::from(1);
            for _ in 0..(2*mont_power - k) {
                two_in_two_m_minus_k_repr.mul2();
                if two_in_two_m_minus_k_repr >= modulus {
                    two_in_two_m_minus_k_repr.sub_noborrow(&modulus);
                }
            }
            debug_assert!(two_in_two_m_minus_k_repr < modulus);

            {
                let mut r_by_two_in_two_m_minus_k = Self {
                    field: self.field,
                    repr: r
                };
    
                let two_in_two_m_minus_k = Self {
                    field: self.field,
                    repr: two_in_two_m_minus_k_repr
                };

                r_by_two_in_two_m_minus_k.mul_assign(&two_in_two_m_minus_k);

                r = r_by_two_in_two_m_minus_k.repr;
            }

            // now back into montgomery form
            // let el = Fp::from_repr(self.field, r);
            let el = Fp::from_raw_repr(self.field, r);
            if el.is_err() {
                println!("Representation is invalid");
                return None;
            }
            let el = el.expect("guaranteed to exist");

            Some(el)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::traits::FieldElement;
    use crate::field::U256Repr;
    use crate::fp::Fp;
    use num_bigint::BigUint;
    use num_traits::Num;

    #[test]
    fn test_mont_inverse() {
        use crate::field::new_field;
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        // this is 7 in BE form
        let mut be_repr = vec![0u8; 32];
        be_repr[31] = 7u8;
        let element = Fp::from_be_bytes(&field, &be_repr[..], false).unwrap();
        let inverse = element.eea_inverse().unwrap();
        let mont_inverse = element.mont_inverse().unwrap();
        assert_eq!(inverse, mont_inverse);
    }

    #[test]
    fn test_new_mont_inverse() {
        use crate::field::new_field;
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        // this is 7 in BE form
        let mut be_repr = vec![0u8; 32];
        be_repr[31] = 7u8;
        let element = Fp::from_be_bytes(&field, &be_repr[..], false).unwrap();
        let inverse = element.eea_inverse().unwrap();
        let mont_inverse = element.new_mont_inverse().unwrap();
        assert_eq!(inverse, mont_inverse);
    }

    #[test]
    fn test_random_mont_inverse() {
        use rand::thread_rng;
        use rand::RngCore;
        use crate::field::new_field;
        let field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let mut rng = thread_rng();
        let mut be_repr = vec![0u8; 32];
        for _ in 0..1000 {
            rng.fill_bytes(&mut be_repr[..]);
            be_repr[0] = be_repr[0] & 0x1f;
            let element = Fp::from_be_bytes(&field, &be_repr[..], false).unwrap();
            let inverse = element.eea_inverse().unwrap();
            let mont_inverse = element.mont_inverse().unwrap();
            assert_eq!(inverse, mont_inverse);
        }
    }

    #[test]
    fn test_small_number_of_limbs_inverse_0() {
        use crate::field::new_field;
        use crate::traits::ZeroAndOne;

        let field = new_field::<U256Repr>("f98889454fc3", 16).unwrap();
        let value = BigUint::from_str_radix("1fe45623da1", 16).unwrap();
        let value_be_bytes = value.to_bytes_be();
        let element = Fp::from_be_bytes(&field, &value_be_bytes[..], true).unwrap();
        let inverse = element.eea_inverse().expect("inverse must exist");
        println!("EEA inverse = {}", inverse);
        let mont_inverse = element.mont_inverse().expect("montgomery inverse must exist");
        // let mont_inverse = element.new_mont_inverse().expect("montgomery inverse must exist");
        println!("Montgomery form inverse = {}", mont_inverse);

        let one = Fp::one(&field);

        let mut may_be_one = element.clone();
        may_be_one.mul_assign(&mont_inverse);
        assert!(may_be_one == one, "montgomery inverse is not an inverse");

        let mut may_be_one = element.clone();
        may_be_one.mul_assign(&inverse);

        assert!(may_be_one == one, "eea inverse is not an inverse");
        assert_eq!(inverse, mont_inverse);
    }

    #[test]
    fn test_small_number_of_limbs_inverse_1() {
        use crate::field::new_field;
        use crate::traits::ZeroAndOne;

        let field = new_field::<U256Repr>("54872962777895", 10).unwrap();
        let value = BigUint::from_str_radix("54872962777893", 10).unwrap();
        let value_be_bytes = value.to_bytes_be();
        let element = Fp::from_be_bytes(&field, &value_be_bytes[..], true).unwrap();
        let inverse = element.eea_inverse().expect("inverse must exist");
        println!("EEA inverse = {}", inverse);
        let mont_inverse = element.mont_inverse().expect("montgomery inverse must exist");
        // let mont_inverse = element.new_mont_inverse().expect("montgomery inverse must exist");
        println!("Montgomery form inverse = {}", mont_inverse);

        let one = Fp::one(&field);

        let mut may_be_one = element.clone();
        may_be_one.mul_assign(&mont_inverse);
        assert!(may_be_one == one, "montgomery inverse is not an inverse");

        let mut may_be_one = element.clone();
        may_be_one.mul_assign(&inverse);

        assert!(may_be_one == one, "eea inverse is not an inverse");
        assert_eq!(inverse, mont_inverse);
    }

    
}