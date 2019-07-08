use crate::field::SizedPrimeField;
use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, MsbBitIterator, ZeroAndOne};
use crate::weierstrass::{CurveParameters};
use crate::weierstrass::curve::{WeierstrassCurve, CurvePoint};
use crate::extension_towers::fp2::{Fp2, Extension2};
use crate::extension_towers::fp4_as_2_over_2::{Fp4, Extension2Over2};
use crate::pairings::PairingEngine;

pub struct MNT4Instance<
    'a, 
        FE: ElementRepr, 
        F: SizedPrimeField<Repr = FE>, 
        CB: CurveParameters<BaseFieldElement = Fp<'a, FE, F>>,
        CTW: CurveParameters<BaseFieldElement = Fp2<'a, FE, F>>
    > {
    pub(crate) x: Vec<u64>,
    pub(crate) x_is_negative: bool,
    pub(crate) exp_w0: Vec<u64>,
    pub(crate) exp_w1: Vec<u64>,
    pub(crate) exp_w0_is_negative: bool,
    pub(crate) base_field: &'a F,
    pub(crate) curve: &'a WeierstrassCurve<'a, CB>,
    pub(crate) curve_twist: &'a WeierstrassCurve<'a, CTW>,
    pub(crate) twist: Fp2<'a, FE, F>,
    pub(crate) fp2_extension: &'a Extension2<'a, FE, F>,
    pub(crate) fp4_extension: &'a Extension2Over2<'a, FE, F>,
}

struct PrecomputedG1<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> {
    pub x: Fp<'a, FE, F>,
    pub y: Fp<'a, FE, F>,
    pub x_by_twist: Fp2<'a, FE, F>,
    pub y_by_twist: Fp2<'a, FE, F>,
}

struct PrecomputedG2<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> {
    pub x: Fp2<'a, FE, F>,
    pub y: Fp2<'a, FE, F>,
    pub x_over_twist: Fp2<'a, FE, F>,
    pub y_over_twist: Fp2<'a, FE, F>,
    pub double_coefficients: Vec<AteDoubleCoefficients<'a, FE, F>>,
    pub addition_coefficients: Vec<AteAdditionCoefficients<'a, FE, F>>,
}

struct AteDoubleCoefficients<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> {
    pub c_h:  Fp2<'a, FE, F>,
    pub c_4c: Fp2<'a, FE, F>,
    pub c_j:  Fp2<'a, FE, F>,
    pub c_l:  Fp2<'a, FE, F>,
}

struct AteAdditionCoefficients<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> {
    pub c_l1: Fp2<'a, FE, F>,
    pub c_rz: Fp2<'a, FE, F>,
}

struct ExtendedCoordinates<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> {
    pub x: Fp2<'a, FE, F>,
    pub y: Fp2<'a, FE, F>,
    pub z: Fp2<'a, FE, F>,
    pub t: Fp2<'a, FE, F>,
}

impl<
    'a, 
        FE: ElementRepr, 
        F: SizedPrimeField<Repr = FE>, 
        CB: CurveParameters<BaseFieldElement = Fp<'a, FE, F>>,
        CTW: CurveParameters<BaseFieldElement = Fp2<'a, FE, F>>
    > MNT4Instance<'a, FE, F, CB, CTW> {
    fn miller_loop<'b, I>(&self, i: I) -> Fp4<'a, FE, F>
    where 'a: 'b,
        I: IntoIterator<
            Item = &'b (&'b CurvePoint<'a, CB>, 
                &'b CurvePoint<'a, CTW>)
        >
    {
        let mut f = Fp4::one(self.fp4_extension);
        for (p, q) in i.into_iter() {
            f.mul_assign(&self.ate_pairing_loop(p, q));
        }

        f
    }

    fn precompute_g1(&self, g1_point: &CurvePoint<'a, CB>) -> PrecomputedG1<'a, FE, F> {
        // not asserting normalization, it will be asserted in the loop
        let mut x_twist = self.twist.clone();
        x_twist.mul_by_fp(&g1_point.x);

        let mut y_twist = self.twist.clone();
        y_twist.mul_by_fp(&g1_point.y);

        PrecomputedG1 {
            x: g1_point.x.clone(),
            y: g1_point.y.clone(),
            x_by_twist: x_twist,
            y_by_twist: y_twist,
        }
    }

    fn doubling_step(&self, r: &mut ExtendedCoordinates<'a, FE, F>) -> AteDoubleCoefficients<'a, FE, F> {
        let mut a = r.t.clone();
        a.square();
        let mut b = r.x.clone();
        b.square();
        let mut c = r.y.clone();
        c.square();
        let mut d = c.clone();
        d.square();

        let mut e = r.x.clone();
        e.add_assign(&c);
        e.square();
        e.sub_assign(&b);
        e.sub_assign(&d);

        let mut f = self.curve_twist.a.clone();
        f.mul_assign(&a);
        f.add_assign(&b);
        f.add_assign(&b);
        f.add_assign(&b);

        let mut g = f.clone();
        g.square();

        let mut d_eight = d.clone();
        d_eight.double();
        d_eight.double();
        d_eight.double();

        let mut t0 = e.clone();
        t0.double();
        t0.double();

        let mut x = g.clone();
        x.sub_assign(&t0);

        let mut y = e.clone();
        y.double();
        y.sub_assign(&x);
        y.mul_assign(&f);
        y.sub_assign(&d_eight);

        let mut t0 = r.z.clone();
        t0.square();

        let mut z = r.y.clone();
        z.add_assign(&r.z);
        z.square();
        z.sub_assign(&c);
        z.sub_assign(&t0);

        let mut t = z.clone();
        t.square();

        let mut c_h = z.clone();
        c_h.add_assign(&r.t);
        c_h.square();
        c_h.sub_assign(&t);
        c_h.sub_assign(&a);

        let mut c_4c = c.clone();
        c_4c.double();
        c_4c.double();

        let mut c_j = f.clone();
        c_j.add_assign(&r.t);
        c_j.square();
        c_j.sub_assign(&g);
        c_j.sub_assign(&a);

        let mut c_l = f.clone();
        c_l.add_assign(&r.x);
        c_l.square();
        c_l.sub_assign(&g);
        c_l.sub_assign(&b);

        let coeff = AteDoubleCoefficients { c_h, c_4c, c_j, c_l};

        r.x = x;
        r.y = y;
        r.z = z;
        r.t = t;

        coeff
    }

    fn addition_step(&self,         
        x: &Fp2<'a, FE, F>,
        y: &Fp2<'a, FE, F>,
        r: &mut ExtendedCoordinates<'a, FE, F>) 
    -> AteAdditionCoefficients<'a, FE, F> {
        let mut a = y.clone();
        a.square();
        let mut b = r.t.clone();
        b.mul_assign(&x);

        let mut d = r.z.clone();
        d.add_assign(&y);
        d.square();
        d.sub_assign(&a);
        d.sub_assign(&r.t);
        d.mul_assign(&r.t);

        let mut h = b.clone();
        h.sub_assign(&r.x);

        let mut i = h.clone();
        i.square();

        let mut e = i.clone();
        e.double();
        e.double();

        let mut j = h.clone();
        j.mul_assign(&e);

        let mut v = r.x.clone();
        v.mul_assign(&e);

        let mut l1 = d.clone();
        l1.sub_assign(&r.y);
        l1.sub_assign(&r.y);

        let mut x = l1.clone();
        x.square();
        x.sub_assign(&j);
        x.sub_assign(&v);
        x.sub_assign(&v);

        let mut t0 = r.y.clone();
        t0.double();
        t0.mul_assign(&j);

        let mut y = v.clone();
        y.sub_assign(&x);
        y.mul_assign(&l1);
        y.sub_assign(&t0);

        let mut z = r.z.clone();
        z.add_assign(&h);
        z.square();
        z.sub_assign(&r.t);
        z.sub_assign(&i);

        let mut t = z.clone();
        t.square();

        let coeff = AteAdditionCoefficients { c_l1: l1, c_rz: z.clone() };

        r.x = x;
        r.y = y;
        r.z = z;
        r.t = t;

        coeff
    }


    fn precompute_g2(&self, g2_point: &CurvePoint<'a, CTW>, twist_inv: &Fp2<'a, FE, F>) -> PrecomputedG2<'a, FE, F> {
        // not asserting normalization, it will be asserted in the loop
        // precompute addition and doubling coefficients
        let mut x_over_twist = g2_point.x.clone();
        x_over_twist.mul_assign(&twist_inv);

        let mut y_over_twist = g2_point.y.clone();
        y_over_twist.mul_assign(&twist_inv);

        let mut g2_p = PrecomputedG2 {
            x: g2_point.x.clone(),
            y: g2_point.y.clone(),
            x_over_twist: x_over_twist,
            y_over_twist: y_over_twist,
            double_coefficients: vec![],
            addition_coefficients: vec![],
        };

        let mut r = ExtendedCoordinates {
            x: g2_point.x.clone(),
            y: g2_point.y.clone(),
            z: Fp2::one(&self.fp2_extension),
            t: Fp2::one(&self.fp2_extension),
        };

        for bit in MsbBitIterator::new(&self.x).skip(1) {
            let coeff = self.doubling_step(&mut r);
            g2_p.double_coefficients.push(coeff);

            if bit {
                let coeff = self.addition_step(&g2_point.x, &g2_point.y, &mut r);
                g2_p.addition_coefficients.push(coeff);
            }
        }

        if self.x_is_negative {
            let rz_inv = r.z.inverse().unwrap();
            let mut rz2_inv = rz_inv.clone();
            rz2_inv.square();
            let mut rz3_inv = rz_inv.clone();
            rz3_inv.mul_assign(&rz2_inv);

            let mut minus_r_affine_x = rz2_inv;
            minus_r_affine_x.mul_assign(&r.x);
            let mut minus_r_affine_y = rz3_inv;
            minus_r_affine_y.mul_assign(&r.y);
            minus_r_affine_y.negate();

            let coeff = self.addition_step(
                &minus_r_affine_x,
                &minus_r_affine_y,
                &mut r,
            );

            g2_p.addition_coefficients.push(coeff);
        }

        g2_p
    }

    fn ate_pairing_loop(
        &self, 
        point: &CurvePoint<'a, CB>, 
        twist_point: &CurvePoint<'a, CTW> 
    ) -> Fp4<'a, FE, F> {
        debug_assert!(point.is_normalized());
        debug_assert!(twist_point.is_normalized());

        let twist_inv = self.twist.inverse().unwrap();

        let p = self.precompute_g1(&point);
        let q = self.precompute_g2(&twist_point, &twist_inv);
        let mut l1_coeff = Fp2::zero(&self.fp2_extension);
        l1_coeff.c0 = p.x.clone();
        l1_coeff.sub_assign(&q.x_over_twist);

        let mut f = Fp4::one(self.fp4_extension);

        let mut dbl_idx: usize = 0;
        let mut add_idx: usize = 0;

        // The for loop is executed for all bits (EXCEPT the MSB itself) of
        for bit in MsbBitIterator::new(&self.x).skip(1) {

            let dc = &q.double_coefficients[dbl_idx];
            dbl_idx += 1;

            let mut g_rr_at_p = Fp4::zero(&self.fp4_extension);

            let mut t0 = dc.c_j.clone();
            t0.mul_assign(&p.x_by_twist);
            t0.negate();
            t0.add_assign(&dc.c_l);
            t0.sub_assign(&dc.c_4c);

            let mut t1 = dc.c_h.clone();
            t1.mul_assign(&p.y_by_twist);

            g_rr_at_p.c0 = t0;
            g_rr_at_p.c1 = t1;

            f.square();
            f.mul_assign(&g_rr_at_p);

            if bit {
                let ac = &q.addition_coefficients[add_idx];
                add_idx += 1;

                let mut g_rq_at_p = Fp4::zero(&self.fp4_extension);

                let mut t0 = ac.c_rz.clone();
                t0.mul_assign(&p.y_by_twist);

                let mut t = l1_coeff.clone();
                t.mul_assign(&ac.c_l1);

                let mut t1 = q.y_over_twist.clone();
                t1.mul_assign(&ac.c_rz);
                t1.add_assign(&t);
                t1.negate();

                g_rq_at_p.c0 = t0;
                g_rq_at_p.c1 = t1;

                f.mul_assign(&g_rq_at_p);
            }
        }

        if self.x_is_negative {
            let ac = &q.addition_coefficients[add_idx];

            let mut g_rnegr_at_p = Fp4::zero(&self.fp4_extension);

            let mut t0 = ac.c_rz.clone();
            t0.mul_assign(&p.y_by_twist);

            let mut t = l1_coeff.clone();
            t.mul_assign(&ac.c_l1);

            let mut t1 = q.y_over_twist.clone();
            t1.mul_assign(&ac.c_rz);
            t1.add_assign(&t);
            t1.negate();

            g_rnegr_at_p.c0 = t0;
            g_rnegr_at_p.c1 = t1;

            f.mul_assign(&g_rnegr_at_p);
            f = f.inverse().unwrap();
        }

        f
    }

    fn final_exponentiation(&self, f: &Fp4<'a, FE, F>) -> Option<Fp4<'a, FE, F>> {
        let value_inv = f.inverse();
        if value_inv.is_none() {
            return None;
        }
        let value_inv = value_inv.unwrap();
        let value_to_first_chunk = self.final_exponentiation_part_one(f, &value_inv);
        let value_inv_to_first_chunk = self.final_exponentiation_part_one(&value_inv, f);
        
        Some(self.final_exponentiation_part_two(&value_to_first_chunk, &value_inv_to_first_chunk))
    }

    fn final_exponentiation_part_one(&self, elt: &Fp4<'a, FE, F>, elt_inv: &Fp4<'a, FE, F>) -> Fp4<'a, FE, F> {
        /* (q^2-1) */

        /* elt_q2 = elt^(q^2) */
        let mut elt_q2 = elt.clone();
        elt_q2.frobenius_map(2);
        /* elt_q2_over_elt = elt^(q^2-1) */
        let mut elt_q2_over_elt = elt_q2;
        elt_q2_over_elt.mul_assign(&elt_inv);

        elt_q2_over_elt
    }

    fn final_exponentiation_part_two(&self, elt: &Fp4<'a, FE, F>, elt_inv: &Fp4<'a, FE, F>) -> Fp4<'a, FE, F> {
        let mut elt_q = elt.clone();
        elt_q.frobenius_map(1);

        let mut w1_part = elt_q.cyclotomic_exp(&self.exp_w1);
        let w0_part = match self.exp_w0_is_negative {
            true => elt_inv.cyclotomic_exp(&self.exp_w0),
            false => elt.cyclotomic_exp(&self.exp_w0),
        };

        w1_part.mul_assign(&w0_part);

        w1_part
    }
}


impl<
    'a, 
        FE: ElementRepr, 
        F: SizedPrimeField<Repr = FE>, 
        CB: CurveParameters<BaseFieldElement = Fp<'a, FE, F>>,
        CTW: CurveParameters<BaseFieldElement = Fp2<'a, FE, F>>
    > PairingEngine for MNT4Instance<'a, FE, F, CB, CTW> {
    type PairingResult = Fp4<'a, FE, F>;
    type G1 = CurvePoint<'a, CB>;
    type G2 = CurvePoint<'a, CTW>;

    fn pair<'b>
        (&self, points: &'b [CurvePoint<'a, CB>], twists: &'b [CurvePoint<'a, CTW>]) -> Option<Self::PairingResult> {
            let mut pairs = vec![];
            for (p, q) in points.iter().zip(twists.iter()) {
                pairs.push((p, q));
            }
            let loop_result = self.miller_loop(&pairs[..]);

            self.final_exponentiation(&loop_result)
        }   
}

#[cfg(test)]
mod tests {
    use num_bigint::BigUint;
    use crate::field::{U320Repr, new_field, biguint_to_u64_vec};
    use crate::fp::Fp;
    use crate::traits::{FieldElement, ZeroAndOne};
    use crate::extension_towers::fp2::{Fp2, Extension2};
    use crate::extension_towers::fp4_as_2_over_2::{Extension2Over2};
    use num_traits::Num;
    use crate::pairings::{frobenius_calculator_fp2, frobenius_calculator_fp4_as_2_over_2};
    use crate::weierstrass::{Group, CurveParameters, CurveOverFpParameters, CurveOverFp2Parameters};
    use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
    use crate::pairings::{PairingEngine};

    #[test]
    fn test_mnt4_pairing() {
        let modulus = BigUint::from_str_radix("475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081", 10).unwrap();
        let base_field = new_field::<U320Repr>("475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081", 10).unwrap();
        let nonres_repr = U320Repr::from(17);
        let fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

        let mut extension_2 = Extension2::new(fp_non_residue.clone());
        extension_2.calculate_frobenius_coeffs(modulus.clone());

        let one = Fp::one(&base_field);

        let mut fp2_non_residue = Fp2::zero(&extension_2); // non-residue is 13 + 0*u + 0*u^2
        fp2_non_residue.c0 = fp_non_residue;

        let f_c1 = [Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field),
                    Fp::zero(&base_field)];

        let mut extension_4 = Extension2Over2::new(fp2_non_residue);

        let coeffs = frobenius_calculator_fp4_as_2_over_2(modulus, &extension_4).unwrap();
        extension_4.frobenius_coeffs_c1 = coeffs;
        extension_4.frobenius_coeffs_are_calculated = true;

        let b_fp = BigUint::from_str_radix("423894536526684178289416011533888240029318103673896002803341544124054745019340795360841685", 10).unwrap().to_bytes_be();
        let b_fp = Fp::from_be_bytes(&base_field, &b_fp, true).unwrap();

        let a_fp = Fp::from_repr(&base_field, U320Repr::from(2)).unwrap();

        let mut twist = Fp2::zero(&extension_2);
        twist.c1 = one.clone();

        let mut twist_squared = twist.clone();
        twist_squared.square();

        let mut twist_cubed = twist_squared.clone();
        twist_cubed.mul_assign(&twist);

        let mut a_fp2 = twist_squared.clone();
        a_fp2.mul_by_fp(&a_fp);

        let mut b_fp2 = twist_cubed.clone();
        b_fp2.mul_by_fp(&b_fp);

        // let scalar_field = new_field::<U320Repr>("475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137", 10).unwrap();
        let group_order = BigUint::from_str_radix("475922286169261325753349249653048451545124878552823515553267735739164647307408490559963137", 10).unwrap();
        let group_order = biguint_to_u64_vec(group_order);

        let fp_params = CurveOverFpParameters::new(&base_field);
        let fp2_params = CurveOverFp2Parameters::new(&extension_2);

        let curve = WeierstrassCurve::new(group_order.clone(), a_fp, b_fp, &fp_params);
        let curve_twist = WeierstrassCurve::new(group_order.clone(), a_fp2, b_fp2, &fp2_params);

        let p_x = BigUint::from_str_radix("60760244141852568949126569781626075788424196370144486719385562369396875346601926534016838", 10).unwrap().to_bytes_be();
        let p_y = BigUint::from_str_radix("363732850702582978263902770815145784459747722357071843971107674179038674942891694705904306", 10).unwrap().to_bytes_be();

        let q_x_0 = BigUint::from_str_radix("438374926219350099854919100077809681842783509163790991847867546339851681564223481322252708", 10).unwrap().to_bytes_be();
        let q_x_1 = BigUint::from_str_radix("37620953615500480110935514360923278605464476459712393277679280819942849043649216370485641", 10).unwrap().to_bytes_be();
  
        let q_y_0 = BigUint::from_str_radix("37437409008528968268352521034936931842973546441370663118543015118291998305624025037512482", 10).unwrap().to_bytes_be();
        let q_y_1 = BigUint::from_str_radix("424621479598893882672393190337420680597584695892317197646113820787463109735345923009077489", 10).unwrap().to_bytes_be();

        let p_x = Fp::from_be_bytes(&base_field, &p_x, true).unwrap();
        let p_y = Fp::from_be_bytes(&base_field, &p_y, true).unwrap();

        let q_x_0 = Fp::from_be_bytes(&base_field, &q_x_0, true).unwrap();
        let q_x_1 = Fp::from_be_bytes(&base_field, &q_x_1, true).unwrap();

        let q_y_0 = Fp::from_be_bytes(&base_field, &q_y_0, true).unwrap();
        let q_y_1 = Fp::from_be_bytes(&base_field, &q_y_1, true).unwrap();

        let mut q_x = Fp2::zero(&extension_2);
        q_x.c0 = q_x_0;
        q_x.c1 = q_x_1;

        let mut q_y = Fp2::zero(&extension_2);
        q_y.c0 = q_y_0;
        q_y.c1 = q_y_1;

        let p = CurvePoint::point_from_xy(&curve, p_x, p_y);
        let q = CurvePoint::point_from_xy(&curve_twist, q_x, q_y);

        assert!(p.check_on_curve());
        assert!(q.check_on_curve());

        let engine = super::MNT4Instance {
            x: biguint_to_u64_vec(BigUint::from_str_radix("689871209842287392837045615510547309923794944", 10).unwrap()),
            x_is_negative: false,
            exp_w0: biguint_to_u64_vec(BigUint::from_str_radix("689871209842287392837045615510547309923794945", 10).unwrap()),
            exp_w1: vec![1u64],
            exp_w0_is_negative: false,
            base_field: &base_field,
            curve: &curve,
            curve_twist: &curve_twist,
            twist: twist,
            fp2_extension: &extension_2,
            fp4_extension: &extension_4,
        };

        let mut p2 = p.mul(vec![12345678]);
        p2.normalize();

        let mut q2 = q.mul(vec![12345678]);
        q2.normalize();

        let _pairing_result = engine.pair(&[p.clone()], &[q.clone()]).unwrap();

        let ans1 = engine.pair(&[p.clone()], &[q2]).unwrap();
        let ans2 = engine.pair(&[p2], &[q.clone()]).unwrap();
        let ans3 = engine.pair(&[p], &[q]).unwrap();
        let ans3 = ans3.pow(&vec![12345678]);

        assert!(ans1 == ans2);
        assert!(ans1 == ans3);

        // // let expected_c0_c0_c0 = BigUint::from_str_radix("2819105605953691245277803056322684086884703000473961065716485506033588504203831029066448642358042597501014294104502", 10).unwrap();
        
        // let expected_c0_c0_c0 = BigUint::from_str_radix("1250ebd871fc0a92a7b2d83168d0d727272d441befa15c503dd8e90ce98db3e7b6d194f60839c508a84305aaca1789b6", 16).unwrap();
        
        // println!("Res = {}", pairing_result);

        // assert!(format!("{}", pairing_result.c0.c0) == "0x0000014ac12149eebffe74a1c75a7225deb91ca243c49eef01392080ff519ab6209431f81b50ec03");
    }

    #[test]
    fn test_mnt4_753_pairing() {
        use crate::field::U768Repr;

        let modulus = BigUint::from_str_radix("41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888253786114353726529584385201591605722013126468931404347949840543007986327743462853720628051692141265303114721689601", 10).unwrap();
        let base_field = new_field::<U768Repr>("41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888253786114353726529584385201591605722013126468931404347949840543007986327743462853720628051692141265303114721689601", 10).unwrap();
        let nonres_repr = U768Repr::from(13);
        let fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

        let mut extension_2 = Extension2::new(fp_non_residue.clone());
        extension_2.calculate_frobenius_coeffs(modulus.clone()).expect("must work");

        let one = Fp::one(&base_field);

        let mut fp2_non_residue = Fp2::zero(&extension_2); // non-residue is 13 + 0*u + 0*u^2
        fp2_non_residue.c0 = fp_non_residue;

        let mut extension_4 = Extension2Over2::new(fp2_non_residue);

        let coeffs = frobenius_calculator_fp4_as_2_over_2(modulus, &extension_4).unwrap();
        extension_4.frobenius_coeffs_c1 = coeffs;
        extension_4.frobenius_coeffs_are_calculated = true;

        let b_fp = BigUint::from_str_radix("28798803903456388891410036793299405764940372360099938340752576406393880372126970068421383312482853541572780087363938442377933706865252053507077543420534380486492786626556269083255657125025963825610840222568694137138741554679540", 10).unwrap().to_bytes_be();
        let b_fp = Fp::from_be_bytes(&base_field, &b_fp, true).unwrap();

        let a_fp = Fp::from_repr(&base_field, U768Repr::from(2)).unwrap();

        let mut twist = Fp2::zero(&extension_2);
        twist.c1 = one.clone();

        let mut twist_squared = twist.clone();
        twist_squared.square();

        let mut twist_cubed = twist_squared.clone();
        twist_cubed.mul_assign(&twist);

        let mut a_fp2 = twist_squared.clone();
        a_fp2.mul_by_fp(&a_fp);

        let mut b_fp2 = twist_cubed.clone();
        b_fp2.mul_by_fp(&b_fp);

        let group_order = BigUint::from_str_radix("41898490967918953402344214791240637128170709919953949071783502921025352812571106773058893763790338921418070971888458477323173057491593855069696241854796396165721416325350064441470418137846398469611935719059908164220784476160001", 10).unwrap();
        let group_order = biguint_to_u64_vec(group_order.clone());

        let fp_params = CurveOverFpParameters::new(&base_field);
        let fp2_params = CurveOverFp2Parameters::new(&extension_2);

        let curve = WeierstrassCurve::new(group_order.clone(), a_fp, b_fp, &fp_params);
        let curve_twist = WeierstrassCurve::new(group_order.clone(), a_fp2, b_fp2, &fp2_params);

        let p_x = BigUint::from_str_radix("23803503838482697364219212396100314255266282256287758532210460958670711284501374254909249084643549104668878996224193897061976788052185662569738774028756446662400954817676947337090686257134874703224133183061214213216866019444443", 10).unwrap().to_bytes_be();
        let p_y = BigUint::from_str_radix("21091012152938225813050540665280291929032924333518476279110711148670464794818544820522390295209715531901248676888544060590943737249563733104806697968779796610374994498702698840169538725164956072726942500665132927942037078135054", 10).unwrap().to_bytes_be();

        let q_x_0 = BigUint::from_str_radix("22367666623321080720060256844679369841450849258634485122226826668687008928557241162389052587294939105987791589807198701072089850184203060629036090027206884547397819080026926412256978135536735656049173059573120822105654153939204", 10).unwrap().to_bytes_be();
        let q_x_1 = BigUint::from_str_radix("19674349354065582663569886390557105215375764356464013910804136534831880915742161945711267871023918136941472003751075703860943205026648847064247080124670799190998395234694182621794580160576822167228187443851233972049521455293042", 10).unwrap().to_bytes_be();
  
        let q_y_0 = BigUint::from_str_radix("6945425020677398967988875731588951175743495235863391886533295045397037605326535330657361771765903175481062759367498970743022872494546449436815843306838794729313050998681159000579427733029709987073254733976366326071957733646574", 10).unwrap().to_bytes_be();
        let q_y_1 = BigUint::from_str_radix("17406100775489352738678485154027036191618283163679980195193677896785273172506466216232026037788788436442188057889820014276378772936042638717710384987239430912364681046070625200474931975266875995282055499803236813013874788622488", 10).unwrap().to_bytes_be();

        let p_x = Fp::from_be_bytes(&base_field, &p_x, true).unwrap();
        let p_y = Fp::from_be_bytes(&base_field, &p_y, true).unwrap();

        let q_x_0 = Fp::from_be_bytes(&base_field, &q_x_0, true).unwrap();
        let q_x_1 = Fp::from_be_bytes(&base_field, &q_x_1, true).unwrap();

        let q_y_0 = Fp::from_be_bytes(&base_field, &q_y_0, true).unwrap();
        let q_y_1 = Fp::from_be_bytes(&base_field, &q_y_1, true).unwrap();

        let mut q_x = Fp2::zero(&extension_2);
        q_x.c0 = q_x_0;
        q_x.c1 = q_x_1;

        let mut q_y = Fp2::zero(&extension_2);
        q_y.c0 = q_y_0;
        q_y.c1 = q_y_1;

        let p = CurvePoint::point_from_xy(&curve, p_x, p_y);
        let q = CurvePoint::point_from_xy(&curve_twist, q_x, q_y);

        let mut p0 = p.clone();
        p0.double();
        p0.negate();
        p0.normalize();

        let q0 = q.clone();

        let p1 = p.clone();
        let mut q1 = q.clone();
        q1.double();
        q1.normalize();

        assert!(p0.check_on_curve());
        assert!(p1.check_on_curve());
        assert!(q0.check_on_curve());
        assert!(q1.check_on_curve());

        let engine = super::MNT4Instance {
            x: biguint_to_u64_vec(BigUint::from_str_radix("204691208819330962009469868104636132783269696790011977400223898462431810102935615891307667367766898917669754470400", 10).unwrap()),
            x_is_negative: true,
            exp_w0: biguint_to_u64_vec(BigUint::from_str_radix("204691208819330962009469868104636132783269696790011977400223898462431810102935615891307667367766898917669754470399", 10).unwrap()),
            exp_w1: vec![1u64],
            exp_w0_is_negative: true,
            base_field: &base_field,
            curve: &curve,
            curve_twist: &curve_twist,
            twist: twist,
            fp2_extension: &extension_2,
            fp4_extension: &extension_4,
        };

        let mut p2 = p.mul(vec![12345678]);
        p2.normalize();

        let mut q2 = q.mul(vec![12345678]);
        q2.normalize();

        let pairing_result = engine.pair(&[p.clone()], &[q.clone()]).unwrap();

        let ans1 = engine.pair(&[p.clone()], &[q2]).unwrap();
        let ans2 = engine.pair(&[p2], &[q.clone()]).unwrap();
        let ans3 = engine.pair(&[p], &[q]).unwrap();
        let ans3 = ans3.pow(&vec![12345678]);

        assert!(ans1 == ans2);
        assert!(ans1 == ans3);
    }
}