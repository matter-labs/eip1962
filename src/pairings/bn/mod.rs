use crate::field::SizedPrimeField;
use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, MsbBitIterator, ZeroAndOne};
use crate::weierstrass::Group;
use crate::weierstrass::{CurveParameters};
use crate::weierstrass::curve::{WeierstrassCurve, CurvePoint};
use crate::extension_towers::fp2::{Fp2, Extension2};
use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};
use crate::extension_towers::fp6_as_3_over_2::{Extension3Over2};
use crate::pairings::{PairingEngine};
use crate::pairings::TwistType;

pub struct PreparedTwistPoint<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> {
    is_infinity: bool,
    pub ell_coeffs: Vec<(Fp2<'a, FE, F>, Fp2<'a, FE, F>, Fp2<'a, FE, F>)>
}

impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> PreparedTwistPoint<'a, FE, F> {
    pub(crate) fn is_zero(&self) -> bool {
        self.is_infinity
    }
}

pub struct BnInstance<
    'a, 
        FE: ElementRepr, 
        F: SizedPrimeField<Repr = FE>, 
        CB: CurveParameters<BaseFieldElement = Fp<'a, FE, F>>,
        CTW: CurveParameters<BaseFieldElement = Fp2<'a, FE, F>>
    > {
    pub(crate) u: Vec<u64>,
    pub(crate) six_u_plus_2: Vec<u64>,
    pub(crate) u_is_negative: bool,
    pub(crate) twist_type: TwistType,
    pub(crate) base_field: &'a F,
    pub(crate) curve: &'a WeierstrassCurve<'a, CB>,
    pub(crate) curve_twist: &'a WeierstrassCurve<'a, CTW>,
    fp2_extension: &'a Extension2<'a, FE, F>,
    fp6_extension: &'a Extension3Over2<'a, FE, F>,
    fp12_extension: &'a Extension2Over3Over2<'a, FE, F>,
    non_residue_in_p_minus_one_over_2: Fp2<'a, FE, F>
}

impl<
    'a, 
        FE: ElementRepr, 
        F: SizedPrimeField<Repr = FE>, 
        CB: CurveParameters<BaseFieldElement = Fp<'a, FE, F>>,
        CTW: CurveParameters<BaseFieldElement = Fp2<'a, FE, F>>
    > BnInstance<'a, FE, F, CB, CTW> {
    fn ell(
        &self,
        f: &mut Fp12<'a, FE, F>,
        coeffs: &(Fp2<'a, FE, F>, Fp2<'a, FE, F>, Fp2<'a, FE, F>),
        p: & CurvePoint<'a, CB>,
    ) {
        debug_assert!(p.is_normalized());
        let mut c0 = coeffs.0.clone();
        let mut c1 = coeffs.1.clone();
        let mut c2 = coeffs.2.clone();

        match self.twist_type {
            TwistType::M => {
                c2.mul_by_fp(&p.y);
                c1.mul_by_fp(&p.x);
                f.mul_by_014(&c0, &c1, &c2);
            },
            TwistType::D => {
                c0.mul_by_fp(&p.y);
                c1.mul_by_fp(&p.x);
                f.mul_by_034(&c0, &c1, &c2);
            },
        }
    }

    fn exp_by_x(&self, f: &mut Fp12<'a, FE, F>) {
        *f = f.cyclotomic_exp(&self.u);
        if self.u_is_negative {
            f.conjugate();
        }
    }

    fn doubling_step(
        &self,
        r: &mut CurvePoint<'a, CTW>,
        two_inv: &Fp<'a, FE, F>,
    ) -> (Fp2<'a, FE, F>, Fp2<'a, FE, F>, Fp2<'a, FE, F>) {
        // Use adapted formulas from ZEXE instead

        // X*Y/2
        let mut a = r.x.clone();
        a.mul_assign(&r.y);
        a.mul_by_fp(two_inv);
        
        // Y^2
        let mut b = r.y.clone();
        b.square();

        // Z^2
        let mut c = r.z.clone();
        c.square();

        let mut e = self.curve_twist.b.clone();

        // 3*Z^2
        let mut t0 = c.clone();
        t0.double();
        t0.add_assign(&c);

        // 3*b*Z^2
        e.mul_assign(&t0);

        // 9*b*Z^2
        let mut f = e.clone();
        f.double();
        f.add_assign(&e);

        // (Y^2 + 9*b*Z^2)/2
        let mut g = b.clone();
        g.add_assign(&f);
        g.mul_by_fp(two_inv);
        
        // (Y + Z)^2
        let mut h = r.y.clone();
        h.add_assign(&r.z);
        h.square();

        // (Y^2 + Z^2)
        let mut t1 = b.clone();
        t1.add_assign(&c);

        // 2*Y*Z
        h.sub_assign(&t1);

        // 3*b*Z^2 - Y^2
        let mut i = e.clone();
        i.sub_assign(&b);

        // X^2
        let mut j = r.x.clone();
        j.square();

        // (3*b*Z^2)^2
        let mut e_square = e.clone();
        e_square.square();

        // X = (Y^2 - 9*b*Z^2)*X*Y/2
        r.x = b.clone();
        r.x.sub_assign(&f);
        r.x.mul_assign(&a);

        // 27*b^2*Z^4
        let mut e_square_by_3 = e_square.clone();
        e_square_by_3.double();
        e_square_by_3.add_assign(&e_square);

        // Y = ((Y^2 + 9*b*Z^2)/2)^2 - 27*b^2*Z^4
        r.y = g;
        r.y.square();
        r.y.sub_assign(&e_square_by_3);

        // Z = 2*Y^3*Z
        r.z = b.clone();
        r.z.mul_assign(&h);

        // 3*X^2
        let mut j_by_three = j.clone();
        j_by_three.double();
        j_by_three.add_assign(&j);

        // - 2*Y*Z
        h.negate();

        // i.mul_by_nonresidue(self.fp6_extension);
        match self.twist_type {
            TwistType::M => {
                (i, j_by_three, h)
            },
            TwistType::D => {
                // (0, 3, 4) = (-2*Y*Z, 3*X^2, 3*b*Z^2 - Y^2)
                (h, j_by_three, i)
            },
        }
    }

    fn addition_step(
        &self,
        r: &mut CurvePoint<'a, CTW>,
        q: & CurvePoint<'a, CTW>,
    ) -> (Fp2<'a, FE, F>, Fp2<'a, FE, F>, Fp2<'a, FE, F>) {
        debug_assert!(q.is_normalized());
        // use adapted zexe formulas too instead of ones from pairing crate
        // Capitals are coors of R (homogenious), normals are coordinates of Q (affine)
        // Y - y*Z
        let mut theta = q.y.clone();
        theta.mul_assign(&r.z);
        theta.negate();
        theta.add_assign(&r.y);

        // X - x*Z
        let mut lambda = q.x.clone();
        lambda.mul_assign(&r.z);
        lambda.negate();
        lambda.add_assign(&r.x);

        // Theta^2
        let mut c = theta.clone();
        c.square();
        
        // Lambda^2
        let mut d = lambda.clone();
        d.square();

        // Lambda^3
        let mut e = lambda.clone();
        e.mul_assign(&d);

        // Theta^2 * Z
        let mut f = r.z.clone();
        f.mul_assign(&c);

        // Lambda^2 * X
        let mut g = r.x.clone();
        g.mul_assign(&d);

        // Lambda^3 + Theta^2 * Z - 2*Lambda^2 * X
        let mut h = g.clone();
        h.double();
        h.negate();
        h.add_assign(&e);
        h.add_assign(&f);
        
        r.x = lambda.clone();
        r.x.mul_assign(&h);

        // (Lambda^2 * X - H)*Theta
        let mut t0 = g.clone();
        t0.sub_assign(&h);
        t0.mul_assign(&theta);

        // Y = (Lambda^2 * X - H)*Theta - Lambda^3 * Y
        r.y.mul_assign(&e);
        r.y.negate();
        r.y.add_assign(&t0);

        // Z = Lambda^3 * Z
        r.z.mul_assign(&e);

        // Lambda*y
        let mut t1 = lambda.clone();
        t1.mul_assign(&q.y);
        
        // Theta*x - Lambda*y
        let mut j = theta.clone();
        j.mul_assign(&q.x);
        j.sub_assign(&t1);

        theta.negate();

        // lambda.negate();
        // j.mul_by_nonresidue(self.fp6_extension);
        match self.twist_type {
            TwistType::M => (j, theta, lambda),
            TwistType::D => 
            {
                // (0, 3, 4) = (lambda, -theta, Theta*x - Lambda*y)
                (lambda, theta, j)
            },
        }
    }

    pub fn prepare(&self, twist_point: & CurvePoint<'a, CTW>) -> PreparedTwistPoint<'a, FE, F> {
        debug_assert!(twist_point.is_normalized());

        let mut two_inv = Fp::one(self.base_field);
        two_inv.double();
        let two_inv = two_inv.inverse().unwrap();

        if twist_point.is_zero() {
            return PreparedTwistPoint {
                ell_coeffs: vec![],
                is_infinity:   true,
            };
        }

        let mut ell_coeffs = vec![];
        let mut r = CurvePoint::point_from_xy(&self.curve_twist, twist_point.x.clone(), twist_point.y.clone());

        for i in MsbBitIterator::new(&self.six_u_plus_2).skip(1) {
            ell_coeffs.push(self.doubling_step(&mut r, &two_inv));

            if i {
                ell_coeffs.push(self.addition_step(&mut r, &twist_point));
            }
        }

        if self.u_is_negative {
            r.negate();
        }

        let mut q = twist_point.clone();

        q.x.c1.negate();
        q.x.mul_assign(&self.fp6_extension.frobenius_coeffs_c1[1]);

        q.y.c1.negate();
        q.y.mul_assign(&self.non_residue_in_p_minus_one_over_2);

        ell_coeffs.push(self.addition_step(&mut r, &q));

        let mut minusq2 = twist_point.clone();
        minusq2.x.mul_assign(&self.fp6_extension.frobenius_coeffs_c1[2]);

        ell_coeffs.push(self.addition_step(&mut r, &minusq2));

        PreparedTwistPoint {
            ell_coeffs,
            is_infinity: false,
        }
    }

    fn miller_loop<'b, I>(&self, i: I) -> Fp12<'a, FE, F>
    where 'a: 'b,
        I: IntoIterator<
            Item = &'b (&'b CurvePoint<'a, CB>, 
                &'b CurvePoint<'a, CTW>)
        >
    {
        let mut g1_references = vec![];
        let mut prepared_coeffs = vec![];

        for (p, q) in i.into_iter() {
            if !p.is_zero() && !q.is_zero() {
                let coeffs = self.prepare(&q.clone());
                let ell_coeffs = coeffs.ell_coeffs;
                prepared_coeffs.push(ell_coeffs);
                g1_references.push(p);
            }
        }

        let mut prepared_coeffs: Vec<_> = prepared_coeffs.into_iter().map(|el| el.into_iter()).collect();

        let mut f = Fp12::one(self.fp12_extension);

        for i in MsbBitIterator::new(&self.six_u_plus_2).skip(1) {    
            f.square();

            for (p, coeffs) in g1_references.iter().zip(prepared_coeffs.iter_mut()) {
                self.ell(&mut f, &coeffs.next().unwrap(), p);
            }

            if i {
                for (p, coeffs) in g1_references.iter().zip(prepared_coeffs.iter_mut()) {
                    self.ell(&mut f, &coeffs.next().unwrap(), p);
                }
            }
        }

        if self.u_is_negative {
            f.conjugate();
        }

        for (p, coeffs) in g1_references.iter().zip(prepared_coeffs.iter_mut()) {
            self.ell(&mut f, &coeffs.next().unwrap(), p);
        }

        for (p, coeffs) in g1_references.iter().zip(prepared_coeffs.iter_mut()) {
            self.ell(&mut f, &coeffs.next().unwrap(), p);
        }

        for (_p, coeffs) in g1_references.iter().zip(prepared_coeffs.iter_mut()) {
            debug_assert!(coeffs.next().is_none());
        }

        f
    }

    pub fn final_exponentiation(&self, f: &Fp12<'a, FE, F>) -> Option<Fp12<'a, FE, F>> {
        // use Zexe and pairing crate fused
        // https://eprint.iacr.org/2012/232.pdf

        // f1 = r.conjugate() = f^(p^6)
        let mut f1 = f.clone();
        f1.frobenius_map(6);

        match f.inverse() {
            Some(mut f2) => {
                let mut r = f1.clone();
                r.mul_assign(&f2);

                f2 = r.clone();

                r.frobenius_map(2);
                r.mul_assign(&f2);

                let mut fp = r.clone();
                fp.frobenius_map(1);

	            let mut fp2 = r.clone();
                fp2.frobenius_map(2);
                let mut fp3 = fp2.clone();
                fp3.frobenius_map(1);

                let mut fu = r.clone();
                self.exp_by_x(&mut fu);
                // exp_by_x(&mut fu, x);

                let mut fu2 = fu.clone();
                self.exp_by_x(&mut fu2);
                // exp_by_x(&mut fu2, x);

                let mut fu3 = fu2.clone();
                self.exp_by_x(&mut fu3);
                // exp_by_x(&mut fu3, x);

                let mut y3 = fu.clone();
                y3.frobenius_map(1);

                let mut fu2p = fu2.clone();
                fu2p.frobenius_map(1);

                let mut fu3p = fu3.clone();
                fu3p.frobenius_map(1);

                let mut y2 = fu2.clone();
                y2.frobenius_map(2);

                let mut y0 = fp.clone();
                y0.mul_assign(&fp2);
                y0.mul_assign(&fp3);

                let mut y1 = r.clone();
                y1.conjugate();

                let mut y5 = fu2.clone();
                y5.conjugate();

                y3.conjugate();

                let mut y4 = fu.clone();
                y4.mul_assign(&fu2p);
                y4.conjugate();

                let mut y6 = fu3.clone();
                y6.mul_assign(&fu3p);
                y6.conjugate();

                y6.square();
                y6.mul_assign(&y4);
                y6.mul_assign(&y5);

                let mut t1 = y3.clone();
                t1.mul_assign(&y5);
                t1.mul_assign(&y6);

                y6.mul_assign(&y2);

                t1.square();
                t1.mul_assign(&y6);
                t1.square();

                let mut t0 = t1.clone();
                t0.mul_assign(&y1);

                t1.mul_assign(&y0);

                t0.square();
                t0.mul_assign(&t1);

                Some(t0)
            }
            None => None,
        }
    }
}


impl<
    'a, 
        FE: ElementRepr, 
        F: SizedPrimeField<Repr = FE>, 
        CB: CurveParameters<BaseFieldElement = Fp<'a, FE, F>>,
        CTW: CurveParameters<BaseFieldElement = Fp2<'a, FE, F>>
    > PairingEngine for BnInstance<'a, FE, F, CB, CTW> {
    type PairingResult = Fp12<'a, FE, F>;
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
    use crate::field::{U256Repr, new_field};
    use crate::fp::Fp;
    use crate::traits::{FieldElement, ZeroAndOne};
    use crate::extension_towers::fp2::{Fp2, Extension2};
    use crate::extension_towers::fp6_as_3_over_2::{Fp6, Extension3Over2};
    use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};
    use num_traits::Num;
    use crate::pairings::{frobenius_calculator_fp2, frobenius_calculator_fp6_as_3_over_2, frobenius_calculator_fp12};
    use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
    use crate::weierstrass::{CurveParameters, CurveOverFpParameters, CurveOverFp2Parameters};
    use crate::pairings::{PairingEngine};
    use crate::representation::ElementRepr;
    use crate::field::{biguint_to_u64_vec};

    #[test]
    fn test_bn254_pairing() {
        let modulus = BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let base_field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        // let scalar_field = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let group_order = BigUint::from_str_radix("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let group_order = biguint_to_u64_vec(group_order);
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let mut extension_2 = Extension2::new(fp_non_residue);
        extension_2.calculate_frobenius_coeffs(modulus.clone());

        let one = Fp::one(&base_field);

        // non-residue is u+9
        let mut fp2_non_residue = Fp2::zero(&extension_2);
        let fp_9_repr = U256Repr::from(9u64);
        let fp_9 = Fp::from_repr(&base_field, fp_9_repr).unwrap(); 
        fp2_non_residue.c0 = fp_9.clone();
        fp2_non_residue.c1 = one.clone();

        let mut extension_6 = Extension3Over2::new(fp2_non_residue.clone());

        let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6_as_3_over_2(modulus.clone(), &extension_6).unwrap();

        extension_6.frobenius_coeffs_c1 = coeffs_c1;
        extension_6.frobenius_coeffs_c2 = coeffs_c2;
        extension_6.frobenius_coeffs_are_calculated = true;

        let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(&extension_6));

        let coeffs = frobenius_calculator_fp12(modulus.clone(), &extension_12).unwrap();
        extension_12.frobenius_coeffs_c1 = coeffs;
        extension_12.frobenius_coeffs_are_calculated = true;

        let b_fp = Fp::from_repr(&base_field, U256Repr::from(3)).unwrap();
        // here it's b/(u+9)
        let mut b_fp2 = fp2_non_residue.inverse().unwrap();
        b_fp2.mul_by_fp(&b_fp);

        let a_fp = Fp::zero(&base_field);
        let a_fp2 = Fp2::zero(&extension_2);

        let fp_params = CurveOverFpParameters::new(&base_field);
        let fp2_params = CurveOverFp2Parameters::new(&extension_2);

        let curve = WeierstrassCurve::new(group_order.clone(), a_fp, b_fp, &fp_params);
        let twist = WeierstrassCurve::new(group_order.clone(), a_fp2, b_fp2, &fp2_params);

        let p_x = BigUint::from_str_radix("1", 10).unwrap().to_bytes_be();
        let p_y = BigUint::from_str_radix("2", 10).unwrap().to_bytes_be();

        let q_x_0 = BigUint::from_str_radix("10857046999023057135944570762232829481370756359578518086990519993285655852781", 10).unwrap().to_bytes_be();
        let q_x_1 = BigUint::from_str_radix("11559732032986387107991004021392285783925812861821192530917403151452391805634", 10).unwrap().to_bytes_be();
        let q_y_0 = BigUint::from_str_radix("8495653923123431417604973247489272438418190587263600148770280649306958101930", 10).unwrap().to_bytes_be();
        let q_y_1 = BigUint::from_str_radix("4082367875863433681332203403145435568316851327593401208105741076214120093531", 10).unwrap().to_bytes_be();

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
        // println!("P.x = {}", p.x.into_repr());
        let q = CurvePoint::point_from_xy(&twist, q_x, q_y);
        // println!("Q.x = {}", q.x.c0.repr);
        // println!("Q.y = {}", q.y.c0.repr);

        // let x = BigUint::from_str_radix("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507", 10).unwrap();
        // println!("x = {}", x);
        // println!("x = {:x}", x);

        assert!(p.check_on_curve());
        assert!(q.check_on_curve());

        let mut minus_one_over_2 = Fp::one(&base_field);
        minus_one_over_2.negate();
        let mut two = Fp::one(&base_field);
        two.double();
        let two_inv = two.inverse().unwrap();
        minus_one_over_2.mul_assign(&two_inv);

        let non_residue_in_p_minus_one_over_2 = fp2_non_residue.pow(&minus_one_over_2.into_repr());

        let u = U256Repr::from(4965661367192848881);
        let mut six_u_plus_2 = u;
        six_u_plus_2.mul2();
        let two_u = six_u_plus_2;
        six_u_plus_2.mul2();
        six_u_plus_2.add_nocarry(&two_u);
        let mut two = U256Repr::from(1);
        two.mul2();
        six_u_plus_2.add_nocarry(&two);

        // println!("Expected coeff = {:x}", BigUint::from_str_radix("827617134098165717451808940080463277390770457691666780560712143809003953598", 10).unwrap());
        // println!("Expected coeff = {:x}", BigUint::from_str_radix("987776078024262725561041258416387561158070255475504730561661362421251696401", 10).unwrap());
        // println!("Expected coeff = {:x}", BigUint::from_str_radix("2813988028633040066320201189843971639620433430176492766961373503539074898364", 10).unwrap());

        let engine = super::BnInstance {
            u: vec![4965661367192848881],
            u_is_negative: false,
            six_u_plus_2: six_u_plus_2.0[..2].to_vec(),
            twist_type: super::TwistType::D,
            base_field: &base_field,
            curve: &curve,
            curve_twist: &twist,
            fp2_extension: &extension_2,
            fp6_extension: &extension_6,
            fp12_extension: &extension_12,
            non_residue_in_p_minus_one_over_2: non_residue_in_p_minus_one_over_2
        };

        let pairing_result = engine.pair(&[p], &[q]).unwrap();
        assert!(format!("{}", pairing_result.c0.c0.c0) == "0x12c70e90e12b7874510cd1707e8856f71bf7f61d72631e268fca81000db9a1f5");
    }

    #[test]
    fn test_final_exponentiation() {
        let modulus = BigUint::from_str_radix("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        let base_field = new_field::<U256Repr>("21888242871839275222246405745257275088696311157297823662689037894645226208583", 10).unwrap();
        // let scalar_field = new_field::<U256Repr>("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let group_order = BigUint::from_str_radix("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let group_order = biguint_to_u64_vec(group_order);
        let r = BigUint::from_str_radix("21888242871839275222246405745257275088548364400416034343698204186575808495617", 10).unwrap();
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let mut extension_2 = Extension2::new(fp_non_residue);
        extension_2.calculate_frobenius_coeffs(modulus.clone());

        let one = Fp::one(&base_field);

        // non-residue is u+9
        let mut fp2_non_residue = Fp2::zero(&extension_2);
        let fp_9_repr = U256Repr::from(9u64);
        let fp_9 = Fp::from_repr(&base_field, fp_9_repr).unwrap(); 
        fp2_non_residue.c0 = fp_9.clone();
        fp2_non_residue.c1 = one.clone();

        let mut extension_6 = Extension3Over2::new(fp2_non_residue.clone());

        let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6_as_3_over_2(modulus.clone(), &extension_6).unwrap();

        extension_6.frobenius_coeffs_c1 = coeffs_c1;
        extension_6.frobenius_coeffs_c2 = coeffs_c2;
        extension_6.frobenius_coeffs_are_calculated = true;

        let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(&extension_6));

        let coeffs = frobenius_calculator_fp12(modulus.clone(), &extension_12).unwrap();
        extension_12.frobenius_coeffs_c1 = coeffs;
        extension_12.frobenius_coeffs_are_calculated = true;

        let b_fp = Fp::from_repr(&base_field, U256Repr::from(3)).unwrap();
        // here it's b/(u+9)
        let mut b_fp2 = fp2_non_residue.inverse().unwrap();
        b_fp2.mul_by_fp(&b_fp);

        let a_fp = Fp::zero(&base_field);
        let a_fp2 = Fp2::zero(&extension_2);

        let fp_params = CurveOverFpParameters::new(&base_field);
        let fp2_params = CurveOverFp2Parameters::new(&extension_2);

        let curve = WeierstrassCurve::new(group_order.clone(), a_fp, b_fp, &fp_params);
        let twist = WeierstrassCurve::new(group_order.clone(), a_fp2, b_fp2, &fp2_params);

        let p_x = BigUint::from_str_radix("1", 10).unwrap().to_bytes_be();
        let p_y = BigUint::from_str_radix("2", 10).unwrap().to_bytes_be();

        let q_x_0 = BigUint::from_str_radix("10857046999023057135944570762232829481370756359578518086990519993285655852781", 10).unwrap().to_bytes_be();
        let q_x_1 = BigUint::from_str_radix("11559732032986387107991004021392285783925812861821192530917403151452391805634", 10).unwrap().to_bytes_be();
        let q_y_0 = BigUint::from_str_radix("8495653923123431417604973247489272438418190587263600148770280649306958101930", 10).unwrap().to_bytes_be();
        let q_y_1 = BigUint::from_str_radix("4082367875863433681332203403145435568316851327593401208105741076214120093531", 10).unwrap().to_bytes_be();

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
        let q = CurvePoint::point_from_xy(&twist, q_x, q_y);

        assert!(p.check_on_curve());
        assert!(q.check_on_curve());

        let mut minus_one_over_2 = Fp::one(&base_field);
        minus_one_over_2.negate();
        let mut two = Fp::one(&base_field);
        two.double();
        let two_inv = two.inverse().unwrap();
        minus_one_over_2.mul_assign(&two_inv);

        let non_residue_in_p_minus_one_over_2 = fp2_non_residue.pow(&minus_one_over_2.into_repr());

        let u = U256Repr::from(4965661367192848881);
        let mut six_u_plus_2 = u;
        six_u_plus_2.mul2();
        let two_u = six_u_plus_2;
        six_u_plus_2.mul2();
        six_u_plus_2.add_nocarry(&two_u);
        let mut two = U256Repr::from(1);
        two.mul2();
        six_u_plus_2.add_nocarry(&two);

        let engine = super::BnInstance {
            u: vec![4965661367192848881],
            u_is_negative: false,
            six_u_plus_2: six_u_plus_2.0[..2].to_vec(),
            twist_type: super::TwistType::D,
            base_field: &base_field,
            curve: &curve,
            curve_twist: &twist,
            fp2_extension: &extension_2,
            fp6_extension: &extension_6,
            fp12_extension: &extension_12,
            non_residue_in_p_minus_one_over_2: non_residue_in_p_minus_one_over_2
        };

        let mut fp12 = Fp12::zero(&extension_12);
        fp12.c0.c0 = fp2_non_residue.clone();
        fp12.c0.c1 = fp2_non_residue.clone();
        fp12.c0.c2 = fp2_non_residue.clone();
        fp12.c1.c0 = fp2_non_residue.clone();
        fp12.c1.c1 = fp2_non_residue.clone();
        fp12.c1.c2 = fp2_non_residue.clone();

        let final_exp = engine.final_exponentiation(&fp12).unwrap();

        use crate::field::biguint_to_u64_vec;
        let one = BigUint::from(1u64);
        let mut power = BigUint::from(1u64);
        for _ in 0..12 {
            power = power * modulus.clone();
        }
        power = power - one;
        power = power / r;
        let naive = fp12.pow(&biguint_to_u64_vec(power.clone()));

        assert!(naive == final_exp);
    }
}