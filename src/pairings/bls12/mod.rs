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
use crate::pairings::PairingEngine;
use crate::pairings::TwistType;

pub struct PreparedTwistPoint<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> {
    is_infinity: bool,
    pub ell_coeffs: Vec<(Fp2<'a, FE, F>, Fp2<'a, FE, F>, Fp2<'a, FE, F>)>
}

impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> PreparedTwistPoint<'a, FE, F> {
    pub fn is_zero(&self) -> bool {
        self.is_infinity
    }
}

pub struct Bls12Instance<
    'a, 
        FE: ElementRepr, 
        F: SizedPrimeField<Repr = FE>, 
        CB: CurveParameters<BaseFieldElement = Fp<'a, FE, F>>,
        CTW: CurveParameters<BaseFieldElement = Fp2<'a, FE, F>>
    > {
    pub(crate) x: Vec<u64>,
    pub(crate) x_is_negative: bool,
    pub(crate) twist_type: TwistType,
    pub(crate) base_field: &'a F,
    pub(crate) curve: &'a WeierstrassCurve<'a, CB>,
    pub(crate) curve_twist: &'a WeierstrassCurve<'a, CTW>,
    pub(crate) fp2_extension: &'a Extension2<'a, FE, F>,
    pub(crate) fp6_extension: &'a Extension3Over2<'a, FE, F>,
    pub(crate) fp12_extension: &'a Extension2Over3Over2<'a, FE, F>,
}

impl<
    'a, 
        FE: ElementRepr, 
        F: SizedPrimeField<Repr = FE>, 
        CB: CurveParameters<BaseFieldElement = Fp<'a, FE, F>>,
        CTW: CurveParameters<BaseFieldElement = Fp2<'a, FE, F>>
    > Bls12Instance<'a, FE, F, CB, CTW> {
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
        *f = f.cyclotomic_exp(&self.x);
        if self.x_is_negative {
            f.conjugate();
        }
    }

    fn doubling_step(
        &self,
        r: &mut CurvePoint<'a, CTW>,
        two_inv: &Fp<'a, FE, F>,
    ) -> (Fp2<'a, FE, F>, Fp2<'a, FE, F>, Fp2<'a, FE, F>) {
        // Use adapted formulas from ZEXE instead
        let mut a = r.x.clone();
        a.mul_assign(&r.y);
        a.mul_by_fp(two_inv);
        let mut b = r.y.clone();
        b.square();
        let mut c = r.z.clone();
        c.square();

        let mut e = self.curve_twist.b.clone();
        let mut t0 = c.clone();
        t0.double();
        t0.add_assign(&c);

        e.mul_assign(&t0);

        let mut f = e.clone();
        f.double();
        f.add_assign(&e);

        let mut g = b.clone();
        g.add_assign(&f);
        g.mul_by_fp(two_inv);

        let mut h = r.y.clone();
        h.add_assign(&r.z);
        h.square();

        let mut t1 = b.clone();
        t1.add_assign(&c);

        h.sub_assign(&t1);

        let mut i = e.clone();
        i.sub_assign(&b);

        let mut j = r.x.clone();
        j.square();

        let mut e_square = e.clone();
        e_square.square();

        r.x = b.clone();
        r.x.sub_assign(&f);
        r.x.mul_assign(&a);

        let mut e_square_by_3 = e_square.clone();
        e_square_by_3.double();
        e_square_by_3.add_assign(&e_square);

        r.y = g;
        r.y.square();
        r.y.sub_assign(&e_square_by_3);

        r.z = b.clone();
        r.z.mul_assign(&h);

        let mut j_by_three = j.clone();
        j_by_three.double();
        j_by_three.add_assign(&j);
        h.negate();

        match self.twist_type {
            TwistType::M => {
                (i, j_by_three, h)
            },
            TwistType::D => {
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
        let mut theta = q.y.clone();
        theta.mul_assign(&r.z);
        theta.negate();
        theta.add_assign(&r.y);

        let mut lambda = q.x.clone();
        lambda.mul_assign(&r.z);
        lambda.negate();
        lambda.add_assign(&r.x);

        let mut c = theta.clone();
        c.square();
        let mut d = lambda.clone();
        d.square();
        let mut e = lambda.clone();
        e.mul_assign(&d);
        let mut f = r.z.clone();
        f.mul_assign(&c);
        let mut g = r.x.clone();
        g.mul_assign(&d);

        let mut h = g.clone();
        h.double();
        h.negate();
        h.add_assign(&e);
        h.add_assign(&f);
        

        r.x = lambda.clone();
        r.x.mul_assign(&h);

        let mut t0 = g.clone();
        t0.sub_assign(&h);
        t0.mul_assign(&theta);

        r.y.mul_assign(&e);
        r.y.negate();
        r.y.add_assign(&t0);

        r.z.mul_assign(&e);

        let mut t1 = lambda.clone();
        t1.mul_assign(&q.y);
        
        let mut j = theta.clone();
        j.mul_assign(&q.x);
        j.sub_assign(&t1);

        theta.negate();
        match self.twist_type {
            TwistType::M => (j, theta, lambda),
            TwistType::D => (lambda, theta, j),
        }
    }

    pub fn prepare(&self, twist_point: & CurvePoint<'a, CTW>) -> PreparedTwistPoint<'a, FE, F> {
        debug_assert!(twist_point.is_normalized());

        let mut two_inv = Fp::one(self.base_field);
        two_inv.double();
        let two_inv = two_inv.inverse().expect("inverse of 2 is guaranteed to exist");

        if twist_point.is_zero() {
            return PreparedTwistPoint {
                ell_coeffs: vec![],
                is_infinity:   true,
            };
        }

        let mut ell_coeffs = vec![];
        let mut r = CurvePoint::<CTW>::point_from_xy(&self.curve_twist, twist_point.x.clone(), twist_point.y.clone());

        for i in MsbBitIterator::new(&self.x).skip(1) {
            ell_coeffs.push(self.doubling_step(&mut r, &two_inv));

            if i {
                ell_coeffs.push(self.addition_step(&mut r, &twist_point));
            }
        }

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

        for i in MsbBitIterator::new(&self.x).skip(1) {
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

        if self.x_is_negative {
            f.conjugate();
        }

        f
    }

    fn final_exponentiation(&self, f: &Fp12<'a, FE, F>) -> Option<Fp12<'a, FE, F>> {
        // Computing the final exponentation following
        // https://eprint.iacr.org/2016/130.pdf.
        // We don't use their "faster" formula because it is difficult to make
        // it work for curves with odd `P::X`.
        // Hence we implement the algorithm from Table 1 below.

        // f1 = r.conjugate() = f^(p^6)
        let mut f1 = f.clone();
        // f1.conjugate();
        f1.frobenius_map(6);

        match f.inverse() {
            Some(mut f2) => {
                // f2 = f^(-1);
                // r = f^(p^6 - 1)
                let mut r = f1.clone();
                r.mul_assign(&f2);

                // f2 = f^(p^6 - 1)
                f2 = r.clone();
                // r = f^((p^6 - 1)(p^2))
                r.frobenius_map(2);

                // r = f^((p^6 - 1)(p^2) + (p^6 - 1))
                // r = f^((p^6 - 1)(p^2 + 1))
                r.mul_assign(&f2);

                // Hard part of the final exponentation is below:
                // From https://eprint.iacr.org/2016/130.pdf, Table 1
                let mut y0 = r.clone();
                y0.cyclotomic_square();
                y0.conjugate();

                let mut y5 = r.clone();
                self.exp_by_x(&mut y5);

                let mut y1 = y5.clone();
                y1.cyclotomic_square();

                let mut y3 = y0.clone();
                y3.mul_assign(&y5);

                let mut y0 = y3.clone();
                self.exp_by_x(&mut y0);
            
                let mut y2 = y0.clone();
                self.exp_by_x(&mut y2);

                let mut y4 = y2.clone();
                self.exp_by_x(&mut y4);
                y4.mul_assign(&y1);

                let mut y1 = y4.clone();
                self.exp_by_x(&mut y1);

                y3.conjugate();
                y1.mul_assign(&y3);
                y1.mul_assign(&r);

                let mut y3 = r.clone();
                y3.conjugate();
                y0.mul_assign(&r);
                y0.frobenius_map(3);

                y4.mul_assign(&y3);
                y4.frobenius_map(1);
                
                y5.mul_assign(&y2);
                y5.frobenius_map(2);

                y5.mul_assign(&y0);
                y5.mul_assign(&y4);
                y5.mul_assign(&y1);

                Some(y5)
            },
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
    > PairingEngine for Bls12Instance<'a, FE, F, CB, CTW> {
    type PairingResult = Fp12<'a, FE, F>;
    type G1 = CurvePoint<'a, CB>;
    type G2 = CurvePoint<'a, CTW>;

    fn pair<'b>
        (&self, points: &'b [CurvePoint<'a, CB>], twists: &'b [CurvePoint<'a, CTW>]) -> Option<Self::PairingResult> {
            if points.len() != twists.len() {
                return None;
            }

            if !std::option_env!("GAS_METERING").is_some() {
                if points.len() == 0 || twists.len() == 0 {
                    return None;
                }
            }
            
            let mut pairs = Vec::with_capacity(points.len());
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
    use crate::field::{U384Repr, new_field};
    use crate::fp::Fp;
    use crate::traits::{FieldElement, ZeroAndOne};
    use crate::extension_towers::fp2::{Fp2, Extension2};
    use crate::extension_towers::fp6_as_3_over_2::{Fp6, Extension3Over2};
    use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};
    use num_traits::Num;
    use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
    use crate::weierstrass::{CurveParameters, CurveOverFpParameters, CurveOverFp2Parameters};
    use crate::pairings::{PairingEngine};
    use crate::field::{biguint_to_u64_vec};
    use crate::sliding_window_exp::WindowExpBase;

    #[test]
    fn test_bls12_381_pairing() {
        let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        // let scalar_field = new_field::<U256Repr>("52435875175126190479447740508185965837690552500527637822603658699938581184513", 10).unwrap();
        let group_order = BigUint::from_str_radix("52435875175126190479447740508185965837690552500527637822603658699938581184513", 10).unwrap();
        let group_order = biguint_to_u64_vec(group_order);
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let mut extension_2 = Extension2::new(fp_non_residue);
        extension_2.calculate_frobenius_coeffs(modulus.clone()).expect("must work");

        let one = Fp::one(&base_field);

        let mut fp2_non_residue = Fp2::zero(&extension_2);
        fp2_non_residue.c0 = one.clone();
        fp2_non_residue.c1 = one.clone();

        let exp_base = WindowExpBase::new(&fp2_non_residue, Fp2::one(&extension_2), 8, 7);

        let mut extension_6 = Extension3Over2::new(fp2_non_residue);
        extension_6.calculate_frobenius_coeffs(modulus.clone(), &exp_base).expect("must work");

        let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(&extension_6));
        extension_12.calculate_frobenius_coeffs(modulus.clone(), &exp_base).expect("must work");

        let b_fp = Fp::from_repr(&base_field, U384Repr::from(4)).unwrap();
        let mut b_fp2 = Fp2::zero(&extension_2);
        b_fp2.c0 = b_fp.clone();
        b_fp2.c1 = b_fp.clone();

        let a_fp = Fp::zero(&base_field);
        let a_fp2 = Fp2::zero(&extension_2);

        let fp_params = CurveOverFpParameters::new(&base_field);
        let fp2_params = CurveOverFp2Parameters::new(&extension_2);

        let curve = WeierstrassCurve::new(group_order.clone(), a_fp, b_fp, &fp_params).unwrap();
        let twist = WeierstrassCurve::new(group_order.clone(), a_fp2, b_fp2, &fp2_params).unwrap();

        let p_x = BigUint::from_str_radix("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507", 10).unwrap().to_bytes_be();
        let p_y = BigUint::from_str_radix("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569", 10).unwrap().to_bytes_be();

        let q_x_0 = BigUint::from_str_radix("352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160", 10).unwrap().to_bytes_be();
        let q_x_1 = BigUint::from_str_radix("3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758", 10).unwrap().to_bytes_be();
        let q_y_0 = BigUint::from_str_radix("1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905", 10).unwrap().to_bytes_be();
        let q_y_1 = BigUint::from_str_radix("927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582", 10).unwrap().to_bytes_be();

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

        // let x = BigUint::from_str_radix("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507", 10).unwrap();
        // println!("x = {}", x);
        // println!("x = {:x}", x);

        assert!(p.is_on_curve());
        assert!(q.is_on_curve());

        let bls12_engine = super::Bls12Instance {
            x: vec![0xd201000000010000],
            x_is_negative: true,
            twist_type: super::TwistType::M,
            base_field: &base_field,
            curve: &curve,
            curve_twist: &twist,
            fp2_extension: &extension_2,
            fp6_extension: &extension_6,
            fp12_extension: &extension_12,
        };

        let pairing_result = bls12_engine.pair(&[p], &[q]).unwrap();

        // let expected_c0_c0_c0 = BigUint::from_str_radix("1250ebd871fc0a92a7b2d83168d0d727272d441befa15c503dd8e90ce98db3e7b6d194f60839c508a84305aaca1789b6", 16).unwrap();
        
        assert!(format!("{}",pairing_result.c0.c0.c0) == "0x1250ebd871fc0a92a7b2d83168d0d727272d441befa15c503dd8e90ce98db3e7b6d194f60839c508a84305aaca1789b6");
        // println!("Res = {}", pairing_result);
    }

    #[test]
    fn test_bls12_377_pairing() {
        let modulus = BigUint::from_str_radix("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177", 10).unwrap();
        let base_field = new_field::<U384Repr>("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177", 10).unwrap();
        // let scalar_field = new_field::<U256Repr>("8444461749428370424248824938781546531375899335154063827935233455917409239041", 10).unwrap();
        let group_order = BigUint::from_str_radix("8444461749428370424248824938781546531375899335154063827935233455917409239041", 10).unwrap();
        let group_order = biguint_to_u64_vec(group_order);
        let fp_nonres_repr = U384Repr::from(5);
        let mut fp_non_residue = Fp::from_repr(&base_field, fp_nonres_repr).unwrap();
        fp_non_residue.negate();

        let mut extension_2 = Extension2::new(fp_non_residue);
        extension_2.calculate_frobenius_coeffs(modulus.clone()).expect("must work");

        let one = Fp::one(&base_field);

        // it's just 0 + u
        let mut fp2_non_residue = Fp2::zero(&extension_2);
        fp2_non_residue.c1 = one.clone();

        let exp_base = WindowExpBase::new(&fp2_non_residue, Fp2::one(&extension_2), 8, 7);

        let mut extension_6 = Extension3Over2::new(fp2_non_residue.clone());
        extension_6.calculate_frobenius_coeffs(modulus.clone(), &exp_base).expect("must work");

        let mut extension_12 = Extension2Over3Over2::new(Fp6::zero(&extension_6));
        extension_12.calculate_frobenius_coeffs(modulus.clone(), &exp_base).expect("must work");

        let b_fp = Fp::from_repr(&base_field, U384Repr::from(1)).unwrap();
        let mut b_fp2 = fp2_non_residue.clone().inverse().unwrap();
        b_fp2.mul_by_fp(&b_fp);

        let a_fp = Fp::zero(&base_field);
        let a_fp2 = Fp2::zero(&extension_2);

        let fp_params = CurveOverFpParameters::new(&base_field);
        let fp2_params = CurveOverFp2Parameters::new(&extension_2);

        let curve = WeierstrassCurve::new(group_order.clone(), a_fp, b_fp, &fp_params).unwrap();
        let twist = WeierstrassCurve::new(group_order.clone(), a_fp2, b_fp2, &fp2_params).unwrap();

        let p_x = BigUint::from_str_radix("008848defe740a67c8fc6225bf87ff5485951e2caa9d41bb188282c8bd37cb5cd5481512ffcd394eeab9b16eb21be9ef", 16).unwrap().to_bytes_be();
        let p_y = BigUint::from_str_radix("01914a69c5102eff1f674f5d30afeec4bd7fb348ca3e52d96d182ad44fb82305c2fe3d3634a9591afd82de55559c8ea6", 16).unwrap().to_bytes_be();

        let q_x_0 = BigUint::from_str_radix("018480be71c785fec89630a2a3841d01c565f071203e50317ea501f557db6b9b71889f52bb53540274e3e48f7c005196", 16).unwrap().to_bytes_be();
        let q_x_1 = BigUint::from_str_radix("00ea6040e700403170dc5a51b1b140d5532777ee6651cecbe7223ece0799c9de5cf89984bff76fe6b26bfefa6ea16afe", 16).unwrap().to_bytes_be();
        let q_y_0 = BigUint::from_str_radix("00690d665d446f7bd960736bcbb2efb4de03ed7274b49a58e458c282f832d204f2cf88886d8c7c2ef094094409fd4ddf", 16).unwrap().to_bytes_be();
        let q_y_1 = BigUint::from_str_radix("00f8169fd28355189e549da3151a70aa61ef11ac3d591bf12463b01acee304c24279b83f5e52270bd9a1cdd185eb8f93", 16).unwrap().to_bytes_be();

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

        assert!(p.is_on_curve());
        assert!(q.is_on_curve());

        let bls12_engine = super::Bls12Instance {
            x: vec![0x8508c00000000001],
            x_is_negative: false,
            twist_type: super::TwistType::D,
            base_field: &base_field,
            curve: &curve,
            curve_twist: &twist,
            fp2_extension: &extension_2,
            fp6_extension: &extension_6,
            fp12_extension: &extension_12,
        };

        let pairing_result = bls12_engine.pair(&[p], &[q]).unwrap();

        assert!(format!("{}",pairing_result.c0.c0.c0) == "0x00b718ff624a95f189bfb44bcd6d6556226837c1f74d1afbf4bea573b71c17d3a243cae41d966e2164aad0991fd790cc");
    }
}