use crate::field::SizedPrimeField;
use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, BitIterator};
use crate::weierstrass::Group;
use crate::weierstrass::curve::{WeierstrassCurve, CurvePoint};
use crate::weierstrass::twist::{WeierstrassCurveTwist, TwistPoint};
use crate::extension_towers::fp2::{Fp2, Extension2};
use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};
use crate::extension_towers::fp6_as_3_over_2::{Fp6, Extension3Over2};
use crate::pairings::PairingEngine;

#[derive(Eq, PartialEq, Clone, Copy)]
pub enum TwistType {
    D,
    M
}

pub struct PreparedTwistPoint<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> {
    is_infinity: bool,
    pub ell_coeffs: Vec<(Fp2<'a, FE, F>, Fp2<'a, FE, F>, Fp2<'a, FE, F>)>
}

impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> PreparedTwistPoint<'a, FE, F> {
    pub fn is_zero(&self) -> bool {
        self.is_infinity
    }
}

pub struct Bls12Instance<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> {
    pub x: Vec<u64>,
    pub x_is_negative: bool,
    pub twist_type: TwistType,
    pub base_field: &'a F,
    pub curve: &'a WeierstrassCurve<'a, FE, F, GE, G>,
    pub curve_twist: &'a WeierstrassCurveTwist<'a, FE, F, GE, G>,
    fp2_extension: &'a Extension2<'a, FE, F>,
    fp6_extension: &'a Extension3Over2<'a, FE, F>,
    fp12_extension: &'a Extension2Over3Over2<'a, FE, F>,
}

impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> Bls12Instance<'a, FE, F, GE, G> {
    fn ell(
        &self,
        f: &mut Fp12<'a, FE, F>,
        coeffs: &(Fp2<'a, FE, F>, Fp2<'a, FE, F>, Fp2<'a, FE, F>),
        p: & CurvePoint<'a, FE, F, GE, G>,
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
        r: &mut TwistPoint<'a, FE, F, GE, G>,
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
        r: &mut TwistPoint<'a, FE, F, GE, G>,
        q: &TwistPoint<'a, FE, F, GE, G>,
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

    pub fn prepare(&self, twist_point: & TwistPoint<'a, FE, F, GE, G>) -> PreparedTwistPoint<'a, FE, F> {
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
        let mut r = TwistPoint::point_from_xy(&self.curve_twist, twist_point.x.clone(), twist_point.y.clone());

        for i in BitIterator::new(&self.x).skip(1) {
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
            Item = &'b (&'b CurvePoint<'a, FE, F, GE, G>, 
                &'b TwistPoint<'a, FE, F, GE, G>)
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

        for i in BitIterator::new(&self.x).skip(1) {
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


impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> PairingEngine for Bls12Instance<'a, FE, F, GE, G> {
    type PairingResult = Fp12<'a, FE, F>;
    type G1 = CurvePoint<'a, FE, F, GE, G>;
    type G2 = TwistPoint<'a, FE, F, GE, G>;

    fn pair<'b>
        (&self, points: &'b [CurvePoint<'a, FE, F, GE, G>], twists: &'b [TwistPoint<'a, FE, F, GE, G>]) -> Option<Self::PairingResult> {
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
    use num_traits::FromPrimitive;
    use num_integer::Integer;
    use num_traits::Zero;
    use crate::field::{U384Repr, U256Repr, new_field};
    use crate::fp::Fp;
    use crate::traits::{FieldElement};
    use crate::extension_towers::fp2::{Fp2, Extension2};
    use crate::extension_towers::fp6_as_3_over_2::{Fp6, Extension3Over2};
    use crate::extension_towers::fp12_as_2_over3_over_2::{Fp12, Extension2Over3Over2};
    use num_traits::Num;
    use crate::pairings::{frobenius_calculator_fp2, frobenius_calculator_fp6_as_3_over_2, frobenius_calculator_fp12};
    use crate::weierstrass::{Group};
    use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
    use crate::weierstrass::twist::{TwistPoint, WeierstrassCurveTwist};
    use crate::pairings::{PairingEngine};
    use rust_test::Bencher;

    #[test]
    fn test_bls12_381_pairing() {
        let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let scalar_field = new_field::<U256Repr>("52435875175126190479447740508185965837690552500527637822603658699938581184513", 10).unwrap();
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let mut extension_2 = Extension2 {
            field: &base_field,
            non_residue: fp_non_residue,
            frobenius_coeffs_c1: [Fp::zero(&base_field), Fp::zero(&base_field)]
        };

        let coeffs = frobenius_calculator_fp2(&extension_2).unwrap();
        extension_2.frobenius_coeffs_c1 = coeffs;

        let one = Fp::one(&base_field);

        let mut fp2_non_residue = Fp2::zero(&extension_2);
        fp2_non_residue.c0 = one.clone();
        fp2_non_residue.c1 = one.clone();

        let f_c1 = [Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2)];

        let mut extension_6 = Extension3Over2 {
            non_residue: fp2_non_residue,
            field: &extension_2,
            frobenius_coeffs_c1: f_c1.clone(),
            frobenius_coeffs_c2: f_c1,
        };

        let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6_as_3_over_2(modulus.clone(), &extension_6).unwrap();

        extension_6.frobenius_coeffs_c1 = coeffs_c1;
        extension_6.frobenius_coeffs_c2 = coeffs_c2;

        let mut fp2_non_residue = Fp2::zero(&extension_2);

         let f_c1 = [Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2)];

        let mut extension_12 = Extension2Over3Over2 {
            non_residue: Fp6::zero(&extension_6),
            field: &extension_6,
            frobenius_coeffs_c1: f_c1,
        };


        let coeffs = frobenius_calculator_fp12(modulus, &extension_12).unwrap();
        extension_12.frobenius_coeffs_c1 = coeffs;

        let b_fp = Fp::from_repr(&base_field, U384Repr::from(4)).unwrap();
        let mut b_fp2 = Fp2::zero(&extension_2);
        b_fp2.c0 = b_fp.clone();
        b_fp2.c1 = b_fp.clone();

        let a_fp = Fp::zero(&base_field);
        let a_fp2 = Fp2::zero(&extension_2);

        let curve = WeierstrassCurve::new(&scalar_field, a_fp, b_fp);
        let twist = WeierstrassCurveTwist::new(&scalar_field, &extension_2, a_fp2, b_fp2);

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
        let q = TwistPoint::point_from_xy(&twist, q_x, q_y);

        // let x = BigUint::from_str_radix("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507", 10).unwrap();
        // println!("x = {}", x);
        // println!("x = {:x}", x);

        assert!(p.check_on_curve());
        assert!(q.check_on_curve());

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

        // let expected_c0_c0_c0 = BigUint::from_str_radix("2819105605953691245277803056322684086884703000473961065716485506033588504203831029066448642358042597501014294104502", 10).unwrap();
        
        let expected_c0_c0_c0 = BigUint::from_str_radix("1250ebd871fc0a92a7b2d83168d0d727272d441befa15c503dd8e90ce98db3e7b6d194f60839c508a84305aaca1789b6", 16).unwrap();
        
        
        // println!("Res = {}", pairing_result);
    }

    #[bench]
    fn bench_bls12_381_pairing(b: &mut Bencher) {
        let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let base_field = new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let scalar_field = new_field::<U256Repr>("52435875175126190479447740508185965837690552500527637822603658699938581184513", 10).unwrap();
        let mut fp_non_residue = Fp::one(&base_field);
        fp_non_residue.negate(); // non-residue is -1

        let mut extension_2 = Extension2 {
            field: &base_field,
            non_residue: fp_non_residue,
            frobenius_coeffs_c1: [Fp::zero(&base_field), Fp::zero(&base_field)]
        };

        let coeffs = frobenius_calculator_fp2(&extension_2).unwrap();
        extension_2.frobenius_coeffs_c1 = coeffs;

        let one = Fp::one(&base_field);

        let mut fp2_non_residue = Fp2::zero(&extension_2);
        fp2_non_residue.c0 = one.clone();
        fp2_non_residue.c1 = one.clone();

        let f_c1 = [Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2)];

        let mut extension_6 = Extension3Over2 {
            non_residue: fp2_non_residue,
            field: &extension_2,
            frobenius_coeffs_c1: f_c1.clone(),
            frobenius_coeffs_c2: f_c1,
        };

        let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6_as_3_over_2(modulus.clone(), &extension_6).unwrap();

        extension_6.frobenius_coeffs_c1 = coeffs_c1;
        extension_6.frobenius_coeffs_c2 = coeffs_c2;

        let mut fp2_non_residue = Fp2::zero(&extension_2);

         let f_c1 = [Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2),
                    Fp2::zero(&extension_2), Fp2::zero(&extension_2), Fp2::zero(&extension_2)];

        let mut extension_12 = Extension2Over3Over2 {
            non_residue: Fp6::zero(&extension_6),
            field: &extension_6,
            frobenius_coeffs_c1: f_c1,
        };


        let coeffs = frobenius_calculator_fp12(modulus, &extension_12).unwrap();
        extension_12.frobenius_coeffs_c1 = coeffs;

        let b_fp = Fp::from_repr(&base_field, U384Repr::from(4)).unwrap();
        let mut b_fp2 = Fp2::zero(&extension_2);
        b_fp2.c0 = b_fp.clone();
        b_fp2.c1 = b_fp.clone();

        let a_fp = Fp::zero(&base_field);
        let a_fp2 = Fp2::zero(&extension_2);

        let curve = WeierstrassCurve::new(&scalar_field, a_fp, b_fp);
        let twist = WeierstrassCurveTwist::new(&scalar_field, &extension_2, a_fp2, b_fp2);

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
        let q = TwistPoint::point_from_xy(&twist, q_x, q_y);

        // let x = BigUint::from_str_radix("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507", 10).unwrap();
        // println!("x = {}", x);
        // println!("x = {:x}", x);

        assert!(p.check_on_curve());
        assert!(q.check_on_curve());

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

        b.iter(|| {
            bls12_engine.pair(&[p.clone()], &[q.clone()]).unwrap();
        });
    }
}