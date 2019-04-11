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
        let mut c0 = coeffs.0.clone();
        let mut c1 = coeffs.1.clone();
        let mut c2 = coeffs.2.clone();

        let mut p = p.clone();
        p.normalize();

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
        f.cyclotomic_exp(&self.x);
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

        match self.twist_type {
            TwistType::M => {
                let mut j_by_three = j.clone();
                j_by_three.double();
                j_by_three.add_assign(&j);
                h.negate();

                (i, j_by_three, h)
            },
            TwistType::D => {
                let mut j_by_three = j.clone();
                j_by_three.double();
                j_by_three.add_assign(&j);
                h.negate();
                (h, j_by_three, i)
            },
        }
    }

    fn addition_step(
        &self,
        r: &mut TwistPoint<'a, FE, F, GE, G>,
        q: &TwistPoint<'a, FE, F, GE, G>,
    ) -> (Fp2<'a, FE, F>, Fp2<'a, FE, F>, Fp2<'a, FE, F>) {
        debug_assert!(r.is_normalized());
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
}


impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> PairingEngine for Bls12Instance<'a, FE, F, GE, G> {
    type PairingResult = Fp12<'a, FE, F>;
    type G1 = CurvePoint<'a, FE, F, GE, G>;
    type G2 = TwistPoint<'a, FE, F, GE, G>;

    fn pair<'b>
        (&self, point: &'b CurvePoint<'a, FE, F, GE, G>, twist_point: &'b TwistPoint<'a, FE, F, GE, G>) -> Self::PairingResult {
            Fp12::one(self.fp12_extension)
        }   
}