use crate::field::SizedPrimeField;
use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, BitIterator};
use crate::weierstrass::Group;
use crate::weierstrass::curve::{WeierstrassCurve, CurvePoint};
use crate::weierstrass::cubic_twist::{WeierstrassCurveTwist, TwistPoint};
use crate::extension_towers::fp3::{Fp3, Extension3};
use crate::extension_towers::fp6_as_2_over_3::{Fp6, Extension2Over3};
use crate::pairings::PairingEngine;

#[derive(Eq, PartialEq, Clone, Copy)]
pub enum TwistType {
    D,
    M
}

pub struct PreparedTwistPoint<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> {
    is_infinity: bool,
    pub ell_coeffs: Vec<(Fp3<'a, FE, F>, Fp3<'a, FE, F>, Fp3<'a, FE, F>)>
}

impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>> PreparedTwistPoint<'a, FE, F> {
    pub fn is_zero(&self) -> bool {
        self.is_infinity
    }
}

pub struct CPInstance6<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> {
    pub x: Vec<u64>,
    pub x_is_negative: bool,
    pub twist_type: TwistType,
    pub base_field: &'a F,
    pub curve: &'a WeierstrassCurve<'a, FE, F, GE, G>,
    pub curve_twist: &'a WeierstrassCurveTwist<'a, FE, F, GE, G>,
    pub twist: Fp3<'a, FE, F>,
    fp3_extension: &'a Extension3<'a, FE, F>,
    fp6_extension: &'a Extension2Over3<'a, FE, F>,
}

impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> CPInstance6<'a, FE, F, GE, G> {
    fn miller_loop<'b, I>(&self, i: I) -> Fp6<'a, FE, F>
    where 'a: 'b,
        I: IntoIterator<
            Item = &'b (&'b CurvePoint<'a, FE, F, GE, G>, 
                &'b TwistPoint<'a, FE, F, GE, G>)
        >
    {

    }

    fn ate_pairing_loop(
        &self, 
        point: &CurvePoint<'a, FE, F, GE, G>, 
        twist_point: &TwistPoint<'a, FE, F, GE, G> 
    ) -> Fp6<'a, FE, F> {
        debug_assert!(point.is_normalized());
        debug_assert!(twist_point.is_normalized());
        let px = point.x.clone();
        let py = point.y.clone();
        let qx = twist_point.x.clone();
        let qy = twist_point.y.clone();
        let mut py_twist_squared = self.twist.clone();
        py_twist_squared.square();
        py_twist_squared.mul_by_fp(&py);

        let mut old_rx;
        let mut old_ry;
        let mut rx = qx;
        let mut ry = qy;

        let mut f = Fp6::one(self.fp6_extension);

        // The for loop is executed for all bits (EXCEPT the MSB itself) of
        // sw6_param_p (skipping leading zeros) in MSB to LSB order
        let mut found_one = false;
        for bit in BitIterator::new(self.x) {
            if !found_one && bit {
                found_one = true;
                continue;
            } else if !found_one {
                continue;
            }

            old_rx = rx;
            old_ry = ry;

            let mut old_rx_square = old_rx.clone();
            old_rx_square.square();
            let mut old_rx_square_3 = old_rx_square.clone();
            old_rx_square_3.double();
            old_rx_square_3.add_assign(&old_rx_square);
            let mut old_rx_square_3_a = old_rx_square_3.clone();
            old_rx_square_3_a.add_assign(&self.curve_twist.a);
            let mut old_ry_double_inverse = old_ry.clone();
            old_ry_double_inverse.double();
            let old_ry_double_inverse = old_ry_double_inverse.inverse().unwrap();

            let mut gamma = old_rx_square_3_a.clone();
            gamma.mul_assign(&old_ry_double_inverse);

            let mut gamma_twist = gamma.clone();
            gamma_twist.mul_assign(&self.twist);

            let mut gamma_old_rx = gamma.clone();
            gamma_old_rx.mul_assign(&old_rx);

            let mut gamma_twist_px = gamma_twist.clone();
            gamma_twist_px.mul_by_fp(&px);

            let x = py_twist_squared;

            let mut y = gamma_old_rx.clone();
            y.sub_assign(&old_ry);
            y.sub_assign(&gamma_twist_px);

            let ell_rr_at_p = Fp6 {
                c0: x,
                c1: y,
                extension_field: self.fp6_extension
            };

            rx = gamma.square() - &old_rx.double();
            ry = gamma * &(old_rx - &rx) - &old_ry;
            f = f.square() * &ell_rr_at_p;

            if bit {
                old_rx = rx;
                old_ry = ry;

                let gamma = (old_ry - &qy) * &((old_rx - &qx).inverse().unwrap());
                let gamma_twist = gamma * &TWIST;
                let gamma_qx = gamma * &qx;
                let mut gamma_twist_px = gamma_twist;
                gamma_twist_px.mul_assign_by_fp(&px);

                let x = py_twist_squared;
                let y = gamma_qx - &qy - &gamma_twist_px;
                let ell_rq_at_p = Fq6::new(x, y);

                rx = gamma.square() - &old_rx - &qx;
                ry = gamma * &(old_rx - &rx) - &old_ry;
                f = f * &ell_rq_at_p;
            }
        }

        f
    }

    fn final_exponentiation(&self, f: &Fp6<'a, FE, F>) -> Option<Fp6<'a, FE, F>> {
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
    use crate::pairings::{frobenius_calculator_fp2, frobenius_calculator_fp6, frobenius_calculator_fp12};
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

        let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6(modulus.clone(), &extension_6).unwrap();

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

        let (coeffs_c1, coeffs_c2) = frobenius_calculator_fp6(modulus.clone(), &extension_6).unwrap();

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