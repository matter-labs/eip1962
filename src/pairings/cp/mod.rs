use crate::field::SizedPrimeField;
use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::traits::{FieldElement, BitIterator, MsbBitIterator};
use crate::weierstrass::Group;
use crate::weierstrass::curve::{WeierstrassCurve, CurvePoint};
use crate::weierstrass::cubic_twist::{WeierstrassCurveTwist, TwistPoint};
use crate::extension_towers::fp3::{Fp3, Extension3};
use crate::extension_towers::fp6_as_2_over_3::{Fp6, Extension2Over3};
use crate::pairings::PairingEngine;

pub struct CPInstance6<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> {
    pub x: Vec<u64>,
    pub x_is_negative: bool,
    pub exp_w0: Vec<u64>,
    pub exp_w1: Vec<u64>,
    pub exp_w0_is_negative: bool,
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
        let mut f = Fp6::one(self.fp6_extension);
        for (p, q) in i.into_iter() {
            f.mul_assign(&self.ate_pairing_loop(p, q));
        }

        f
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
        let mut rx = qx.clone();
        let mut ry = qy.clone();

        let mut f = Fp6::one(self.fp6_extension);

        // The for loop is executed for all bits (EXCEPT the MSB itself) of
        // sw6_param_p (skipping leading zeros) in MSB to LSB order
        // let mut found_one = false;
        for bit in MsbBitIterator::new(&self.x).skip(1) {
        // for bit in BitIterator::new(&self.x) {
        //     if !found_one && bit {
        //         found_one = true;
        //         continue;
        //     } else if !found_one {
        //         continue;
        //     }

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

            let x = py_twist_squared.clone();

            let mut y = gamma_old_rx.clone();
            y.sub_assign(&old_ry);
            y.sub_assign(&gamma_twist_px);

            let ell_rr_at_p = Fp6 {
                c0: x,
                c1: y,
                extension_field: self.fp6_extension
            };

            rx = gamma.clone();
            rx.square();
            let mut t0 = old_rx.clone();
            t0.double();
            rx.sub_assign(&t0);

            let mut t0 = old_rx.clone();
            t0.sub_assign(&rx);

            ry = gamma.clone();
            ry.mul_assign(&t0);
            ry.sub_assign(&old_ry);

            f.square();
            f.mul_assign(&ell_rr_at_p);

            if bit {
                old_rx = rx.clone();
                old_ry = ry.clone();

                let mut t0 = old_ry.clone();
                t0.sub_assign(&qy);

                let mut t1 = old_rx.clone();
                t1.sub_assign(&qx);
                let t1 = t1.inverse().unwrap();

                let mut gamma = t0;
                gamma.mul_assign(&t1);
                let mut gamma_twist = gamma.clone();
                gamma_twist.mul_assign(&self.twist);
                let mut gamma_qx = gamma.clone();
                gamma_qx.mul_assign(&qx);
                let mut gamma_twist_px = gamma_twist.clone();
                gamma_twist_px.mul_by_fp(&px);

                let x = py_twist_squared.clone();
                let mut y = gamma_qx.clone();
                y.sub_assign(&qy);
                y.sub_assign(&gamma_twist_px);

                let ell_rq_at_p = Fp6 {
                    c0: x,
                    c1: y,
                    extension_field: self.fp6_extension
                };

                rx = gamma.clone();
                rx.square();
                rx.sub_assign(&old_rx);
                rx.sub_assign(&qx);

                ry = old_rx.clone();
                ry.sub_assign(&rx);
                ry.mul_assign(&gamma);
                ry.sub_assign(&old_ry);

                f.mul_assign(&ell_rq_at_p);
            }
        }

        f
    }

    fn final_exponentiation(&self, f: &Fp6<'a, FE, F>) -> Option<Fp6<'a, FE, F>> {
        let value_inv = f.inverse();
        if value_inv.is_none() {
            return None;
        }
        let value_inv = value_inv.unwrap();
        let value_to_first_chunk = self.final_exponentiation_part_one(f, &value_inv);
        let value_inv_to_first_chunk = self.final_exponentiation_part_one(&value_inv, f);
        
        Some(self.final_exponentiation_part_two(&value_to_first_chunk, &value_inv_to_first_chunk))
    }

    fn final_exponentiation_part_one(&self, elt: &Fp6<'a, FE, F>, elt_inv: &Fp6<'a, FE, F>) -> Fp6<'a, FE, F> {
        // (q^3-1)*(q+1)

        // elt_q3 = elt^(q^3)
        let mut elt_q3 = elt.clone();
        elt_q3.frobenius_map(3);
        // elt_q3_over_elt = elt^(q^3-1)
        let mut elt_q3_over_elt = elt_q3;
        elt_q3_over_elt.mul_assign(&elt_inv);
        // alpha = elt^((q^3-1) * q)
        let mut alpha = elt_q3_over_elt.clone();
        alpha.frobenius_map(1);
        // beta = elt^((q^3-1)*(q+1)
        alpha.mul_assign(&elt_q3_over_elt);

        alpha
    }

    fn final_exponentiation_part_two(&self, elt: &Fp6<'a, FE, F>, elt_inv: &Fp6<'a, FE, F>) -> Fp6<'a, FE, F> {
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


impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> PairingEngine for CPInstance6<'a, FE, F, GE, G> {
    type PairingResult = Fp6<'a, FE, F>;
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
    use crate::field::{U832Repr, new_field, biguint_to_u64_vec};
    use crate::fp::Fp;
    use crate::traits::{FieldElement};
    use crate::extension_towers::fp3::{Fp3, Extension3};
    use crate::extension_towers::fp6_as_2_over_3::{Fp6, Extension2Over3};
    use num_traits::Num;
    use crate::pairings::{frobenius_calculator_fp3, frobenius_calculator_fp6_as_2_over_3};
    use crate::weierstrass::{Group};
    use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
    use crate::weierstrass::cubic_twist::{TwistPoint, WeierstrassCurveTwist};
    use crate::pairings::{PairingEngine};
    use rust_test::Bencher;

    #[test]
    fn test_cp6_pairing() {
        let modulus = BigUint::from_str_radix("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
        let base_field = new_field::<U832Repr>("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
        let nonres_repr = U832Repr::from(13);
        let mut fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

        let mut extension_3 = Extension3 {
            field: &base_field,
            non_residue: fp_non_residue.clone(),
            frobenius_coeffs_c1: [Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field)],
            frobenius_coeffs_c2: [Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field)]
        };

        let (coeffs_1, coeffs_2) = frobenius_calculator_fp3(modulus.clone(), &extension_3).unwrap();
        extension_3.frobenius_coeffs_c1 = coeffs_1;
        extension_3.frobenius_coeffs_c2 = coeffs_2;

        let one = Fp::one(&base_field);

        let mut fp3_non_residue = Fp3::zero(&extension_3); // non-residue is 13 + 0*u + 0*u^2
        fp3_non_residue.c0 = fp_non_residue;

        let f_c1 = [Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field),
                    Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field)];

        let mut extension_6 = Extension2Over3 {
            non_residue: fp3_non_residue,
            field: &extension_3,
            frobenius_coeffs_c1: f_c1
        };

        let [c0, c1, c2, c3, c4, c5] = frobenius_calculator_fp6_as_2_over_3(modulus, &extension_6).unwrap();
        extension_6.frobenius_coeffs_c1 = [c0.c0, c1.c0, c2.c0, c3.c0, c4.c0, c5.c0];

        let b_fp = BigUint::from_str_radix("17764315118651679038286329069295091506801468118146712649886336045535808055361274148466772191243305528312843236347777260247138934336850548243151534538734724191505953341403463040067571652261229308333392040104884438208594329793895206056414", 10).unwrap().to_bytes_be();
        let b_fp = Fp::from_be_bytes(&base_field, &b_fp, true).unwrap();

        let a_fp = Fp::from_repr(&base_field, U832Repr::from(5)).unwrap();

        let mut twist = Fp3::zero(&extension_3);
        twist.c1 = one.clone();

        let mut twist_squared = twist.clone();
        twist_squared.square();

        let mut twist_cubed = twist_squared.clone();
        twist_cubed.mul_assign(&twist);

        let mut a_fp3 = twist_squared.clone();
        a_fp3.mul_by_fp(&a_fp);

        let mut b_fp3 = Fp3::zero(&extension_3);
        b_fp3.c0 = one.clone();
        b_fp3.c0.mul_assign(&b_fp);

        let mut b_fp3 = twist_cubed.clone();
        b_fp3.mul_by_fp(&b_fp);

        let scalar_field = new_field::<U832Repr>("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();

        let curve = WeierstrassCurve::new(&scalar_field, a_fp, b_fp);
        let curve_twist = WeierstrassCurveTwist::new(&scalar_field, &extension_3, a_fp3, b_fp3);

        let p_x = BigUint::from_str_radix("5511163824921585887915590525772884263960974614921003940645351443740084257508990841338974915037175497689287870585840954231884082785026301437744745393958283053278991955159266640440849940136976927372133743626748847559939620888818486853646", 10).unwrap().to_bytes_be();
        let p_y = BigUint::from_str_radix("7913123550914612057135582061699117755797758113868200992327595317370485234417808273674357776714522052694559358668442301647906991623400754234679697332299689255516547752391831738454121261248793568285885897998257357202903170202349380518443", 10).unwrap().to_bytes_be();

        let q_x_0 = BigUint::from_str_radix("13426761183630949215425595811885033211332897733228446437546263564078445562454176776915160094418980045665397361295624472103734543457352048745726512354895954850428989867542989474136256025045975283415690491751906307188562464175510373683338", 10).unwrap().to_bytes_be();
        let q_x_1 = BigUint::from_str_radix("20471601555918880743198170952645906008198510944268658573129351735028343217532386920456705632337352161031960990613816401042894531220068552819818037605513359562118363589199569321421558696125646867661360498323171027455638052943806292028610", 10).unwrap().to_bytes_be();
        let q_x_2 = BigUint::from_str_radix("3905053196875761830053608605277158152930144841844497593936739534395003062685449846381431331169369910535935138116320442345524758217411779027270883193856999691582831339845600938304719916501940381093815781408183227875600753651697934495980", 10).unwrap().to_bytes_be();
        
        let q_y_0 = BigUint::from_str_radix("8567517639523571619872938228644013584947463594196306323477160496987712111576624702939472765993995586889532559039169098780892505598589581147768095093536988446010255611523736706017580686335404469207486594272103717837888228343074699140243", 10).unwrap().to_bytes_be();
        let q_y_1 = BigUint::from_str_radix("3890537069205870914984502594450293167889863914413852788876350245583932846980126025043974070704295857226211547108005650399870458089721518559480870503159804530091559886149680718531004778697982910253701559194337987238111062202037698927752", 10).unwrap().to_bytes_be();
        let q_y_2 = BigUint::from_str_radix("10936269922612615564271188303104593362724754284143779051599749016735041389483971486958818324356025479751246744831831158558101688599198721653921723013062333636402617118847009085485166284126970598561393411916461254016145116183331671450721", 10).unwrap().to_bytes_be();

        let p_x = Fp::from_be_bytes(&base_field, &p_x, true).unwrap();
        let p_y = Fp::from_be_bytes(&base_field, &p_y, true).unwrap();

        let q_x_0 = Fp::from_be_bytes(&base_field, &q_x_0, true).unwrap();
        let q_x_1 = Fp::from_be_bytes(&base_field, &q_x_1, true).unwrap();
        let q_x_2 = Fp::from_be_bytes(&base_field, &q_x_2, true).unwrap();

        let q_y_0 = Fp::from_be_bytes(&base_field, &q_y_0, true).unwrap();
        let q_y_1 = Fp::from_be_bytes(&base_field, &q_y_1, true).unwrap();
        let q_y_2 = Fp::from_be_bytes(&base_field, &q_y_2, true).unwrap();

        let mut q_x = Fp3::zero(&extension_3);
        q_x.c0 = q_x_0;
        q_x.c1 = q_x_1;
        q_x.c2 = q_x_2;

        let mut q_y = Fp3::zero(&extension_3);
        q_y.c0 = q_y_0;
        q_y.c1 = q_y_1;
        q_y.c2 = q_y_2;

        let p = CurvePoint::point_from_xy(&curve, p_x, p_y);
        let q = TwistPoint::point_from_xy(&curve_twist, q_x, q_y);

        // let x = BigUint::from_str_radix("506464946133393486072777102926336625944849939610982267859828541006717966526573193706126370441346337661774335955699621", 10).unwrap();
        // println!("X len = {}", biguint_to_u64_vec(x.clone()).len());
        // println!("{:x}", biguint_to_u64_vec(x.clone())[0]);
        let w0 = BigUint::from_str_radix("7000705447348627246181409558336018323010329260726930841638672011287206690002601216854775649561085256265269640040570922609783227469279331691880282815325569032149343779036142830666859805506518426649197067288711084398033", 10).unwrap();
        let w1 = BigUint::from_str_radix("86482221941698704497288378992285180119495364068003923046442785886272123124361700722982503222189455144364945735564951562986", 10).unwrap();
        
        let x: Vec<u64> = vec![
            0x55c5b9b57b942ae8,
            0x3d52287d3dfd424a,
            0xcf1ff9d6a543deb7,
            0x820c9c5711ceeebc,
            0x549a2d44305d20fe,
            0x50f5c131afd70235,
            0xab3596c8617c5792,
            0x830c728d80f9d78b,
            0x6a7223ee72023d07,
            0xbc5d176b746af026,
            0xe959283d8f526663,
            0xc4d2263babf8941f,
            0x3848,
        ];

        assert!(p.check_on_curve());
        assert!(q.check_on_curve());

        let engine = super::CPInstance6 {
            x: x,
            x_is_negative: false,
            exp_w0: biguint_to_u64_vec(w0),
            exp_w1: biguint_to_u64_vec(w1),
            exp_w0_is_negative: true,
            base_field: &base_field,
            curve: &curve,
            curve_twist: &curve_twist,
            twist: twist,
            fp3_extension: &extension_3,
            fp6_extension: &extension_6,
        };

        let pairing_result = engine.pair(&[p], &[q]).unwrap();

        // // let expected_c0_c0_c0 = BigUint::from_str_radix("2819105605953691245277803056322684086884703000473961065716485506033588504203831029066448642358042597501014294104502", 10).unwrap();
        
        // let expected_c0_c0_c0 = BigUint::from_str_radix("1250ebd871fc0a92a7b2d83168d0d727272d441befa15c503dd8e90ce98db3e7b6d194f60839c508a84305aaca1789b6", 16).unwrap();
        
        // println!("Res = {}", pairing_result);

        assert!(format!("{}", pairing_result.c0.c0) == "0x0000000000003621057ef03d9de232637cd9d965a98d2670da76793d38546678fb8f5f148ddc5b1cbcd74c66f8d8462197406e42e9713aa158c0a4b02e4f26c785f73a0b9027c2ff50282278d2e91afbfa6dfa584f55987f960ce39228cf56ab169ba8932adc28df");
    }

    #[bench]
    fn bench_cp6_pairing(b: &mut Bencher) {
        let modulus = BigUint::from_str_radix("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
        let base_field = new_field::<U832Repr>("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();
        let nonres_repr = U832Repr::from(13);
        let mut fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

        let mut extension_3 = Extension3 {
            field: &base_field,
            non_residue: fp_non_residue.clone(),
            frobenius_coeffs_c1: [Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field)],
            frobenius_coeffs_c2: [Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field)]
        };

        let (coeffs_1, coeffs_2) = frobenius_calculator_fp3(modulus.clone(), &extension_3).unwrap();
        extension_3.frobenius_coeffs_c1 = coeffs_1;
        extension_3.frobenius_coeffs_c2 = coeffs_2;

        let one = Fp::one(&base_field);

        let mut fp3_non_residue = Fp3::zero(&extension_3); // non-residue is 13 + 0*u + 0*u^2
        fp3_non_residue.c0 = fp_non_residue;

        let f_c1 = [Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field),
                    Fp::zero(&base_field), Fp::zero(&base_field), Fp::zero(&base_field)];

        let mut extension_6 = Extension2Over3 {
            non_residue: fp3_non_residue,
            field: &extension_3,
            frobenius_coeffs_c1: f_c1
        };

        let [c0, c1, c2, c3, c4, c5] = frobenius_calculator_fp6_as_2_over_3(modulus, &extension_6).unwrap();
        extension_6.frobenius_coeffs_c1 = [c0.c0, c1.c0, c2.c0, c3.c0, c4.c0, c5.c0];

        let b_fp = BigUint::from_str_radix("17764315118651679038286329069295091506801468118146712649886336045535808055361274148466772191243305528312843236347777260247138934336850548243151534538734724191505953341403463040067571652261229308333392040104884438208594329793895206056414", 10).unwrap().to_bytes_be();
        let b_fp = Fp::from_be_bytes(&base_field, &b_fp, true).unwrap();

        let a_fp = Fp::from_repr(&base_field, U832Repr::from(5)).unwrap();

        let mut twist = Fp3::zero(&extension_3);
        twist.c1 = one.clone();

        let mut twist_squared = twist.clone();
        twist_squared.square();

        let mut twist_cubed = twist_squared.clone();
        twist_cubed.mul_assign(&twist);

        let mut a_fp3 = twist_squared.clone();
        a_fp3.mul_by_fp(&a_fp);

        let mut b_fp3 = Fp3::zero(&extension_3);
        b_fp3.c0 = one.clone();
        b_fp3.c0.mul_assign(&b_fp);

        let mut b_fp3 = twist_cubed.clone();
        b_fp3.mul_by_fp(&b_fp);

        let scalar_field = new_field::<U832Repr>("22369874298875696930346742206501054934775599465297184582183496627646774052458024540232479018147881220178054575403841904557897715222633333372134756426301062487682326574958588001132586331462553235407484089304633076250782629492557320825577", 10).unwrap();

        let curve = WeierstrassCurve::new(&scalar_field, a_fp, b_fp);
        let curve_twist = WeierstrassCurveTwist::new(&scalar_field, &extension_3, a_fp3, b_fp3);

        let p_x = BigUint::from_str_radix("5511163824921585887915590525772884263960974614921003940645351443740084257508990841338974915037175497689287870585840954231884082785026301437744745393958283053278991955159266640440849940136976927372133743626748847559939620888818486853646", 10).unwrap().to_bytes_be();
        let p_y = BigUint::from_str_radix("7913123550914612057135582061699117755797758113868200992327595317370485234417808273674357776714522052694559358668442301647906991623400754234679697332299689255516547752391831738454121261248793568285885897998257357202903170202349380518443", 10).unwrap().to_bytes_be();

        let q_x_0 = BigUint::from_str_radix("13426761183630949215425595811885033211332897733228446437546263564078445562454176776915160094418980045665397361295624472103734543457352048745726512354895954850428989867542989474136256025045975283415690491751906307188562464175510373683338", 10).unwrap().to_bytes_be();
        let q_x_1 = BigUint::from_str_radix("20471601555918880743198170952645906008198510944268658573129351735028343217532386920456705632337352161031960990613816401042894531220068552819818037605513359562118363589199569321421558696125646867661360498323171027455638052943806292028610", 10).unwrap().to_bytes_be();
        let q_x_2 = BigUint::from_str_radix("3905053196875761830053608605277158152930144841844497593936739534395003062685449846381431331169369910535935138116320442345524758217411779027270883193856999691582831339845600938304719916501940381093815781408183227875600753651697934495980", 10).unwrap().to_bytes_be();
        
        let q_y_0 = BigUint::from_str_radix("8567517639523571619872938228644013584947463594196306323477160496987712111576624702939472765993995586889532559039169098780892505598589581147768095093536988446010255611523736706017580686335404469207486594272103717837888228343074699140243", 10).unwrap().to_bytes_be();
        let q_y_1 = BigUint::from_str_radix("3890537069205870914984502594450293167889863914413852788876350245583932846980126025043974070704295857226211547108005650399870458089721518559480870503159804530091559886149680718531004778697982910253701559194337987238111062202037698927752", 10).unwrap().to_bytes_be();
        let q_y_2 = BigUint::from_str_radix("10936269922612615564271188303104593362724754284143779051599749016735041389483971486958818324356025479751246744831831158558101688599198721653921723013062333636402617118847009085485166284126970598561393411916461254016145116183331671450721", 10).unwrap().to_bytes_be();

        let p_x = Fp::from_be_bytes(&base_field, &p_x, true).unwrap();
        let p_y = Fp::from_be_bytes(&base_field, &p_y, true).unwrap();

        let q_x_0 = Fp::from_be_bytes(&base_field, &q_x_0, true).unwrap();
        let q_x_1 = Fp::from_be_bytes(&base_field, &q_x_1, true).unwrap();
        let q_x_2 = Fp::from_be_bytes(&base_field, &q_x_2, true).unwrap();

        let q_y_0 = Fp::from_be_bytes(&base_field, &q_y_0, true).unwrap();
        let q_y_1 = Fp::from_be_bytes(&base_field, &q_y_1, true).unwrap();
        let q_y_2 = Fp::from_be_bytes(&base_field, &q_y_2, true).unwrap();

        let mut q_x = Fp3::zero(&extension_3);
        q_x.c0 = q_x_0;
        q_x.c1 = q_x_1;
        q_x.c2 = q_x_2;

        let mut q_y = Fp3::zero(&extension_3);
        q_y.c0 = q_y_0;
        q_y.c1 = q_y_1;
        q_y.c2 = q_y_2;

        let p = CurvePoint::point_from_xy(&curve, p_x, p_y);
        let q = TwistPoint::point_from_xy(&curve_twist, q_x, q_y);

        // let x = BigUint::from_str_radix("506464946133393486072777102926336625944849939610982267859828541006717966526573193706126370441346337661774335955699621", 10).unwrap();
        // println!("X len = {}", biguint_to_u64_vec(x.clone()).len());
        // println!("{:x}", biguint_to_u64_vec(x.clone())[0]);
        let w0 = BigUint::from_str_radix("7000705447348627246181409558336018323010329260726930841638672011287206690002601216854775649561085256265269640040570922609783227469279331691880282815325569032149343779036142830666859805506518426649197067288711084398033", 10).unwrap();
        let w1 = BigUint::from_str_radix("86482221941698704497288378992285180119495364068003923046442785886272123124361700722982503222189455144364945735564951562986", 10).unwrap();
        
        let x: Vec<u64> = vec![
            0x55c5b9b57b942ae8,
            0x3d52287d3dfd424a,
            0xcf1ff9d6a543deb7,
            0x820c9c5711ceeebc,
            0x549a2d44305d20fe,
            0x50f5c131afd70235,
            0xab3596c8617c5792,
            0x830c728d80f9d78b,
            0x6a7223ee72023d07,
            0xbc5d176b746af026,
            0xe959283d8f526663,
            0xc4d2263babf8941f,
            0x3848,
        ];

        assert!(p.check_on_curve());
        assert!(q.check_on_curve());

        let engine = super::CPInstance6 {
            x: x,
            x_is_negative: false,
            exp_w0: biguint_to_u64_vec(w0),
            exp_w1: biguint_to_u64_vec(w1),
            exp_w0_is_negative: true,
            base_field: &base_field,
            curve: &curve,
            curve_twist: &curve_twist,
            twist: twist,
            fp3_extension: &extension_3,
            fp6_extension: &extension_6,
        };

        b.iter(|| {
            engine.pair(&[p.clone()], &[q.clone()]).unwrap();
        });
    }
}