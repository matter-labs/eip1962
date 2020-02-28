pub mod bls12_381;

#[cfg(test)]
mod test {
    use crate::representation::*;
    use crate::field::*;
    use crate::traits::*;
    use crate::fp::*;
    use crate::extension_towers::fp2::*;
    use num_bigint::BigUint;
    use num_traits::Num;
    use crate::pairings::TwistType;
    use crate::integers::MaxFieldUint;

    fn print_single(r: &[u64]) {
        let mut limb_string = vec![];
        for limb in r {
            let t = format!("0x{:016x}", limb);
            limb_string.push(t);
        }

        println!("{}", limb_string.join(","));
    }

    fn pretty_print_field_constants<E: ElementRepr, F: SizedPrimeField<Repr = E>>(field: &F) {
        println!("Modulus:");
        print_single(field.modulus().as_ref());
        println!("R:");
        print_single(field.mont_r().as_ref());
        println!("R2:");
        print_single(field.mont_r2().as_ref());
        println!("Inv:");
        println!("0x{:016x}", field.mont_inv())
    }

    #[test]
    fn test_print_bls12_381_params() {
        let field = crate::field::new_field::<U384Repr>("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        pretty_print_field_constants(&field);
    }

    fn print_bls12_constants<FE: ElementRepr>(
        modulus: BigUint, 
        twist_type: TwistType,
        b: BigUint,
        fp_non_residue: BigUint,
        fp2_non_residue_c0: BigUint,
        fp2_non_residue_c1: BigUint,
        generator_g1_x: BigUint,
        generator_g1_y: BigUint,
        generator_g2_x_0: BigUint,
        generator_g2_x_1: BigUint,
        generator_g2_y_0: BigUint,
        generator_g2_y_1: BigUint,
    ) {
        let modulus_uint = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());
        let field = field_from_modulus::<FE>(&modulus_uint).unwrap();
        let b_fp = Fp::from_be_bytes(&field, &b.to_bytes_be(), true).unwrap();
        let fp_non_residue = Fp::from_be_bytes(&field, &fp_non_residue.to_bytes_be(), true).unwrap();
        let extension_2 = Extension2::new(fp_non_residue.clone());
        let mut fp2_non_residue = Fp2::zero(&extension_2);
        let fp2_non_residue_c0 = Fp::from_be_bytes(&field, &fp2_non_residue_c0.to_bytes_be(), true).unwrap();
        let fp2_non_residue_c1 = Fp::from_be_bytes(&field, &fp2_non_residue_c1.to_bytes_be(), true).unwrap();
        fp2_non_residue.c0 = fp2_non_residue_c0;
        fp2_non_residue.c1 = fp2_non_residue_c1;
    
        let fp2_non_residue_inv = fp2_non_residue.inverse().unwrap();
        let b_fp2 = match twist_type {
            TwistType::D => {
                let mut b_fp2 = fp2_non_residue_inv.clone();
                b_fp2.mul_by_fp(&b_fp);
    
                b_fp2
            },
            TwistType::M => {
                let mut b_fp2 = fp2_non_residue.clone();
                b_fp2.mul_by_fp(&b_fp);
    
                b_fp2
            },
        };
    
        let g1_generator_x = Fp::from_be_bytes(&field, &generator_g1_x.to_bytes_be(), true).unwrap();
        let g1_generator_y = Fp::from_be_bytes(&field, &generator_g1_y.to_bytes_be(), true).unwrap();
    
        let g1_generator_x_c0 = Fp::from_be_bytes(&field, &generator_g2_x_0.to_bytes_be(), true).unwrap();
        let g1_generator_x_c1 = Fp::from_be_bytes(&field, &generator_g2_x_1.to_bytes_be(), true).unwrap();
        let g1_generator_y_c0 = Fp::from_be_bytes(&field, &generator_g2_y_0.to_bytes_be(), true).unwrap();
        let g1_generator_y_c1 = Fp::from_be_bytes(&field, &generator_g2_y_1.to_bytes_be(), true).unwrap();

        println!("B:");
        print_single(&b_fp.repr.as_ref());

        println!("B for twist:");
        println!("C0:");
        print_single(&b_fp2.c0.repr.as_ref());
        println!("C1:");
        print_single(&b_fp2.c1.repr.as_ref());

        println!("Fp non-residue:");
        print_single(&fp_non_residue.repr.as_ref());

        println!("Fp2 non-residue:");
        println!("C0:");
        print_single(&fp2_non_residue.c0.repr.as_ref());
        println!("C1:");
        print_single(&fp2_non_residue.c1.repr.as_ref());

        println!("G1 generator");
        println!("X:");
        print_single(&g1_generator_x.repr.as_ref());
        println!("Y:");
        print_single(&g1_generator_y.repr.as_ref());

        println!("G2 generator");
        println!("X:");
        println!("C0:");
        print_single(&g1_generator_x_c0.repr.as_ref());
        println!("C1:");
        print_single(&g1_generator_x_c1.repr.as_ref());
        println!("Y:");
        println!("C0:");
        print_single(&g1_generator_y_c0.repr.as_ref());
        println!("C1:");
        print_single(&g1_generator_y_c1.repr.as_ref());
    }

    #[test]
    fn print_bls12_381_constants() {
        let modulus = BigUint::from_str_radix("4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787", 10).unwrap();
        let group_order = BigUint::from_str_radix("52435875175126190479447740508185965837690552500527637822603658699938581184513", 10).unwrap();

        // non-residue is -1
        let mut fp_non_residue = modulus.clone();
        fp_non_residue -= BigUint::from(1u64);

        // fp2 non-residue is (1, 1)

        let fp2_non_residue_c0 = BigUint::from(1u64);
        let fp2_non_residue_c1 = BigUint::from(1u64);

        let b = BigUint::from(4u64);

        let p_x = BigUint::from_str_radix("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507", 10).unwrap();
        let p_y = BigUint::from_str_radix("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569", 10).unwrap();

        let q_x_0 = BigUint::from_str_radix("352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160", 10).unwrap();
        let q_x_1 = BigUint::from_str_radix("3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758", 10).unwrap();
        let q_y_0 = BigUint::from_str_radix("1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905", 10).unwrap();
        let q_y_1 = BigUint::from_str_radix("927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582", 10).unwrap();

        let twist_type = TwistType::M;

        print_bls12_constants::<U384Repr>(
            modulus,
            twist_type,
            b,
            fp_non_residue,
            fp2_non_residue_c0,
            fp2_non_residue_c1,
            p_x,
            p_y,
            q_x_0,
            q_x_1,
            q_y_0,
            q_y_1,
        );
    }

    #[test]
    fn calculate_bls12_381_constants() {
        let mut ext_2 = super::bls12_381::BLS12_381_EXTENSION_2_FIELD.clone();
        let modulus = super::bls12_381::BLS12_381_MODULUS_UINT;
        ext_2.calculate_frobenius_coeffs(&modulus).expect("must calcualte frobenius for Fp2");

        println!("Frobenius coeffs for Fp2 c1");
        for (idx, c) in ext_2.frobenius_coeffs_c1.iter().enumerate() {
            println!("Frobenius coeff {}", idx);
            print_single(c.repr.as_ref());
        }

        let mut ext_6 = super::bls12_381::BLS12_381_EXTENSION_6_FIELD.clone();
        ext_6.calculate_frobenius_coeffs_optimized(&modulus).expect("must calcualte frobenius for Fp6");

        println!("Frobenius coeffs for Fp6 c1");
        for (idx, c) in ext_6.frobenius_coeffs_c1.iter().enumerate() {
            println!("Frobenius coeff {}", idx);
            println!("C0:");
            print_single(c.c0.repr.as_ref());
            println!("C1:");
            print_single(c.c1.repr.as_ref());
        }

        println!("Frobenius coeffs for Fp6 c2");
        for (idx, c) in ext_6.frobenius_coeffs_c2.iter().enumerate() {
            println!("Frobenius coeff {}", idx);
            println!("C0:");
            print_single(c.c0.repr.as_ref());
            println!("C1:");
            print_single(c.c1.repr.as_ref());
        }   
        
        let mut ext_12 = super::bls12_381::BLS12_381_EXTENSION_12_FIELD.clone();
        ext_12.calculate_frobenius_coeffs_optimized(&modulus).expect("must calcualte frobenius for Fp12");

        println!("Frobenius coeffs for Fp12 c1");
        for (idx, c) in ext_12.frobenius_coeffs_c1.iter().enumerate() {
            println!("Frobenius coeff {}", idx);
            println!("C0:");
            print_single(c.c0.repr.as_ref());
            println!("C1:");
            print_single(c.c1.repr.as_ref());
        }
    }
}