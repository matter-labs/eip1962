#[macro_use]
mod convenience;

pub mod generic;
pub mod bls12_381;
pub mod bls12_377;

#[cfg(feature = "eip_196")]
pub mod bn254;

pub(crate) fn print_single(r: &[u64]) {
    let mut limb_string = vec![];
    for limb in r {
        let t = format!("0x{:016x}", limb);
        limb_string.push(t);
    }

    println!("{}", limb_string.join(","));
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::representation::*;
    use crate::field::*;
    use crate::traits::*;
    use crate::fp::*;
    use crate::extension_towers::fp2::*;
    use num_bigint::BigUint;
    use num_traits::Num;
    use crate::pairings::TwistType;
    use crate::integers::MaxFieldUint;

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

    #[test]
    fn test_print_bls12_377_params() {
        let field = crate::field::new_field::<U384Repr>("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177", 10).unwrap();
        pretty_print_field_constants(&field);
    }

    fn print_bls12_engine_constants<FE: ElementRepr>(
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
        // non-residue is -1
        let mut fp_non_residue = modulus.clone();
        fp_non_residue -= BigUint::from(1u64);

        // fp2 non-residue is (1, 1)

        let fp2_non_residue_c0 = BigUint::from(1u64);
        let fp2_non_residue_c1 = BigUint::from(1u64);

        let b = BigUint::from(4u64);

        let beta = BigUint::from_str_radix("5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a00022e01fffffffefffe", 16).unwrap();
        let modulus_uint = MaxFieldUint::from_big_endian(&modulus.to_bytes_be());
        let field = field_from_modulus::<U384Repr>(&modulus_uint).unwrap();
        let beta = Fp::from_be_bytes(&field, &beta.to_bytes_be(), true).unwrap();
        print_single(beta.repr.as_ref());

        let p_x = BigUint::from_str_radix("3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507", 10).unwrap();
        let p_y = BigUint::from_str_radix("1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569", 10).unwrap();

        let q_x_0 = BigUint::from_str_radix("352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160", 10).unwrap();
        let q_x_1 = BigUint::from_str_radix("3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758", 10).unwrap();
        let q_y_0 = BigUint::from_str_radix("1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905", 10).unwrap();
        let q_y_1 = BigUint::from_str_radix("927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582", 10).unwrap();

        let twist_type = TwistType::M;

        print_bls12_engine_constants::<U384Repr>(
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
    fn print_bls12_377_constants() {
        let modulus = BigUint::from_str_radix("258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177", 10).unwrap();
        // non-residue is -5
        let mut fp_non_residue = modulus.clone();
        fp_non_residue -= BigUint::from(5u64);

        // fp2 non-residue is (0, 1)

        let fp2_non_residue_c0 = BigUint::from(0u64);
        let fp2_non_residue_c1 = BigUint::from(1u64);

        let b = BigUint::from(1u64);

        let p_x = BigUint::from_str_radix("008848defe740a67c8fc6225bf87ff5485951e2caa9d41bb188282c8bd37cb5cd5481512ffcd394eeab9b16eb21be9ef", 16).unwrap();
        let p_y = BigUint::from_str_radix("01914a69c5102eff1f674f5d30afeec4bd7fb348ca3e52d96d182ad44fb82305c2fe3d3634a9591afd82de55559c8ea6", 16).unwrap();

        let q_x_0 = BigUint::from_str_radix("018480be71c785fec89630a2a3841d01c565f071203e50317ea501f557db6b9b71889f52bb53540274e3e48f7c005196", 16).unwrap();
        let q_x_1 = BigUint::from_str_radix("00ea6040e700403170dc5a51b1b140d5532777ee6651cecbe7223ece0799c9de5cf89984bff76fe6b26bfefa6ea16afe", 16).unwrap();
        let q_y_0 = BigUint::from_str_radix("00690d665d446f7bd960736bcbb2efb4de03ed7274b49a58e458c282f832d204f2cf88886d8c7c2ef094094409fd4ddf", 16).unwrap();
        let q_y_1 = BigUint::from_str_radix("00f8169fd28355189e549da3151a70aa61ef11ac3d591bf12463b01acee304c24279b83f5e52270bd9a1cdd185eb8f93", 16).unwrap();

        let twist_type = TwistType::D;

        print_bls12_engine_constants::<U384Repr>(
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

    #[test]
    fn calculate_bls12_377_constants() {
        let mut ext_2 = super::bls12_377::BLS12_377_EXTENSION_2_FIELD.clone();
        let modulus = super::bls12_377::BLS12_377_MODULUS_UINT;
        ext_2.calculate_frobenius_coeffs(&modulus).expect("must calcualte frobenius for Fp2");

        println!("Frobenius coeffs for Fp2 c1");
        for (idx, c) in ext_2.frobenius_coeffs_c1.iter().enumerate() {
            println!("Frobenius coeff {}", idx);
            print_single(c.repr.as_ref());
        }

        let mut ext_6 = super::bls12_377::BLS12_377_EXTENSION_6_FIELD.clone();
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
        
        let mut ext_12 = super::bls12_377::BLS12_377_EXTENSION_12_FIELD.clone();
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

    #[test]
    fn calculate_bls12_381_g1_isogeny_constants() {
        let fp_field = &super::bls12_381::BLS12_381_FIELD;
        let a_biguint = BigUint::from_str_radix("144698a3b8e9433d693a02c96d4982b0ea985383ee66a8d8e8981aefd881ac98936f8da0e0f97f5cf428082d584c1d", 16).unwrap();
        let b_biguint = BigUint::from_str_radix("12e2908d11688030018b12e8753eee3b2016c1f0f24f4070a0b9c14fcef35ef55a23215a316ceaa5d1cc48e98e172be0", 16).unwrap();
        let a = Fp::from_be_bytes(fp_field, &a_biguint.to_bytes_be(), true).unwrap();
        let b = Fp::from_be_bytes(fp_field, &b_biguint.to_bytes_be(), true).unwrap();

        println!("G1 isogeny A =");
        print_single(a.repr.as_ref());

        println!("G1 isogeny B =");
        print_single(b.repr.as_ref());
    }

    #[test]
    fn calculate_bls12_381_fp_to_g1_mapping_constants() {
        let fp_field = &super::bls12_381::BLS12_381_FIELD;
        let z_biguint = BigUint::from_str_radix("11", 10).unwrap();
        let z = Fp::from_be_bytes(fp_field, &z_biguint.to_bytes_be(), true).unwrap();

        println!("Z = ");
        print_single(z.repr.as_ref());
    }

    #[test]
    fn calculate_bls12_381_fp2_to_g2_mapping_constants() {
        let fp_field = &super::bls12_381::BLS12_381_FIELD;
        let fp2_field = &super::bls12_381::BLS12_381_EXTENSION_2_FIELD; 

        let z_c0_biguint = BigUint::from_str_radix("2", 10).unwrap();
        let z_c1_biguint = BigUint::from_str_radix("1", 10).unwrap();

        let z_c0 = Fp::from_be_bytes(fp_field, &z_c0_biguint.to_bytes_be(), true).unwrap();
        let z_c1 = Fp::from_be_bytes(fp_field, &z_c1_biguint.to_bytes_be(), true).unwrap();

        let mut z = Fp2::zero(fp2_field);
        z.c0 = z_c0;
        z.c1 = z_c1;

        z.negate();

        println!("Z for G2 mapping:");
        println!("C0:");
        print_single(z.c0.repr.as_ref());
        println!("C1:");
        print_single(z.c1.repr.as_ref());
    }

    #[test]
    fn calculate_bls12_381_g2_isogeny_constants() {
        let fp_field = &super::bls12_381::BLS12_381_FIELD;

        let a_c0_biguint = BigUint::from_str_radix("0", 10).unwrap();
        let a_c1_biguint = BigUint::from_str_radix("240", 10).unwrap();

        let b_c0_biguint = BigUint::from_str_radix("1012", 10).unwrap();
        let b_c1_biguint = BigUint::from_str_radix("1012", 10).unwrap();

        let a_c0 = Fp::from_be_bytes(fp_field, &a_c0_biguint.to_bytes_be(), true).unwrap();
        let a_c1 = Fp::from_be_bytes(fp_field, &a_c1_biguint.to_bytes_be(), true).unwrap();

        let b_c0 = Fp::from_be_bytes(fp_field, &b_c0_biguint.to_bytes_be(), true).unwrap();
        let b_c1 = Fp::from_be_bytes(fp_field, &b_c1_biguint.to_bytes_be(), true).unwrap();

        println!("G2 isogeny A = ");
        println!("C0:");
        print_single(a_c0.repr.as_ref());
        println!("C1:");
        print_single(a_c1.repr.as_ref());

        println!("G2 isogeny B =");
        println!("C0:");
        print_single(b_c0.repr.as_ref());
        println!("C1:");
        print_single(b_c1.repr.as_ref());
    }

    #[cfg(feature = "mappings")]
    #[test]
    fn test_bls12_engine_map_fp_to_g1() {
        use crate::weierstrass::Group;

        let mut x = super::bls12_381::BLS12_381_FP_ONE.clone();
        x.double();
        x.double();
        x.square();

        let mapped = crate::engines::bls12_381::mapping::fp_to_g1(&x).unwrap();

        let should_be_zero = mapped.mul(&super::bls12_381::BLS12_381_SUBGROUP_ORDER[..]);

        assert!(should_be_zero.is_zero());
    }

    #[cfg(feature = "mappings")]
    #[test]
    fn test_bls12_engine_map_fp2_to_g2() {
        use crate::weierstrass::Group;
        
        let mut x = super::bls12_381::BLS12_381_FP2_ONE.clone();
        x.double();
        x.double();
        x.square();

        let mapped = crate::engines::bls12_381::mapping::fp2_to_g2(&x).unwrap();

        let should_be_zero = mapped.mul(&super::bls12_381::BLS12_381_SUBGROUP_ORDER[..]);

        assert!(should_be_zero.is_zero());
    }

    // These test vectors are from https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-07#page-143

    #[cfg(feature = "mappings")]
    #[test]
    fn test_bls12_ietf_vector_fp_to_g1_0() {
        use crate::weierstrass::Group;

        let biguint = BigUint::from_str_radix("07fdf49ea58e96015d61f6b5c9d1c8f277146a533ae7fbca2a8ef4c41055cd961fbc6e26979b5554e4b4f22330c0e16d", 16).unwrap();
        let x = Fp::from_be_bytes(&crate::engines::bls12_381::BLS12_381_FIELD, &biguint.to_bytes_be(), true).unwrap();

        let mapped = crate::engines::bls12_381::mapping::fp_to_g1(&x).unwrap();

        let (x, y) = mapped.clone().into_xy();
        // println!("Mapped to x = {}, y = {}", x, y);
        assert_eq!(format!("{}", x), "0x1223effdbb2d38152495a864d78eee14cb0992d89a241707abb03819a91a6d2fd65854ab9a69e9aacb0cbebfd490732c");
        assert_eq!(format!("{}", y), "0x0f925d61e0b235ecd945cbf0309291878df0d06e5d80d6b84aa4ff3e00633b26f9a7cb3523ef737d90e6d71e8b98b2d5");

        let should_be_zero = mapped.mul(&super::bls12_381::BLS12_381_SUBGROUP_ORDER[..]);

        assert!(should_be_zero.is_zero());
    }

    #[cfg(feature = "mappings")]
    #[test]
    fn test_bls12_ietf_vector_fp_to_g1_1() {
        use crate::weierstrass::Group;

        let biguint = BigUint::from_str_radix("1275ab3adbf824a169ed4b1fd669b49cf406d822f7fe90d6b2f8c601b5348436f89761bb1ad89a6fb1137cd91810e5d2", 16).unwrap();
        let x = Fp::from_be_bytes(&crate::engines::bls12_381::BLS12_381_FIELD, &biguint.to_bytes_be(), true).unwrap();

        let mapped = crate::engines::bls12_381::mapping::fp_to_g1(&x).unwrap();

        let (x, y) = mapped.clone().into_xy();
        // println!("Mapped to x = {}, y = {}", x, y);
        assert_eq!(format!("{}", x), "0x179d3fd0b4fb1da43aad06cea1fb3f828806ddb1b1fa9424b1e3944dfdbab6e763c42636404017da03099af0dcca0fd6");
        assert_eq!(format!("{}", y), "0x0d037cb1c6d495c0f5f22b061d23f1be3d7fe64d3c6820cfcd99b6b36fa69f7b4c1f4addba2ae7aa46fb25901ab483e4");


        let should_be_zero = mapped.mul(&super::bls12_381::BLS12_381_SUBGROUP_ORDER[..]);

        assert!(should_be_zero.is_zero());
    }

    #[cfg(feature = "mappings")]
    #[test]
    fn test_bls12_ietf_vector_fp_to_g1_2() {
        use crate::weierstrass::Group;

        let biguint = BigUint::from_str_radix("0e93d11d30de6d84b8578827856f5c05feef36083eef0b7b263e35ecb9b56e86299614a042e57d467fa20948e8564909", 16).unwrap();
        let x = Fp::from_be_bytes(&crate::engines::bls12_381::BLS12_381_FIELD, &biguint.to_bytes_be(), true).unwrap();

        let mapped = crate::engines::bls12_381::mapping::fp_to_g1(&x).unwrap();

        let (x, y) = mapped.clone().into_xy();
        // println!("Mapped to x = {}, y = {}", x, y);
        assert_eq!(format!("{}", x), "0x15aa66c77eded1209db694e8b1ba49daf8b686733afaa7b68c683d0b01788dfb0617a2e2d04c0856db4981921d3004af");
        assert_eq!(format!("{}", y), "0x0952bb2f61739dd1d201dd0a79d74cda3285403d47655ee886afe860593a8a4e51c5b77a22d2133e3a4280eaaaa8b788");


        let should_be_zero = mapped.mul(&super::bls12_381::BLS12_381_SUBGROUP_ORDER[..]);

        assert!(should_be_zero.is_zero());
    }

    #[cfg(feature = "mappings")]
    #[test]
    fn test_bls12_ietf_vector_fp_to_g1_3() {
        use crate::weierstrass::Group;

        let biguint = BigUint::from_str_radix("015a41481155d17074d20be6d8ec4d46632a51521cd9c916e265bd9b47343b3689979b50708c8546cbc2916b86cb1a3a", 16).unwrap();
        let x = Fp::from_be_bytes(&crate::engines::bls12_381::BLS12_381_FIELD, &biguint.to_bytes_be(), true).unwrap();

        let mapped = crate::engines::bls12_381::mapping::fp_to_g1(&x).unwrap();

        let (x, y) = mapped.clone().into_xy();
        // println!("Mapped to x = {}, y = {}", x, y);
        assert_eq!(format!("{}", x), "0x06328ce5106e837935e8da84bd9af473422e62492930aa5f460369baad9545defa468d9399854c23a75495d2a80487ee");
        assert_eq!(format!("{}", y), "0x094bfdfe3e552447433b5a00967498a3f1314b86ce7a7164c8a8f4131f99333b30a574607e301d5f774172c627fd0bca");


        let should_be_zero = mapped.mul(&super::bls12_381::BLS12_381_SUBGROUP_ORDER[..]);

        assert!(should_be_zero.is_zero());
    }

    // these test vectors are from https://tools.ietf.org/html/draft-irtf-cfrg-hash-to-curve-07#page-151

    #[cfg(feature = "mappings")]
    #[test]
    fn test_bls12_ietf_vector_fp2_to_g2_0() {
        use crate::weierstrass::Group;

        let fp_field = &crate::engines::bls12_381::BLS12_381_FIELD;

        let c0 = BigUint::from_str_radix("0e775d7827adf385b83e20e4445bd3fab21d7b4498426daf3c1d608b9d41e9edb5eda0df022e753b8bb4bc3bb7db4914", 16).unwrap();
        let c1 = BigUint::from_str_radix("025fbc07711ba267b7e70c82caa70a16fbb1d470ae24ceef307f5e2000751677820b7013ad4e25492dcf30052d3e5eca", 16).unwrap();

        let c0 = Fp::from_be_bytes(fp_field, &c0.to_bytes_be(), true).unwrap();
        let c1 = Fp::from_be_bytes(fp_field, &c1.to_bytes_be(), true).unwrap();

        let mut fe = crate::engines::bls12_381::BLS12_381_FP2_ZERO.clone();
        fe.c0 = c0;
        fe.c1 = c1;

        println!("Mapping {}", &fe);

        let mapped = crate::engines::bls12_381::mapping::fp2_to_g2(&fe).unwrap();

        let (x, y) = mapped.clone().into_xy();
        // println!("Mapped to x = {}, y = {}", x, y);
        assert_eq!(format!("{}", x.c0), "0x027e4bfada0b47f9f07e04aec463c7371e68f2fd0c738cd517932ea3801a35acf09db018deda57387b0f270f7a219e4d");
        assert_eq!(format!("{}", x.c1), "0x0d4333b77becbf9f9dfa3ca928002233d1ecc854b1447e5a71f751c9042d000f42db91c1d6649a5e0ad22bd7bf7398b8");
        assert_eq!(format!("{}", y.c0), "0x053674cba9ef516ddc218fedb37324e6c47de27f88ab7ef123b006127d738293c0277187f7e2f80a299a24d84ed03da7");
        assert_eq!(format!("{}", y.c1), "0x0cc76dc777ea0d447e02a41004f37a0a7b1fafb6746884e8d9fc276716ccf47e4e0899548a2ec71c2bdf1a2a50e876db");

        let should_be_zero = mapped.mul(&super::bls12_381::BLS12_381_SUBGROUP_ORDER[..]);

        assert!(should_be_zero.is_zero());
    }

    #[cfg(feature = "mappings")]
    #[test]
    fn test_bls12_ieft_vector_fp2_to_g2_1() {
        use crate::weierstrass::Group;

        let fp_field = &crate::engines::bls12_381::BLS12_381_FIELD;

        let c0 = BigUint::from_str_radix("045ab31ce4b5a8ba7c4b2851b64f063a66cd1223d3c85005b78e1beee65e33c90ceef0244e45fc45a5e1d6eab6644fdb", 16).unwrap();
        let c1 = BigUint::from_str_radix("1870a7dbfd2a1deb74015a3546b20f598041bf5d5202997956a94a368d30d3f70f18cdaa1d33ce970a4e16af961cbdcb", 16).unwrap();

        let c0 = Fp::from_be_bytes(fp_field, &c0.to_bytes_be(), true).unwrap();
        let c1 = Fp::from_be_bytes(fp_field, &c1.to_bytes_be(), true).unwrap();

        let mut fe = crate::engines::bls12_381::BLS12_381_FP2_ZERO.clone();
        fe.c0 = c0;
        fe.c1 = c1;

        println!("Mapping {}", &fe);

        let mapped = crate::engines::bls12_381::mapping::fp2_to_g2(&fe).unwrap();

        let (x, y) = mapped.clone().into_xy();
        // println!("Mapped to x = {}, y = {}", x, y);
        assert_eq!(format!("{}", x.c0), "0x09349f1cb5b2e55489dcd45a38545343451cc30a1681c57acd4fb0a6db125f8352c09f4a67eb7d1d8242cb7d3405f97b");
        assert_eq!(format!("{}", x.c1), "0x18f0f87b40af67c056915dbaf48534c592524e82c1c2b50c3734d02c0172c80df780a60b5683759298a3303c5d942778");
        assert_eq!(format!("{}", y.c0), "0x02f2d9deb2c7742512f5b8230bf0fd83ea42279d7d39779543c1a43b61c885982b611f6a7a24b514995e8a098496b811");
        assert_eq!(format!("{}", y.c1), "0x10a2ba341bc689ab947b7941ce6ef39be17acaab067bd32bd652b471ab0792c53a2bd03bdac47f96aaafe96e441f63c0");

        let should_be_zero = mapped.mul(&super::bls12_381::BLS12_381_SUBGROUP_ORDER[..]);

        assert!(should_be_zero.is_zero());
    }

    #[cfg(feature = "mappings")]
    #[test]
    fn test_bls12_ieft_vector_fp2_to_g2_2() {
        use crate::weierstrass::Group;

        let fp_field = &crate::engines::bls12_381::BLS12_381_FIELD;

        let c0 = BigUint::from_str_radix("0b6e6135a4cd31ba980ddbd115ac48abef7ec60e226f264d7befe002c165f3a496f36f76dd524efd75d17422558d10b4", 16).unwrap();
        let c1 = BigUint::from_str_radix("088fe329b054db8a6474f21a7fbfdf17b4c18044db299d9007af582c3d5f17d00e56d99921d4b5640fce44b05219b5de", 16).unwrap();

        let c0 = Fp::from_be_bytes(fp_field, &c0.to_bytes_be(), true).unwrap();
        let c1 = Fp::from_be_bytes(fp_field, &c1.to_bytes_be(), true).unwrap();

        let mut fe = crate::engines::bls12_381::BLS12_381_FP2_ZERO.clone();
        fe.c0 = c0;
        fe.c1 = c1;

        println!("Mapping {}", &fe);

        let mapped = crate::engines::bls12_381::mapping::fp2_to_g2(&fe).unwrap();

        let (x, y) = mapped.clone().into_xy();
        // println!("Mapped to x = {}, y = {}", x, y);
        assert_eq!(format!("{}", x.c0), "0x149fe43777d34f0d25430dea463889bd9393bdfb4932946db23671727081c629ebb98a89604f3433fba1c67d356a4af7");
        assert_eq!(format!("{}", x.c1), "0x19808ec5930a53c7cf5912ccce1cc33f1b3dcff24a53ce1cc4cba41fd6996dbed4843ccdd2eaf6a0cd801e562718d163");
        assert_eq!(format!("{}", y.c0), "0x04c0d6793a766233b2982087b5f4a254f261003ccb3262ea7c50903eecef3e871d1502c293f9e063d7d293f6384f4551");
        assert_eq!(format!("{}", y.c1), "0x04783e391c30c83f805ca271e353582fdf19d159f6a4c39b73acbb637a9b8ac820cfbe2738d683368a7c07ad020e3e33");

        let should_be_zero = mapped.mul(&super::bls12_381::BLS12_381_SUBGROUP_ORDER[..]);

        assert!(should_be_zero.is_zero());
    }

    #[cfg(feature = "mappings")]
    #[test]
    fn test_bls12_ieft_vector_fp2_to_g2_3() {
        use crate::weierstrass::Group;

        let fp_field = &crate::engines::bls12_381::BLS12_381_FIELD;

        let c0 = BigUint::from_str_radix("0f45b50647d67485295aa9eb2d91a877b44813677c67c8d35b2173ff3ba95f7bd0806f9ca8a1436b8b9d14ee81da4d7e", 16).unwrap();
        let c1 = BigUint::from_str_radix("03df16a66a05e4c1188c234788f43896e0565bfb64ac49b9639e6b284cc47dad73c47bb4ea7e677db8d496beb907fbb6", 16).unwrap();

        let c0 = Fp::from_be_bytes(fp_field, &c0.to_bytes_be(), true).unwrap();
        let c1 = Fp::from_be_bytes(fp_field, &c1.to_bytes_be(), true).unwrap();

        let mut fe = crate::engines::bls12_381::BLS12_381_FP2_ZERO.clone();
        fe.c0 = c0;
        fe.c1 = c1;

        println!("Mapping {}", &fe);

        let mapped = crate::engines::bls12_381::mapping::fp2_to_g2(&fe).unwrap();

        let (x, y) = mapped.clone().into_xy();
        // println!("Mapped to x = {}, y = {}", x, y);
        assert_eq!(format!("{}", x.c0), "0x0804152cbf8474669ad7d1796ab92d7ca21f32d8bed70898a748ed4e4e0ec557069003732fc86866d938538a2ae95552");
        assert_eq!(format!("{}", x.c1), "0x0b8e0094c886487870372eb6264613a6a087c7eb9804fab789be4e47a57b29eb19b1983a51165a1b5eb025865e9fc63a");
        assert_eq!(format!("{}", y.c0), "0x09e5c8242dd7281ad32c03fe4af3f19167770016255fb25ad9b67ec51d62fade31a1af101e8f6172ec2ee8857662be3a");
        assert_eq!(format!("{}", y.c1), "0x14c80f068ece15a3936bb00c3c883966f75b4e8d9ddde809c11f781ab92d23a2d1d103ad48f6f3bb158bf3e3a4063449");

        let should_be_zero = mapped.mul(&super::bls12_381::BLS12_381_SUBGROUP_ORDER[..]);

        assert!(should_be_zero.is_zero());
    }
}