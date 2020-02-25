use crate::pairings::TwistType;
use crate::representation::*;
use crate::field::*;
use num_bigint::BigUint;
use crate::integers::MaxFieldUint;
use crate::fp::*;
use crate::extension_towers::fp2::*;
use crate::traits::*;
use num_traits::Num;

fn generate_bls12_spec_params<FE: ElementRepr>(
    modulus: BigUint, 
    twist_type: TwistType,
    b: BigUint,
    main_subgroup_order: BigUint,
    fp_non_residue: BigUint,
    fp2_non_residue_c0: BigUint,
    fp2_non_residue_c1: BigUint,
    generator_g1_x: BigUint,
    generator_g1_y: BigUint,
    generator_g2_x_0: BigUint,
    generator_g2_x_1: BigUint,
    generator_g2_y_0: BigUint,
    generator_g2_y_1: BigUint,
    x: BigUint,
    x_is_negative: bool,
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
    println!("BLS12 Parameters");
    println!("Base field modulus = {}", field.modulus());
    println!("B coefficient = {}", b_fp);
    println!("Main subgroup order = 0x{}", main_subgroup_order.to_str_radix(16));

    println!("Extension tower:");
    println!("Fp2 construction:");
    println!("Fp quadratic non-residue = {}", fp_non_residue);
    println!("Fp6/Fp12 construction:");
    println!("Fp2 cubic non-residue c0 = {}", fp2_non_residue.c0);
    println!("Fp2 cubic non-residue c1 = {}", fp2_non_residue.c1);

    println!("Twist parameters:");
    match twist_type {
        TwistType::D => {
            println!("Twist type: D");
        },
        TwistType::M => {
            println!("Twist type: M");
        },
    };
    println!("B coefficient for twist c0 = {}", b_fp2.c0);
    println!("B coefficient for twist c1 = {}", b_fp2.c1);

    println!("Generators:");
    println!("G1:");
    println!("X = {}", g1_generator_x);
    println!("Y = {}", g1_generator_y);

    println!("G2:");
    println!("X c0 = {}", g1_generator_x_c0);
    println!("X c1 = {}", g1_generator_x_c1);
    println!("Y c0 = {}", g1_generator_y_c0);
    println!("Y c1 = {}", g1_generator_y_c1);

    println!("Pairing parameters:");
    println!("|x| (miller loop scalar) = 0x{}", x.to_str_radix(16));
    println!("x is negative = {}", x_is_negative);
}


#[test]
fn print_bls12_381_parameters() {
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

    let x = BigUint::from(0xd201000000010000u64);
    let x_is_negative = true;

    let twist_type = TwistType::M;

    generate_bls12_spec_params::<U384Repr>(
        modulus,
        twist_type,
        b,
        group_order,
        fp_non_residue,
        fp2_non_residue_c0,
        fp2_non_residue_c1,
        p_x,
        p_y,
        q_x_0,
        q_x_1,
        q_y_0,
        q_y_1,
        x,
        x_is_negative
    );
}