use crate::fp::*;
use crate::representation::*;
use crate::field::*;
use crate::traits::*;
use crate::extension_towers::fp2::Fp2;

pub mod simple_swu;
pub mod isogeny;
pub mod parameters;

#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug)]
pub enum Sign {
    Zero,
    SignPlus,
    SignMinus
}

fn sign_of_fp<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >(
    el: &Fp<'a, E, F>
) -> Sign {
    if el.is_zero() {
        return Sign::Zero;
    }

    let el_repr = el.into_repr();

    if el_repr.as_ref()[0] & 1 == 0 {
        // 0 mod 2
        return Sign::SignPlus;
    } else {
        return Sign::SignMinus;
    }
}

fn sign_of_fp2<'a, E: ElementRepr, F: SizedPrimeField<Repr = E> >(
    el: &Fp2<'a, E, F>
) -> Sign {
    // compare c_0 and then c_1
    let sign_from_c0 = sign_of_fp(&el.c0);

    match sign_from_c0 {
        s @ Sign::SignMinus | 
        s @ Sign::SignPlus => {
            return s;
        },
        Sign::Zero => {
            return sign_of_fp(&el.c1);
        }
    }
}

mod constants {
    use super::isogeny::*;
    use super::simple_swu::*;
    use crate::engines::bls12_381::*;
    use crate::weierstrass::*;
    use crate::field::*;
    use crate::traits::*;
    use num_bigint::BigUint;
    use num_traits::Num;
    use crate::fp::*;
    use crate::extension_towers::fp2::{Extension2, Fp2};

    pub(crate) fn str_radix_into_field<'a>(s: &str, radix: u32, field: &'a PrimeField<U384Repr>) -> Fp<'a, U384Repr, PrimeField<U384Repr>> {
        let biguint = BigUint::from_str_radix(s, radix).unwrap();
        let el = Fp::from_be_bytes(field, &biguint.to_bytes_be(), true).unwrap();

        el
    }

    pub(crate) fn str_radix_into_ext2<'a>(c0: &str, c1: &str, radix: u32, extension: &'a Extension2<'a, U384Repr, PrimeField<U384Repr>>) -> Fp2<'a, U384Repr, PrimeField<U384Repr>> {
        let biguint_c0 = BigUint::from_str_radix(c0, radix).unwrap();
        let biguint_c1 = BigUint::from_str_radix(c1, radix).unwrap();
        let c0 = Fp::from_be_bytes(extension.field(), &biguint_c0.to_bytes_be(), true).unwrap();
        let c1 = Fp::from_be_bytes(extension.field(), &biguint_c1.to_bytes_be(), true).unwrap();

        let mut el = Fp2::zero(extension);
        el.c0 = c0;
        el.c1 = c1;

        el
    }

    pub(crate) fn calculate_bls12_381_g1_mapping_params<'a>(field: &'a PrimeField<U384Repr>) -> (
        SwuParameters<CurveOverFpParameters<'a, U384Repr, PrimeField<U384Repr>>>,
        IsogenyParameters<CurveOverFpParameters<'a, U384Repr, PrimeField<U384Repr>>>
    ) {
        let mut minus_z_inv = BLS12_381_G1_SWU_Z.inverse().unwrap();
        minus_z_inv.negate();

        let mut minus_b_by_a = BLS12_381_G1_ISOGENY_B;
        minus_b_by_a.negate();

        let a_inv = BLS12_381_G1_ISOGENY_A.inverse().unwrap();
        minus_b_by_a.mul_assign(&a_inv);

        let swu = SwuParameters::<_> {
            z: BLS12_381_G1_SWU_Z.clone(),
            minus_b_by_a,
            minus_z_inv
        };

        let mut iso = IsogenyParameters::<_> {
            map_degree: 15,
            k1: vec![
                str_radix_into_field("11a05f2b1e833340b809101dd99815856b303e88a2d7005ff2627b56cdb4e2c85610c2d5f2e62d6eaeac1662734649b7", 16, field),
                str_radix_into_field("17294ed3e943ab2f0588bab22147a81c7c17e75b2f6a8417f565e33c70d1e86b4838f2a6f318c356e834eef1b3cb83bb", 16, field),
                str_radix_into_field("d54005db97678ec1d1048c5d10a9a1bce032473295983e56878e501ec68e25c958c3e3d2a09729fe0179f9dac9edcb0", 16, field),
                str_radix_into_field("1778e7166fcc6db74e0609d307e55412d7f5e4656a8dbf25f1b33289f1b330835336e25ce3107193c5b388641d9b6861", 16, field),
                str_radix_into_field("e99726a3199f4436642b4b3e4118e5499db995a1257fb3f086eeb65982fac18985a286f301e77c451154ce9ac8895d9", 16, field),
                str_radix_into_field("1630c3250d7313ff01d1201bf7a74ab5db3cb17dd952799b9ed3ab9097e68f90a0870d2dcae73d19cd13c1c66f652983", 16, field),
                str_radix_into_field("d6ed6553fe44d296a3726c38ae652bfb11586264f0f8ce19008e218f9c86b2a8da25128c1052ecaddd7f225a139ed84", 16, field),
                str_radix_into_field("17b81e7701abdbe2e8743884d1117e53356de5ab275b4db1a682c62ef0f2753339b7c8f8c8f475af9ccb5618e3f0c88e", 16, field),
                str_radix_into_field("80d3cf1f9a78fc47b90b33563be990dc43b756ce79f5574a2c596c928c5d1de4fa295f296b74e956d71986a8497e317", 16, field),
                str_radix_into_field("169b1f8e1bcfa7c42e0c37515d138f22dd2ecb803a0c5c99676314baf4bb1b7fa3190b2edc0327797f241067be390c9e", 16, field),
                str_radix_into_field("10321da079ce07e272d8ec09d2565b0dfa7dccdde6787f96d50af36003b14866f69b771f8c285decca67df3f1605fb7b", 16, field),
                str_radix_into_field("6e08c248e260e70bd1e962381edee3d31d79d7e22c837bc23c0bf1bc24c6b68c24b1b80b64d391fa9c8ba2e8ba2d229", 16, field),
            ],
            k2: vec![
                str_radix_into_field("8ca8d548cff19ae18b2e62f4bd3fa6f01d5ef4ba35b48ba9c9588617fc8ac62b558d681be343df8993cf9fa40d21b1c", 16, field),
                str_radix_into_field("12561a5deb559c4348b4711298e536367041e8ca0cf0800c0126c2588c48bf5713daa8846cb026e9e5c8276ec82b3bff", 16, field),
                str_radix_into_field("b2962fe57a3225e8137e629bff2991f6f89416f5a718cd1fca64e00b11aceacd6a3d0967c94fedcfcc239ba5cb83e19", 16, field),
                str_radix_into_field("3425581a58ae2fec83aafef7c40eb545b08243f16b1655154cca8abc28d6fd04976d5243eecf5c4130de8938dc62cd8", 16, field),
                str_radix_into_field("13a8e162022914a80a6f1d5f43e7a07dffdfc759a12062bb8d6b44e833b306da9bd29ba81f35781d539d395b3532a21e", 16, field),
                str_radix_into_field("e7355f8e4e667b955390f7f0506c6e9395735e9ce9cad4d0a43bcef24b8982f7400d24bc4228f11c02df9a29f6304a5", 16, field),
                str_radix_into_field("772caacf16936190f3e0c63e0596721570f5799af53a1894e2e073062aede9cea73b3538f0de06cec2574496ee84a3a", 16, field),
                str_radix_into_field("14a7ac2a9d64a8b230b3f5b074cf01996e7f63c21bca68a81996e1cdf9822c580fa5b9489d11e2d311f7d99bbdcc5a5e", 16, field),
                str_radix_into_field("a10ecf6ada54f825e920b3dafc7a3cce07f8d1d7161366b74100da67f39883503826692abba43704776ec3a79a1d641", 16, field),
                str_radix_into_field("95fc13ab9e92ad4476d6e3eb3a56680f682b4ee96f7d03776df533978f31c1593174e4b4b7865002d6384d168ecdd0a", 16, field),
                str_radix_into_field("1", 16, field),
            ],
            k3: vec![
                str_radix_into_field("90d97c81ba24ee0259d1f094980dcfa11ad138e48a869522b52af6c956543d3cd0c7aee9b3ba3c2be9845719707bb33", 16, field),
                str_radix_into_field("134996a104ee5811d51036d776fb46831223e96c254f383d0f906343eb67ad34d6c56711962fa8bfe097e75a2e41c696", 16, field),
                str_radix_into_field("cc786baa966e66f4a384c86a3b49942552e2d658a31ce2c344be4b91400da7d26d521628b00523b8dfe240c72de1f6", 16, field),
                str_radix_into_field("1f86376e8981c217898751ad8746757d42aa7b90eeb791c09e4a3ec03251cf9de405aba9ec61deca6355c77b0e5f4cb", 16, field),
                str_radix_into_field("8cc03fdefe0ff135caf4fe2a21529c4195536fbe3ce50b879833fd221351adc2ee7f8dc099040a841b6daecf2e8fedb", 16, field),
                str_radix_into_field("16603fca40634b6a2211e11db8f0a6a074a7d0d4afadb7bd76505c3d3ad5544e203f6326c95a807299b23ab13633a5f0", 16, field),
                str_radix_into_field("4ab0b9bcfac1bbcb2c977d027796b3ce75bb8ca2be184cb5231413c4d634f3747a87ac2460f415ec961f8855fe9d6f2", 16, field),
                str_radix_into_field("987c8d5333ab86fde9926bd2ca6c674170a05bfe3bdd81ffd038da6c26c842642f64550fedfe935a15e4ca31870fb29", 16, field),
                str_radix_into_field("9fc4018bd96684be88c9e221e4da1bb8f3abd16679dc26c1e8b6e6a1f20cabe69d65201c78607a360370e577bdba587", 16, field),
                str_radix_into_field("e1bba7a1186bdb5223abde7ada14a23c42a0ca7915af6fe06985e7ed1e4d43b9b3f7055dd4eba6f2bafaaebca731c30", 16, field),
                str_radix_into_field("19713e47937cd1be0dfd0b8f1d43fb93cd2fcbcb6caf493fd1183e416389e61031bf3a5cce3fbafce813711ad011c132", 16, field),
                str_radix_into_field("18b46a908f36f6deb918c143fed2edcc523559b8aaf0c2462e6bfe7f911f643249d9cdf41b44d606ce07c8a4d0074d8e", 16, field),
                str_radix_into_field("b182cac101b9399d155096004f53f447aa7b12a3426b08ec02710e807b4633f06c851c1919211f20d4c04f00b971ef8", 16, field),
                str_radix_into_field("245a394ad1eca9b72fc00ae7be315dc757b3b080d4c158013e6632d3c40659cc6cf90ad1c232a6442d9d3f5db980133", 16, field),
                str_radix_into_field("5c129645e44cf1102a159f748c4a3fc5e673d81d7e86568d9ab0f5d396a7ce46ba1049b6579afb7866b1e715475224b", 16, field),
                str_radix_into_field("15e6be4e990f03ce4ea50b3b42df2eb5cb181d8f84965a3957add4fa95af01b2b665027efec01c7704b456be69c8b604", 16, field),
            ],
            k4: vec![
                str_radix_into_field("16112c4c3a9c98b252181140fad0eae9601a6de578980be6eec3232b5be72e7a07f3688ef60c206d01479253b03663c1", 16, field),
                str_radix_into_field("1962d75c2381201e1a0cbd6c43c348b885c84ff731c4d59ca4a10356f453e01f78a4260763529e3532f6102c2e49a03d", 16, field),
                str_radix_into_field("58df3306640da276faaae7d6e8eb15778c4855551ae7f310c35a5dd279cd2eca6757cd636f96f891e2538b53dbf67f2", 16, field),
                str_radix_into_field("16b7d288798e5395f20d23bf89edb4d1d115c5dbddbcd30e123da489e726af41727364f2c28297ada8d26d98445f5416", 16, field),
                str_radix_into_field("be0e079545f43e4b00cc912f8228ddcc6d19c9f0f69bbb0542eda0fc9dec916a20b15dc0fd2ededda39142311a5001d", 16, field),
                str_radix_into_field("8d9e5297186db2d9fb266eaac783182b70152c65550d881c5ecd87b6f0f5a6449f38db9dfa9cce202c6477faaf9b7ac", 16, field),
                str_radix_into_field("166007c08a99db2fc3ba8734ace9824b5eecfdfa8d0cf8ef5dd365bc400a0051d5fa9c01a58b1fb93d1a1399126a775c", 16, field),
                str_radix_into_field("16a3ef08be3ea7ea03bcddfabba6ff6ee5a4375efa1f4fd7feb34fd206357132b920f5b00801dee460ee415a15812ed9", 16, field),
                str_radix_into_field("1866c8ed336c61231a1be54fd1d74cc4f9fb0ce4c6af5920abc5750c4bf39b4852cfe2f7bb9248836b233d9d55535d4a", 16, field),
                str_radix_into_field("167a55cda70a6e1cea820597d94a84903216f763e13d87bb5308592e7ea7d4fbc7385ea3d529b35e346ef48bb8913f55", 16, field),
                str_radix_into_field("4d2f259eea405bd48f010a01ad2911d9c6dd039bb61a6290e591b36e636a5c871a5c29f4f83060400f8b49cba8f6aa8", 16, field),
                str_radix_into_field("accbb67481d033ff5852c1e48c50c477f94ff8aefce42d28c0f9a88cea7913516f968986f7ebbea9684b529e2561092", 16, field),
                str_radix_into_field("ad6b9514c767fe3c3613144b45f1496543346d98adf02267d5ceef9a00d9b8693000763e3b90ac11e99b138573345cc", 16, field),
                str_radix_into_field("2660400eb2e4f3b628bdd0d53cd76f2bf565b94e72927c1cb748df27942480e420517bd8714cc80d1fadc1326ed06f7", 16, field),
                str_radix_into_field("e0fa1d816ddc03e6b24255e0d7819c171c40f65e273b853324efcd6356caa205ca2f570f13497804415473a1d634b8f", 16, field),
                str_radix_into_field("1", 16, field),
            ],
        };

        let zero = Fp::zero(field);

        iso.k1.resize(iso.map_degree + 1, zero.clone());
        iso.k2.resize(iso.map_degree + 1, zero.clone());
        iso.k3.resize(iso.map_degree + 1, zero.clone());
        iso.k4.resize(iso.map_degree + 1, zero.clone());

        (swu, iso)
    }

    pub(crate) fn calculate_bls12_381_g2_mapping_params<'a>(extension: &'a Extension2<'a, U384Repr, PrimeField<U384Repr>>) -> (
        SwuParameters<CurveOverFp2Parameters<'a, U384Repr, PrimeField<U384Repr>>>,
        IsogenyParameters<CurveOverFp2Parameters<'a, U384Repr, PrimeField<U384Repr>>>
    ) {
        let mut minus_z_inv = BLS12_381_G2_SWU_Z.inverse().unwrap();
        minus_z_inv.negate();

        let mut minus_b_by_a = BLS12_381_G2_ISOGENY_B;
        minus_b_by_a.negate();

        let a_inv = BLS12_381_G2_ISOGENY_A.inverse().unwrap();
        minus_b_by_a.mul_assign(&a_inv);

        let swu = SwuParameters::<_> {
            z: BLS12_381_G2_SWU_Z.clone(),
            minus_b_by_a,
            minus_z_inv
        };

        let mut iso = IsogenyParameters::<_> {
            map_degree: 3,
            k1: vec![
                str_radix_into_ext2(
                    "5c759507e8e333ebb5b7a9a47d7ed8532c52d39fd3a042a88b58423c50ae15d5c2638e343d9c71c6238aaaaaaaa97d6", 
                    "5c759507e8e333ebb5b7a9a47d7ed8532c52d39fd3a042a88b58423c50ae15d5c2638e343d9c71c6238aaaaaaaa97d6",
                    16,
                    extension
                ),
                str_radix_into_ext2(
                    "0", 
                    "11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71a",
                    16,
                    extension
                ),
                str_radix_into_ext2(
                    "11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71e", 
                    "8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffffffffe38d",
                    16,
                    extension
                ),
                str_radix_into_ext2(
                    "171d6541fa38ccfaed6dea691f5fb614cb14b4e7f4e810aa22d6108f142b85757098e38d0f671c7188e2aaaaaaaa5ed1", 
                    "0",
                    16,
                    extension
                ),
            ],
            k2: vec![
                str_radix_into_ext2(
                    "0", 
                    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaa63",
                    16,
                    extension
                ),
                str_radix_into_ext2(
                    "c", 
                    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaa9f",
                    16,
                    extension
                ),
                str_radix_into_ext2(
                    "1", 
                    "0",
                    16,
                    extension
                ),
            ],
            k3: vec![
                str_radix_into_ext2(
                    "1530477c7ab4113b59a4c18b076d11930f7da5d4a07f649bf54439d87d27e500fc8c25ebf8c92f6812cfc71c71c6d706", 
                    "1530477c7ab4113b59a4c18b076d11930f7da5d4a07f649bf54439d87d27e500fc8c25ebf8c92f6812cfc71c71c6d706",
                    16,
                    extension
                ),
                str_radix_into_ext2(
                    "0", 
                    "5c759507e8e333ebb5b7a9a47d7ed8532c52d39fd3a042a88b58423c50ae15d5c2638e343d9c71c6238aaaaaaaa97be",
                    16,
                    extension
                ),
                str_radix_into_ext2(
                    "11560bf17baa99bc32126fced787c88f984f87adf7ae0c7f9a208c6b4f20a4181472aaa9cb8d555526a9ffffffffc71c", 
                    "8ab05f8bdd54cde190937e76bc3e447cc27c3d6fbd7063fcd104635a790520c0a395554e5c6aaaa9354ffffffffe38f",
                    16,
                    extension
                ),
                str_radix_into_ext2(
                    "124c9ad43b6cf79bfbf7043de3811ad0761b0f37a1e26286b0e977c69aa274524e79097a56dc4bd9e1b371c71c718b10", 
                    "0",
                    16,
                    extension
                ),
            ],
            k4: vec![
                str_radix_into_ext2(
                    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffa8fb", 
                    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffa8fb",
                    16,
                    extension
                ),
                str_radix_into_ext2(
                    "0", 
                    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffa9d3",
                    16,
                    extension
                ),
                str_radix_into_ext2(
                    "12", 
                    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaa99",
                    16,
                    extension
                ),
                str_radix_into_ext2(
                    "1", 
                    "0",
                    16,
                    extension
                ),
            ],
        };

        let zero = Fp2::zero(extension);

        iso.k1.resize(iso.map_degree + 1, zero.clone());
        iso.k2.resize(iso.map_degree + 1, zero.clone());
        iso.k3.resize(iso.map_degree + 1, zero.clone());
        iso.k4.resize(iso.map_degree + 1, zero.clone());

        (swu, iso)
    }

    
}

#[cfg(test)]
mod test {
    use super::constants::*;
    use super::isogeny::*;
    use super::simple_swu::*;
    use crate::engines::bls12_381::*;
    use crate::weierstrass::curve::*;
    use super::sign_of_fp2;

    #[test]
    fn test_sign_of_zero() {
        let zero = str_radix_into_ext2("0", "0", 10, &BLS12_381_EXTENSION_2_FIELD);
        let sgn = sign_of_fp2(&zero);
        match sgn {
            super::Sign::Zero => {},
            _ => {panic!("invalid sign of zero")}
        }
    }

    #[test]
    fn test_mapping_g1_0() {
        let random_el = str_radix_into_field("42", 10, &BLS12_381_FIELD);
        let (swu, iso) = calculate_bls12_381_g1_mapping_params(&BLS12_381_FIELD);

        let (x_prime, y_prime) = simplified_swu_fp(
            &random_el, 
            &swu, 
            &BLS12_381_G1_CURVE_ISOGENY,
            None
        );

        let (x, y) = apply_isogeny_map(
            &x_prime,
            &y_prime,
            &iso,
            &BLS12_381_G1_CURVE_PARAMETERS
        );

        let point = CurvePoint::point_from_xy(&BLS12_381_G1_CURVE, x, y);
        assert!(point.is_on_curve());
    }

    #[test]
    fn test_mapping_g2_0() {
        let random_el = str_radix_into_ext2("42", "42", 10, &BLS12_381_EXTENSION_2_FIELD);
        let (swu, iso) = calculate_bls12_381_g2_mapping_params(&BLS12_381_EXTENSION_2_FIELD);

        let (x_prime, y_prime) = simplified_swu_fp2(
            &random_el, 
            &swu, 
            &BLS12_381_G2_CURVE_ISOGENY,
            None,
        );

        let (x, y) = apply_isogeny_map(
            &x_prime,
            &y_prime,
            &iso,
            &BLS12_381_G2_CURVE_PARAMETERS
        );

        let point = CurvePoint::point_from_xy(&BLS12_381_G2_CURVE, x, y);
        assert!(point.is_on_curve());
    }

}