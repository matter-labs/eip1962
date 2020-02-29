use crate::field::*;
use crate::fp::*;
use crate::extension_towers::fp2::*;
use crate::extension_towers::fp6_as_3_over_2::*;
use crate::extension_towers::fp12_as_2_over3_over_2::*;
use crate::weierstrass::*;
use crate::weierstrass::curve::*;
use crate::pairings::bls12::*;
use crate::pairings::TwistType;
use crate::integers::MaxFieldUint;

const REPR_ZERO: U384Repr = U384Repr([0,0,0,0,0,0]);

pub const BLS12_377_MODULUS_UINT: MaxFieldUint = MaxFieldUint::from_limbs(
    [
        0x8508c00000000001,0x170b5d4430000000,0x1ef3622fba094800,0x1a22d9f300f5138f,
        0xc63b05c06ca1493b,0x01ae3a4617c510ea, 0x0, 0x0,
        0x0, 0x0, 0x0, 0x0,
        0x0, 0x0, 0x0, 0x0
    ]
);

pub const BLS12_377_MODULUS: U384Repr = U384Repr([0x8508c00000000001,0x170b5d4430000000,0x1ef3622fba094800,0x1a22d9f300f5138f,0xc63b05c06ca1493b,0x01ae3a4617c510ea]);
const BLS12_377_R: U384Repr = U384Repr([0x02cdffffffffff68,0x51409f837fffffb1,0x9f7db3a98a7d3ff2,0x7b4e97b76e7c6305,0x4cf495bf803c84e8,0x008d6661e2fdf49a]);
const BLS12_377_R2: U384Repr = U384Repr([0xb786686c9400cd22,0x0329fcaab00431b1,0x22a5f11162d6b46d,0xbfdf7d03827dc3ac,0x837e92f041790bf9,0x006dfccb1e914b88]);
const BLS12_377_MONT_INV: u64 = 0x8508bfffffffffff;

pub const BLS12_377_FIELD: PrimeField<U384Repr> = PrimeField::<U384Repr> {
    mont_power: 384,
    modulus_bits: 377,
    modulus: BLS12_377_MODULUS,
    mont_r: BLS12_377_R,
    mont_r2: BLS12_377_R2,
    mont_inv: BLS12_377_MONT_INV,  
};


const BLS12_377_FP_NON_RESIDUE_REPR: U384Repr = U384Repr([0xfc0b8000000002fa,0x97d39cf6e000018b,0x2072420fbfa05044,0xcbbcbd50d97c3802,0x0baf1ec35813f9eb,0x009974a2c0945ad2]);

const BLS12_377_FP_NON_RESIDUE: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_377_FP_NON_RESIDUE_REPR, 
    U384Repr,
    BLS12_377_FIELD
);

pub const BLS12_377_FP_ZERO: decl_fp!(U384Repr) = repr_into_fp!(
    REPR_ZERO, 
    U384Repr,
    BLS12_377_FIELD
);

pub const BLS12_377_FP_ONE: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_377_R, 
    U384Repr,
    BLS12_377_FIELD
);
 
const BLS12_377_EXTENSION_2_FROB_COEFF_0_REPR: U384Repr = U384Repr([0x02cdffffffffff68,0x51409f837fffffb1,0x9f7db3a98a7d3ff2,0x7b4e97b76e7c6305,0x4cf495bf803c84e8,0x008d6661e2fdf49a]);
const BLS12_377_EXTENSION_2_FROB_COEFF_1_REPR: U384Repr = U384Repr([0x823ac00000000099,0xc5cabdc0b000004f,0x7f75ae862f8c080d,0x9ed4423b9278b089,0x79467000ec64c452,0x0120d3e434c71c50]);

const BLS12_377_EXTENSION_2_FROB_COEFF_0: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_377_EXTENSION_2_FROB_COEFF_0_REPR, 
    U384Repr,
    BLS12_377_FIELD
);

const BLS12_377_EXTENSION_2_FROB_COEFF_1: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_377_EXTENSION_2_FROB_COEFF_1_REPR, 
    U384Repr,
    BLS12_377_FIELD
);


pub const BLS12_377_EXTENSION_2_FIELD: Extension2<'static, U384Repr, PrimeField<U384Repr>> = 
    Extension2::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_377_FIELD,
        non_residue: BLS12_377_FP_NON_RESIDUE,
        frobenius_coeffs_c1: [BLS12_377_EXTENSION_2_FROB_COEFF_0, BLS12_377_EXTENSION_2_FROB_COEFF_1],
        frobenius_coeffs_are_calculated: true
    };


const BLS12_377_FP2_NON_RESIDUE_C0_REPR: U384Repr = U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]);
const BLS12_377_FP2_NON_RESIDUE_C1_REPR: U384Repr = U384Repr([0x02cdffffffffff68,0x51409f837fffffb1,0x9f7db3a98a7d3ff2,0x7b4e97b76e7c6305,0x4cf495bf803c84e8,0x008d6661e2fdf49a]);

const BLS12_377_FP2_NON_RESIDUE_C0: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_377_FP2_NON_RESIDUE_C1_REPR, 
    U384Repr,
    BLS12_377_FIELD
);

const BLS12_377_FP2_NON_RESIDUE_C1: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_377_FP2_NON_RESIDUE_C1_REPR, 
    U384Repr,
    BLS12_377_FIELD
);

const BLS12_377_FP2_NON_RESIDUE: decl_fp2!(U384Repr) = repr_into_fp2!(
    BLS12_377_FP2_NON_RESIDUE_C0, 
    BLS12_377_FP2_NON_RESIDUE_C1,
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

pub const BLS12_377_SUBGROUP_ORDER: [u64; 4] = [
    0x0a11800000000001,
    0x59aa76fed0000001,
    0x60b44d1e5c37b001,
    0x12ab655e9a2ca556
];

const BLS12_377_X: [u64; 1] = [0x8508c00000000001];
const BLS12_377_X_IS_NEGATIVE: bool = false;

const BLS12_377_B_FOR_G1_REPR: U384Repr = U384Repr([0x862f3ffffffffd9f,0x2df720c9cffffec3,0x5f036c766febb7c9,0xd31784eab8fc7887,0x6d97513d9450ca66,0x00875f417432c17e]);
const BLS12_377_B_FOR_G1: Fp<'static, U384Repr, PrimeField<U384Repr>> = 
    Fp::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_377_FIELD,
        repr: BLS12_377_B_FOR_G1_REPR
    };    

const BLS12_377_B_FOR_G2_C0_REPR: U384Repr = U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]);
const BLS12_377_B_FOR_G2_C1_REPR: U384Repr = U384Repr([0x01c8999999999a14,0x37d5649a266666a6,0xff91586b593cd33e,0xe5769b62db93c06d,0x2dd1f333f0509d0e,0x00e70fe9c3d27d0d]);

const BLS12_377_B_FOR_G2_C0: Fp<'static, U384Repr, PrimeField<U384Repr>> = 
    Fp::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_377_FIELD,
        repr: BLS12_377_B_FOR_G2_C0_REPR
    };  

const BLS12_377_B_FOR_G2_C1: Fp<'static, U384Repr, PrimeField<U384Repr>> = 
    Fp::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_377_FIELD,
        repr: BLS12_377_B_FOR_G2_C1_REPR
    };    

const BLS12_377_B_FOR_G2: Fp2<'static, U384Repr, PrimeField<U384Repr>> = 
    Fp2::<'static, U384Repr, PrimeField<U384Repr>> {
        c0: BLS12_377_B_FOR_G2_C0,
        c1: BLS12_377_B_FOR_G2_C1,
        extension_field: &BLS12_377_EXTENSION_2_FIELD
    };


pub const BLS12_377_G1_CURVE_PARAMETERS: CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>> = 
    CurveOverFpParameters::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_377_FIELD
    };

pub const BLS12_377_G2_CURVE_PARAMETERS: CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>> = 
    CurveOverFp2Parameters::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_377_EXTENSION_2_FIELD
    };

pub const BLS12_377_FP2_ZERO: decl_fp2!(U384Repr) = repr_into_fp2!(
    BLS12_377_FP_ZERO, 
    BLS12_377_FP_ZERO,
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

pub const BLS12_377_FP2_ONE: decl_fp2!(U384Repr) = repr_into_fp2!(
    BLS12_377_FP_ONE, 
    BLS12_377_FP_ZERO,
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP6_FROB_C1_0: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x02cdffffffffff68,0x51409f837fffffb1,0x9f7db3a98a7d3ff2,0x7b4e97b76e7c6305,0x4cf495bf803c84e8,0x008d6661e2fdf49a]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP6_FROB_C1_1: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x1dab52f567cdf753,0x275ba68225d3fcbf,0xd479033c31671c95,0xbaf5e01f45f5eeb2,0xc124ac0d60ff519a,0x012e21e26d3f27d6]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0xe8175168f936f1a3,0xa1f9faad6e240880,0x0327ea0da32ce732,0xcbc495c1e2fe7cc2,0xbef5b5eec7980fde,0x00450861e49e5fb2]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP6_FROB_C1_2: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x2c766f925a7b8727,0x03d7f6b0253d58b5,0x838ec0deec122131,0xbd5eb3e9f658bb10,0x6942bd126ed3e52e,0x01673786dd04ed6a]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP6_FROB_C1_3: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0xe4a79fc4d4b5361f,0x3f77620e95796e81,0xc331b835dbe63ed4,0x7d9226378f678321,0xe8d237560a8dc3e4,0x008a6a848e739d69]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0xa03b8e62d0fd012f,0x93c031ea435b3d76,0xafc6f5dbe5871aef,0x7cc296b19f5040e7,0xcd2bc9301cdc60ac,0x01310f7ced1c32b7]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP6_FROB_C1_4: decl_fp2!(U384Repr) = BLS12_377_FP2_ZERO;
const BLS12_377_FP6_FROB_C1_5: decl_fp2!(U384Repr) = BLS12_377_FP2_ZERO;

const BLS12_377_FP6_FROB_C2_0: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x02cdffffffffff68,0x51409f837fffffb1,0x9f7db3a98a7d3ff2,0x7b4e97b76e7c6305,0x4cf495bf803c84e8,0x008d6661e2fdf49a]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP6_FROB_C2_1: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x88b5c17a7958b2fe,0x7f5e7e57ba07bc1c,0x6610c4919f5335d4,0x2662d66110d01894,0x528db91fde4912ee,0x01578d053afcb64e]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x7ce7510b8d5e6d45,0x0623598c82c4f51b,0xc595749a3ca4dfb7,0x4089827645e7035a,0x1d4fa3a2051729e1,0x0028dcfa7aab9072]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP6_FROB_C2_2: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0xdacd106da5847973,0xd8fe2454bac2a79a,0x1ada4fd6fd832edc,0xfb9868449d150908,0xd63eb8aeea32285e,0x0167d6a36f873fd0]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP6_FROB_C2_3: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x6693d2d815e1849f,0xb8317b4d8c16fa8f,0x417a052262ed3938,0xbb4a9f4398ca0e30,0xf0a0eca1f774c413,0x0103be639ab8b9b1]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x1e9a7f00446c4415,0xa2adab41fb145979,0x8974112b4fb7fd03,0x7ea657b93a6854e4,0xe5d71e58ba63a9d1,0x009d3c2719419801]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP6_FROB_C2_4: decl_fp2!(U384Repr) = BLS12_377_FP2_ZERO;
const BLS12_377_FP6_FROB_C2_5: decl_fp2!(U384Repr) = BLS12_377_FP2_ZERO;

pub const BLS12_377_EXTENSION_6_FIELD: Extension3Over2<'static, U384Repr, PrimeField<U384Repr>> = 
    Extension3Over2::<'static, U384Repr, PrimeField<U384Repr>> {
        non_residue: BLS12_377_FP2_NON_RESIDUE,
        field: &BLS12_377_EXTENSION_2_FIELD,
        frobenius_coeffs_c1: [BLS12_377_FP6_FROB_C1_0, BLS12_377_FP6_FROB_C1_1, BLS12_377_FP6_FROB_C1_2, BLS12_377_FP6_FROB_C1_3, BLS12_377_FP6_FROB_C1_4, BLS12_377_FP6_FROB_C1_5],
        frobenius_coeffs_c2: [BLS12_377_FP6_FROB_C2_0, BLS12_377_FP6_FROB_C2_1, BLS12_377_FP6_FROB_C2_2, BLS12_377_FP6_FROB_C2_3, BLS12_377_FP6_FROB_C2_4, BLS12_377_FP6_FROB_C2_5],
        frobenius_coeffs_are_calculated: true
    };


const BLS12_377_FP12_FROB_C1_0: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x02cdffffffffff68,0x51409f837fffffb1,0x9f7db3a98a7d3ff2,0x7b4e97b76e7c6305,0x4cf495bf803c84e8,0x008d6661e2fdf49a]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP12_FROB_C1_1: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x8e4881eda9715337,0x68933c1b14707f81,0xabe52b0749129354,0x039d6a8fb0acdc2e,0x77c4c0656b2dc79b,0x0102a8dbf4497b6a]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0xfe29bc6ac0baf987,0xcbd4ec5004bfc592,0x12f93242a0ff802b,0x6fbab242c943d0fe,0x0362fa6af0205e6e,0x013a4205b11a7684]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP12_FROB_C1_2: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0xdacd106da5847973,0xd8fe2454bac2a79a,0x1ada4fd6fd832edc,0xfb9868449d150908,0xd63eb8aeea32285e,0x0167d6a36f873fd0]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP12_FROB_C1_3: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0xc059cd145fe9bd8c,0xdb5711d0ca9b5c20,0xff5b91f6ceecf992,0x2c9e02f131e3fecb,0x6b1d1c9c56f21cd0,0x019d034f7fdacf4e]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0xe5bf2331ea1a4ad9,0xfa109dcb130047d8,0x5322840482579c81,0xe07b400cda4647b3,0xcdbfbc15d4645c07,0x008f9a1e7e02e793]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP12_FROB_C1_4: decl_fp2!(U384Repr) = BLS12_377_FP2_ZERO;
const BLS12_377_FP12_FROB_C1_5: decl_fp2!(U384Repr) = BLS12_377_FP2_ZERO;

const BLS12_377_FP12_FROB_C1_6: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x02cdffffffffff68,0x51409f837fffffb1,0x9f7db3a98a7d3ff2,0x7b4e97b76e7c6305,0x4cf495bf803c84e8,0x008d6661e2fdf49a]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_FP12_FROB_C1_7: decl_fp2!(U384Repr) = BLS12_377_FP2_ZERO;
const BLS12_377_FP12_FROB_C1_8: decl_fp2!(U384Repr) = BLS12_377_FP2_ZERO;
const BLS12_377_FP12_FROB_C1_9: decl_fp2!(U384Repr) = BLS12_377_FP2_ZERO;
const BLS12_377_FP12_FROB_C1_10: decl_fp2!(U384Repr) = BLS12_377_FP2_ZERO;
const BLS12_377_FP12_FROB_C1_11: decl_fp2!(U384Repr) = BLS12_377_FP2_ZERO;

const BLS12_377_FP6_ZERO: Fp6<'static, U384Repr, PrimeField<U384Repr>> = 
    Fp6::<'static, U384Repr, PrimeField<U384Repr>> {
        c0: BLS12_377_FP2_ZERO,
        c1: BLS12_377_FP2_ZERO,
        c2: BLS12_377_FP2_ZERO,
        extension_field: &BLS12_377_EXTENSION_6_FIELD
    };

pub const BLS12_377_EXTENSION_12_FIELD: Extension2Over3Over2<'static, U384Repr, PrimeField<U384Repr>> = 
Extension2Over3Over2::<'static, U384Repr, PrimeField<U384Repr>> {
    non_residue: BLS12_377_FP6_ZERO,
    field: &BLS12_377_EXTENSION_6_FIELD,
    frobenius_coeffs_c1: [
            BLS12_377_FP12_FROB_C1_0, BLS12_377_FP12_FROB_C1_1,
            BLS12_377_FP12_FROB_C1_2, BLS12_377_FP12_FROB_C1_3,
            BLS12_377_FP12_FROB_C1_4, BLS12_377_FP12_FROB_C1_5,
            BLS12_377_FP12_FROB_C1_6, BLS12_377_FP12_FROB_C1_7,
            BLS12_377_FP12_FROB_C1_8, BLS12_377_FP12_FROB_C1_9,
            BLS12_377_FP12_FROB_C1_10, BLS12_377_FP12_FROB_C1_11
        ],
    frobenius_coeffs_are_calculated: true
};   

pub const BLS12_377_G1_CURVE: WeierstrassCurve<'static, CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>> = 
    WeierstrassCurve::<'static, CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>> {
        a: BLS12_377_FP_ZERO,
        b: BLS12_377_B_FOR_G1,
        curve_type: CurveType::AIsZero,
        subgroup_order_repr: &BLS12_377_SUBGROUP_ORDER,
        params: &BLS12_377_G1_CURVE_PARAMETERS
    };   

pub const BLS12_377_G2_CURVE: WeierstrassCurve<'static, CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>> = 
    WeierstrassCurve::<'static, CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>> {
        a: BLS12_377_FP2_ZERO,
        b: BLS12_377_B_FOR_G2,
        curve_type: CurveType::AIsZero,
        subgroup_order_repr: &BLS12_377_SUBGROUP_ORDER,
        params: &BLS12_377_G2_CURVE_PARAMETERS
    };   

const BLS12_377_G1_GENERATOR_X: decl_fp!(U384Repr) = repr_into_fp!(
    U384Repr([0x260f33b9772451f4,0xc54dd773169d5658,0x5c1551c469a510dd,0x761662e4425e1698,0xc97d78cc6f065272,0x00a41206b361fd4d]), 
    U384Repr,
    BLS12_377_FIELD
);
    
const BLS12_377_G1_GENERATOR_Y: decl_fp!(U384Repr) = repr_into_fp!(
    U384Repr([0x8193961fb8cb81f3,0x00638d4c5f44adb8,0xfafaf3dad4daf54a,0xc27849e2d655cd18,0x2ec3ddb401d52814,0x007da93326303c71]), 
    U384Repr,
    BLS12_377_FIELD
);

const BLS12_377_G2_GENERATOR_X: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x68904082f268725b,0x668f2ea74f45328b,0xebca7a65802be84f,0x1e1850f4c1ada3e6,0x830dc22d588ef1e9,0x01862a81767c0982]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x5f02a915c91c7f39,0xf8c553ba388da2a7,0xd51a416dbd198850,0xe943c6f38ae3073a,0xffe24aa8259a4981,0x011853391e73dfdd]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

const BLS12_377_G2_GENERATOR_Y: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0xd5b19b897881430f,0x05be9118a5b371ed,0x6063f91f86c131ee,0x3244a61be8f4ec19,0xa02e425b9f9a3a12,0x018af8c04f3360d2]), 
        U384Repr,
        BLS12_377_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x57601ac71a5b96f5,0xe99acc1714f2440e,0x2339612f10118ea9,0x8321e68a3b1cd722,0x2b543b050cc74917,0x00590182b396c112]), 
        U384Repr,
        BLS12_377_FIELD
    ),
    U384Repr,
    BLS12_377_EXTENSION_2_FIELD
);

pub const BLS12_377_G1_GENERATOR: CurvePoint<'static, CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>> = 
    CurvePoint::<'static, CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>> 
    {
        curve: &BLS12_377_G1_CURVE,
        x: BLS12_377_G1_GENERATOR_X,
        y: BLS12_377_G1_GENERATOR_Y,
        z: BLS12_377_FP_ONE,
    };

pub const BLS12_377_G2_GENERATOR: CurvePoint<'static, CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>> = 
    CurvePoint::<'static, CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>>
    {
        curve: &BLS12_377_G2_CURVE,
        x: BLS12_377_G2_GENERATOR_X,
        y: BLS12_377_G2_GENERATOR_Y,
        z: BLS12_377_FP2_ONE,
    };

pub const BLS12_377_PAIRING_ENGINE: Bls12Instance<
    'static, 
    U384Repr, 
    PrimeField<U384Repr>, 
    CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>,
    CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>
> = Bls12Instance::<
    'static, 
    U384Repr, 
    PrimeField<U384Repr>, 
    CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>,
    CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>
> {
    x: &BLS12_377_X,
    x_is_negative: BLS12_377_X_IS_NEGATIVE,
    twist_type: TwistType::D,
    base_field: &BLS12_377_FIELD,
    curve: &BLS12_377_G1_CURVE,
    curve_twist: &BLS12_377_G2_CURVE,
    fp2_extension: &BLS12_377_EXTENSION_2_FIELD,
    fp6_extension: &BLS12_377_EXTENSION_6_FIELD,
    fp12_extension: &BLS12_377_EXTENSION_12_FIELD,
    prefer_naf: false,
    x_naf: Vec::new()
};

#[cfg(test)]
mod test {
    use crate::traits::FieldElement;
    use super::*;
    use crate::public_interface::decode_g1::serialize_g1_point;

    #[test]
    fn test_engine_bilinearity() {
        use crate::weierstrass::Group;
        use crate::pairings::PairingEngine;

        let p = BLS12_377_G1_GENERATOR.clone();
        let q = BLS12_377_G2_GENERATOR.clone();

        let mut p2 = p.mul(vec![12345678]);
        p2.normalize();

        let mut q2 = q.mul(vec![12345678]);
        q2.normalize();

        let ans1 = BLS12_377_PAIRING_ENGINE.pair(&[p.clone()], &[q2]).unwrap();
        let ans2 = BLS12_377_PAIRING_ENGINE.pair(&[p2], &[q.clone()]).unwrap();
        let ans3 = BLS12_377_PAIRING_ENGINE.pair(&[p], &[q]).unwrap();
        let ans3 = ans3.pow(&vec![12345678]);

        assert!(ans1 == ans2);
        assert!(ans1 == ans3);
    }

    fn output_test_vector(input: &[u8], output: &[u8]) {
        println!("Input: 0x{}", hex::encode(input));
        println!("Output: 0x{}", hex::encode(output));
    }

    fn serialize_scalar(repr: &[u64], len: usize) -> Vec<u8> {
        let a = MaxFieldUint::from(repr);

        let mut res = vec![0u8; a.as_ref().len() * 8];
        a.to_big_endian(&mut res);

        res.reverse();
        res.truncate(len);
        res.reverse();

        res
    }


    #[test]
    fn test_g1_mul_by_zero() {
        let mut input_encoding = serialize_g1_point(48, &BLS12_377_G1_GENERATOR).expect("must serialize a generator");
        let scalar = [0u64];
        let scalar_encoding = serialize_scalar(&scalar, 32);
        input_encoding.extend(scalar_encoding);

        let mut out = BLS12_377_G1_GENERATOR.mul(&scalar);
        out.normalize();

        let output_encoding = serialize_g1_point(48, &out).expect("must serialize a generator");

        output_test_vector(&input_encoding, &output_encoding);
    }
}
