use crate::field::*;
use crate::traits::*;
use crate::fp::*;
use crate::extension_towers::fp2::*;
use crate::extension_towers::fp6_as_3_over_2::*;
use crate::extension_towers::fp12_as_2_over3_over_2::*;
use crate::weierstrass::*;
use crate::weierstrass::curve::*;
use crate::pairings::bls12::*;
use crate::pairings::TwistType;
use crate::integers::MaxFieldUint;

struct Bls12_381Extension2;

impl FieldExtension for Bls12_381Extension2 {
    const EXTENSION_DEGREE: usize = 2;
    
    type Element = Fp<'static, U384Repr, PrimeField<U384Repr>>;

    #[inline(always)]
    fn multiply_by_non_residue(&self, el: &mut Self::Element) {
        // non-residue is equal to -1, so we just negate
        el.negate();
    }
}

struct Bls12_381Extension6;

impl FieldExtension for Bls12_381Extension6 {
    const EXTENSION_DEGREE: usize = 3;
    
    type Element = Fp2<'static, U384Repr, PrimeField<U384Repr>>;

    #[inline(always)]
    fn multiply_by_non_residue(&self, el: &mut Self::Element) {
        // manually unroll multiplication by (1, 1)
        let v0 = el.c0.clone();
        let mut v1 = el.c1.clone();

        el.c1.add_assign(&el.c0);

        // non-residue c0 + c1
        let mut t0 = Fp::one(el.extension_field.field);
        t0.double();

        el.c1.mul_assign(&t0);
        el.c1.sub_assign(&v0);
        el.c1.sub_assign(&v1);
        el.c0 = v0;

        v1.mul_by_nonresidue(&Bls12_381Extension2);
        el.c0.add_assign(&v1);
   
    }
}

const REPR_ZERO: U384Repr = U384Repr([0,0,0,0,0,0]);

pub const BLS12_381_MODULUS_UINT: MaxFieldUint = MaxFieldUint::from_limbs(
    [
        0xb9feffffffffaaab,0x1eabfffeb153ffff,0x6730d2a0f6b0f624,0x64774b84f38512bf,
        0x4b1ba7b6434bacd7,0x1a0111ea397fe69a, 0x0, 0x0,
        0x0, 0x0, 0x0, 0x0,
        0x0, 0x0, 0x0, 0x0
    ]
);

const BLS12_381_MODULUS: U384Repr = U384Repr([0xb9feffffffffaaab,0x1eabfffeb153ffff,0x6730d2a0f6b0f624,0x64774b84f38512bf,0x4b1ba7b6434bacd7,0x1a0111ea397fe69a]);
const BLS12_381_R: U384Repr = U384Repr([0x760900000002fffd,0xebf4000bc40c0002,0x5f48985753c758ba,0x77ce585370525745,0x5c071a97a256ec6d,0x15f65ec3fa80e493]);
const BLS12_381_R2: U384Repr = U384Repr([0xf4df1f341c341746,0x0a76e6a609d104f1,0x8de5476c4c95b6d5,0x67eb88a9939d83c0,0x9a793e85b519952d,0x11988fe592cae3aa]);
const BLS12_381_MONT_INV: u64 = 0x89f3fffcfffcfffd;

const BLS12_381_FIELD: PrimeField<U384Repr> = PrimeField::<U384Repr> {
    mont_power: 384,
    modulus_bits: 381,
    modulus: BLS12_381_MODULUS,
    mont_r: BLS12_381_R,
    mont_r2: BLS12_381_R2,
    mont_inv: BLS12_381_MONT_INV,  
};

// macro_rules! decl_repr_into_fp {
//     ($repr:ident, $repr_type:ty, $field:ident) => ( Fp::<'static, $repr_type, PrimeField<$repr_type>> = Fp::<'static, $repr_type, PrimeField<$repr_type>> {
//             field: &$field,
//             repr: $repr
//         };
//     );
// }

macro_rules! decl_fp {
    ($repr_type:ty) => ( Fp::<'static, $repr_type, PrimeField<$repr_type>> );
}

// macro_rules! repr_into_fp {
//     ($repr:ident, $repr_type:ty, $field:ident) => ({
//         Fp::<'static, $repr_type, PrimeField<$repr_type>> {
//             field: &$field,
//             repr: $repr
//         }
//     });
// }

macro_rules! repr_into_fp {
    ($repr:expr, $repr_type:ty, $field:ident) => ({
        Fp::<'static, $repr_type, PrimeField<$repr_type>> {
            field: &$field,
            repr: $repr
        }
    });
}

macro_rules! decl_fp2 {
    ($repr_type:ty) => ( Fp2::<'static, $repr_type, PrimeField<$repr_type>> );
}

// macro_rules! repr_into_fp2 {
//     ($repr_c0:ident, $repr_c1:ident, $repr_type:ty, $field:ident) => ({
//         Fp2::<'static, $repr_type, PrimeField<$repr_type>> {
//             c0: $repr_c0,
//             c1: $repr_c1,
//             extension_field: &$field
//         }
//     });
// }

macro_rules! repr_into_fp2 {
    ($repr_c0:expr, $repr_c1:expr, $repr_type:ty, $field:ident) => ({
        Fp2::<'static, $repr_type, PrimeField<$repr_type>> {
            c0: $repr_c0,
            c1: $repr_c1,
            extension_field: &$field
        }
    });
}

const BLS12_381_FP_NON_RESIDUE_REPR: U384Repr = U384Repr([0x43f5fffffffcaaae,0x32b7fff2ed47fffd,0x07e83a49a2e99d69,0xeca8f3318332bb7a,0xef148d1ea0f4c069,0x040ab3263eff0206]);

const BLS12_381_FP_NON_RESIDUE: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_381_FP_NON_RESIDUE_REPR, 
    U384Repr,
    BLS12_381_FIELD
);

const BLS12_381_FP_ZERO: decl_fp!(U384Repr) = repr_into_fp!(
    REPR_ZERO, 
    U384Repr,
    BLS12_381_FIELD
);

const BLS12_381_FP_ONE: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_381_R, 
    U384Repr,
    BLS12_381_FIELD
);

// const BLS12_381_FP_NON_RESIDUE: Fp<'static, U384Repr, PrimeField<U384Repr>> = repr_into_fp!(
//     BLS12_381_FP_NON_RESIDUE_REPR, 
//     U384Repr,
//     BLS12_381_FIELD
// );

// const BLS12_381_FP_NON_RESIDUE: Fp<'static, U384Repr, PrimeField<U384Repr>> = 
//     Fp::<'static, U384Repr, PrimeField<U384Repr>> {
//         field: &BLS12_381_FIELD,
//         repr: BLS12_381_FP_NON_RESIDUE_REPR
//     };

// const BLS12_381_FP_ZERO: Fp<'static, U384Repr, PrimeField<U384Repr>> = 
//     Fp::<'static, U384Repr, PrimeField<U384Repr>> {
//         field: &BLS12_381_FIELD,
//         repr: REPR_ZERO
//     };    

const BLS12_381_EXTENSION_2_FROB_COEFF_0_REPR: U384Repr = U384Repr([0x760900000002fffd,0xebf4000bc40c0002,0x5f48985753c758ba,0x77ce585370525745,0x5c071a97a256ec6d,0x15f65ec3fa80e493]);
const BLS12_381_EXTENSION_2_FROB_COEFF_1_REPR: U384Repr = U384Repr([0x43f5fffffffcaaae,0x32b7fff2ed47fffd,0x07e83a49a2e99d69,0xeca8f3318332bb7a,0xef148d1ea0f4c069,0x040ab3263eff0206]);

const BLS12_381_EXTENSION_2_FROB_COEFF_0: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_381_EXTENSION_2_FROB_COEFF_0_REPR, 
    U384Repr,
    BLS12_381_FIELD
);

const BLS12_381_EXTENSION_2_FROB_COEFF_1: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_381_EXTENSION_2_FROB_COEFF_1_REPR, 
    U384Repr,
    BLS12_381_FIELD
);


pub const BLS12_381_EXTENSION_2_FIELD: Extension2<'static, U384Repr, PrimeField<U384Repr>> = 
    Extension2::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_381_FIELD,
        non_residue: BLS12_381_FP_NON_RESIDUE,
        frobenius_coeffs_c1: [BLS12_381_EXTENSION_2_FROB_COEFF_0, BLS12_381_EXTENSION_2_FROB_COEFF_1],
        frobenius_coeffs_are_calculated: true
    };


const BLS12_381_FP2_NON_RESIDUE_C0_REPR: U384Repr = U384Repr([0x760900000002fffd,0xebf4000bc40c0002,0x5f48985753c758ba,0x77ce585370525745,0x5c071a97a256ec6d,0x15f65ec3fa80e493]);
const BLS12_381_FP2_NON_RESIDUE_C1_REPR: U384Repr = U384Repr([0x760900000002fffd,0xebf4000bc40c0002,0x5f48985753c758ba,0x77ce585370525745,0x5c071a97a256ec6d,0x15f65ec3fa80e493]);

const BLS12_381_FP2_NON_RESIDUE_C0: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_381_FP2_NON_RESIDUE_C1_REPR, 
    U384Repr,
    BLS12_381_FIELD
);

const BLS12_381_FP2_NON_RESIDUE_C1: decl_fp!(U384Repr) = repr_into_fp!(
    BLS12_381_FP2_NON_RESIDUE_C1_REPR, 
    U384Repr,
    BLS12_381_FIELD
);

const BLS12_381_FP2_NON_RESIDUE: decl_fp2!(U384Repr) = repr_into_fp2!(
    BLS12_381_FP2_NON_RESIDUE_C0, 
    BLS12_381_FP2_NON_RESIDUE_C1,
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

// const BLS12_381_FP2_NON_RESIDUE_C0: Fp<'static, U384Repr, PrimeField<U384Repr>> = 
// Fp::<'static, U384Repr, PrimeField<U384Repr>> {
//     field: &BLS12_381_FIELD,
//     repr: BLS12_381_FP2_NON_RESIDUE_C0_REPR
// };

// const BLS12_381_FP2_NON_RESIDUE_C1:Fp<'static, U384Repr, PrimeField<U384Repr>> = 
// Fp::<'static, U384Repr, PrimeField<U384Repr>> {
//     field: &BLS12_381_FIELD,
//     repr: BLS12_381_FP2_NON_RESIDUE_C1_REPR
// };

// const BLS12_381_FP2_NON_RESIDUE: Fp2<'static, U384Repr, PrimeField<U384Repr>> = 
//     Fp2::<'static, U384Repr, PrimeField<U384Repr>> {
//         c0: BLS12_381_FP2_NON_RESIDUE_C0,
//         c1: BLS12_381_FP2_NON_RESIDUE_C1,
//         extension_field: &BLS12_381_EXTENSION_2_FIELD
//     };

// const BLS12_381_EXTENSION_2_FIELD: Bls12_381Extension2 = Bls12_381Extension2;

const BLS12_381_SUBGROUP_ORDER: [u64; 4] = [
    0xffffffff00000001,
    0x53bda402fffe5bfe,
    0x3339d80809a1d805,
    0x73eda753299d7d48
];

const BLS12_381_X: [u64; 1] = [0xd201000000010000];
const BLS12_381_X_IS_NEGATIVE: bool = true;

const BLS12_381_B_FOR_G1_REPR: U384Repr = U384Repr([0xaa270000000cfff3,0x53cc0032fc34000a,0x478fe97a6b0a807f,0xb1d37ebee6ba24d7,0x8ec9733bbf78ab2f,0x09d645513d83de7e]);
const BLS12_381_B_FOR_G1: Fp<'static, U384Repr, PrimeField<U384Repr>> = 
    Fp::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_381_FIELD,
        repr: BLS12_381_B_FOR_G1_REPR
    };    

const BLS12_381_B_FOR_G2_C0_REPR: U384Repr = U384Repr([0xaa270000000cfff3,0x53cc0032fc34000a,0x478fe97a6b0a807f,0xb1d37ebee6ba24d7,0x8ec9733bbf78ab2f,0x09d645513d83de7e]);
const BLS12_381_B_FOR_G2_C1_REPR: U384Repr = U384Repr([0xaa270000000cfff3,0x53cc0032fc34000a,0x478fe97a6b0a807f,0xb1d37ebee6ba24d7,0x8ec9733bbf78ab2f,0x09d645513d83de7e]);

const BLS12_381_B_FOR_G2_C0: Fp<'static, U384Repr, PrimeField<U384Repr>> = 
    Fp::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_381_FIELD,
        repr: BLS12_381_B_FOR_G2_C0_REPR
    };  

const BLS12_381_B_FOR_G2_C1: Fp<'static, U384Repr, PrimeField<U384Repr>> = 
    Fp::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_381_FIELD,
        repr: BLS12_381_B_FOR_G2_C1_REPR
    };    

const BLS12_381_B_FOR_G2: Fp2<'static, U384Repr, PrimeField<U384Repr>> = 
    Fp2::<'static, U384Repr, PrimeField<U384Repr>> {
        c0: BLS12_381_B_FOR_G2_C0,
        c1: BLS12_381_B_FOR_G2_C1,
        extension_field: &BLS12_381_EXTENSION_2_FIELD
    };


const BLS12_381_G1_CURVE_PARAMETERS: CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>> = 
    CurveOverFpParameters::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_381_FIELD
    };

const BLS12_381_G2_CURVE_PARAMETERS: CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>> = 
    CurveOverFp2Parameters::<'static, U384Repr, PrimeField<U384Repr>> {
        field: &BLS12_381_EXTENSION_2_FIELD
    };

// const BLS12_381_FP2_ZERO: Fp2<'static, U384Repr, PrimeField<U384Repr>> = 
//     Fp2::<'static, U384Repr, PrimeField<U384Repr>> {
//         c0: BLS12_381_FP_ZERO,
//         c1: BLS12_381_FP_ZERO,
//         extension_field: &BLS12_381_EXTENSION_2_FIELD
//     };

const BLS12_381_FP2_ZERO: decl_fp2!(U384Repr) = repr_into_fp2!(
    BLS12_381_FP_ZERO, 
    BLS12_381_FP_ZERO,
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP2_ONE: decl_fp2!(U384Repr) = repr_into_fp2!(
    BLS12_381_FP_ONE, 
    BLS12_381_FP_ZERO,
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

// const BLS12_381_FP2_NON_RESIDUE_C0_REPR: U384Repr = U384Repr([0x760900000002fffd,0xebf4000bc40c0002,0x5f48985753c758ba,0x77ce585370525745,0x5c071a97a256ec6d,0x15f65ec3fa80e493]);
// const BLS12_381_FP2_NON_RESIDUE_C1_REPR: U384Repr = U384Repr([0x760900000002fffd,0xebf4000bc40c0002,0x5f48985753c758ba,0x77ce585370525745,0x5c071a97a256ec6d,0x15f65ec3fa80e493]);

// const BLS12_381_FP2_NON_RESIDUE_C0: decl_fp!(U384Repr) = repr_into_fp!(
//     BLS12_381_FP2_NON_RESIDUE_C1_REPR, 
//     U384Repr,
//     BLS12_381_FIELD
// );

// const BLS12_381_FP2_NON_RESIDUE_C1: decl_fp!(U384Repr) = repr_into_fp!(
//     BLS12_381_FP2_NON_RESIDUE_C1_REPR, 
//     U384Repr,
//     BLS12_381_FIELD
// );

const BLS12_381_FP6_FROB_C1_0: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x760900000002fffd,0xebf4000bc40c0002,0x5f48985753c758ba,0x77ce585370525745,0x5c071a97a256ec6d,0x15f65ec3fa80e493]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP6_FROB_C1_1: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0xcd03c9e48671f071,0x5dab22461fcda5d2,0x587042afd3851b95,0x8eb60ebe01bacb9e,0x03f97d6e83d050d2,0x18f0206554638741]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP6_FROB_C1_2: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x30f1361b798a64e8,0xf3b8ddab7ece5a2a,0x16a8ca3ac61577f7,0xc26a2ff874fd029b,0x3636b76660701c6e,0x051ba4ab241b6160]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP6_FROB_C1_3: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x760900000002fffd,0xebf4000bc40c0002,0x5f48985753c758ba,0x77ce585370525745,0x5c071a97a256ec6d,0x15f65ec3fa80e493]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP6_FROB_C1_4: decl_fp2!(U384Repr) = BLS12_381_FP2_ZERO;
const BLS12_381_FP6_FROB_C1_5: decl_fp2!(U384Repr) = BLS12_381_FP2_ZERO;

const BLS12_381_FP6_FROB_C2_0: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x760900000002fffd,0xebf4000bc40c0002,0x5f48985753c758ba,0x77ce585370525745,0x5c071a97a256ec6d,0x15f65ec3fa80e493]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP6_FROB_C2_1: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x890dc9e4867545c3,0x2af322533285a5d5,0x50880866309b7e2c,0xa20d1b8c7e881024,0x14e4f04fe2db9068,0x14e56d3f1564853a]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP6_FROB_C2_2: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0xcd03c9e48671f071,0x5dab22461fcda5d2,0x587042afd3851b95,0x8eb60ebe01bacb9e,0x03f97d6e83d050d2,0x18f0206554638741]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP6_FROB_C2_3: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x43f5fffffffcaaae,0x32b7fff2ed47fffd,0x07e83a49a2e99d69,0xeca8f3318332bb7a,0xef148d1ea0f4c069,0x040ab3263eff0206]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP6_FROB_C2_4: decl_fp2!(U384Repr) = BLS12_381_FP2_ZERO;
const BLS12_381_FP6_FROB_C2_5: decl_fp2!(U384Repr) = BLS12_381_FP2_ZERO;

pub const BLS12_381_EXTENSION_6_FIELD: Extension3Over2<'static, U384Repr, PrimeField<U384Repr>> = 
    Extension3Over2::<'static, U384Repr, PrimeField<U384Repr>> {
        non_residue: BLS12_381_FP2_NON_RESIDUE,
        field: &BLS12_381_EXTENSION_2_FIELD,
        frobenius_coeffs_c1: [BLS12_381_FP6_FROB_C1_0, BLS12_381_FP6_FROB_C1_1, BLS12_381_FP6_FROB_C1_2, BLS12_381_FP6_FROB_C1_3, BLS12_381_FP6_FROB_C1_4, BLS12_381_FP6_FROB_C1_5],
        frobenius_coeffs_c2: [BLS12_381_FP6_FROB_C2_0, BLS12_381_FP6_FROB_C2_1, BLS12_381_FP6_FROB_C2_2, BLS12_381_FP6_FROB_C2_3, BLS12_381_FP6_FROB_C2_4, BLS12_381_FP6_FROB_C2_5],
        frobenius_coeffs_are_calculated: true
    };


const BLS12_381_FP12_FROB_C1_0: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x760900000002fffd,0xebf4000bc40c0002,0x5f48985753c758ba,0x77ce585370525745,0x5c071a97a256ec6d,0x15f65ec3fa80e493]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP12_FROB_C1_1: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x07089552b319d465,0xc6695f92b50a8313,0x97e83cccd117228f,0xa35baecab2dc29ee,0x1ce393ea5daace4d,0x08f2220fb0fb66eb]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0xb2f66aad4ce5d646,0x5842a06bfc497cec,0xcf4895d42599d394,0xc11b9cba40a8e8d0,0x2e3813cbe5a0de89,0x110eefda88847faf]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP12_FROB_C1_2: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0xecfb361b798dba3a,0xc100ddb891865a2c,0x0ec08ff1232bda8e,0xd5c13cc6f1ca4721,0x47222a47bf7b5c04,0x0110f184e51c5f59]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP12_FROB_C1_3: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x3e2f585da55c9ad1,0x4294213d86c18183,0x382844c88b623732,0x92ad2afd19103e18,0x1d794e4fac7cf0b9,0x0bd592fc7d825ec8]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x7bcfa7a25aa30fda,0xdc17dec12a927e7c,0x2f088dd86b4ebef1,0xd1ca2087da74d4a7,0x2da2596696cebc1d,0x0e2b7eedbbfd87d2]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP12_FROB_C1_4: decl_fp2!(U384Repr) = BLS12_381_FP2_ZERO;
const BLS12_381_FP12_FROB_C1_5: decl_fp2!(U384Repr) = BLS12_381_FP2_ZERO;

const BLS12_381_FP12_FROB_C1_6: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x43f5fffffffcaaae,0x32b7fff2ed47fffd,0x07e83a49a2e99d69,0xeca8f3318332bb7a,0xef148d1ea0f4c069,0x040ab3263eff0206]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000,0x0000000000000000]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_FP12_FROB_C1_7: decl_fp2!(U384Repr) = BLS12_381_FP2_ZERO;
const BLS12_381_FP12_FROB_C1_8: decl_fp2!(U384Repr) = BLS12_381_FP2_ZERO;
const BLS12_381_FP12_FROB_C1_9: decl_fp2!(U384Repr) = BLS12_381_FP2_ZERO;
const BLS12_381_FP12_FROB_C1_10: decl_fp2!(U384Repr) = BLS12_381_FP2_ZERO;
const BLS12_381_FP12_FROB_C1_11: decl_fp2!(U384Repr) = BLS12_381_FP2_ZERO;

const BLS12_381_FP6_ZERO: Fp6<'static, U384Repr, PrimeField<U384Repr>> = 
    Fp6::<'static, U384Repr, PrimeField<U384Repr>> {
        c0: BLS12_381_FP2_ZERO,
        c1: BLS12_381_FP2_ZERO,
        c2: BLS12_381_FP2_ZERO,
        extension_field: &BLS12_381_EXTENSION_6_FIELD
    };



pub const BLS12_381_EXTENSION_12_FIELD: Extension2Over3Over2<'static, U384Repr, PrimeField<U384Repr>> = 
Extension2Over3Over2::<'static, U384Repr, PrimeField<U384Repr>> {
    non_residue: BLS12_381_FP6_ZERO,
    field: &BLS12_381_EXTENSION_6_FIELD,
    frobenius_coeffs_c1: [
            BLS12_381_FP12_FROB_C1_0, BLS12_381_FP12_FROB_C1_1,
            BLS12_381_FP12_FROB_C1_2, BLS12_381_FP12_FROB_C1_3,
            BLS12_381_FP12_FROB_C1_4, BLS12_381_FP12_FROB_C1_5,
            BLS12_381_FP12_FROB_C1_6, BLS12_381_FP12_FROB_C1_7,
            BLS12_381_FP12_FROB_C1_8, BLS12_381_FP12_FROB_C1_9,
            BLS12_381_FP12_FROB_C1_10, BLS12_381_FP12_FROB_C1_11
        ],
    frobenius_coeffs_are_calculated: true
};   

const BLS12_381_G1_CURVE: WeierstrassCurve<'static, CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>> = 
    WeierstrassCurve::<'static, CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>> {
        a: BLS12_381_FP_ZERO,
        b: BLS12_381_B_FOR_G1,
        curve_type: CurveType::AIsZero,
        subgroup_order_repr: &BLS12_381_SUBGROUP_ORDER,
        params: &BLS12_381_G1_CURVE_PARAMETERS
    };   

const BLS12_381_G2_CURVE: WeierstrassCurve<'static, CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>> = 
    WeierstrassCurve::<'static, CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>> {
        a: BLS12_381_FP2_ZERO,
        b: BLS12_381_B_FOR_G2,
        curve_type: CurveType::AIsZero,
        subgroup_order_repr: &BLS12_381_SUBGROUP_ORDER,
        params: &BLS12_381_G2_CURVE_PARAMETERS
    };   

const BLS12_381_G1_GENERATOR_X: decl_fp!(U384Repr) = repr_into_fp!(
    U384Repr([0x5cb38790fd530c16,0x7817fc679976fff5,0x154f95c7143ba1c1,0xf0ae6acdf3d0e747,0xedce6ecc21dbf440,0x120177419e0bfb75]), 
    U384Repr,
    BLS12_381_FIELD
);
    
const BLS12_381_G1_GENERATOR_Y: decl_fp!(U384Repr) = repr_into_fp!(
    U384Repr([0xbaac93d50ce72271,0x8c22631a7918fd8e,0xdd595f13570725ce,0x51ac582950405194,0x0e1c8c3fad0059c0,0x0bbc3efc5008a26a]), 
    U384Repr,
    BLS12_381_FIELD
);

const BLS12_381_G2_GENERATOR_X: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0xf5f28fa202940a10,0xb3f5fb2687b4961a,0xa1a893b53e2ae580,0x9894999d1a3caee9,0x6f67b7631863366b,0x058191924350bcd7]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0xa5a9c0759e23f606,0xaaa0c59dbccd60c3,0x3bb17e18e2867806,0x1b1ab6cc8541b367,0xc2b6ed0ef2158547,0x11922a097360edf3]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_G2_GENERATOR_Y: decl_fp2!(U384Repr) = repr_into_fp2!(
    repr_into_fp!(
        U384Repr([0x4c730af860494c4a,0x597cfa1f5e369c5a,0xe7e6856caa0a635a,0xbbefb5e96e0d495f,0x07d3a975f0ef25a2,0x0083fd8e7e80dae5]), 
        U384Repr,
        BLS12_381_FIELD
    ), 
    repr_into_fp!(
        U384Repr([0xadc0fc92df64b05d,0x18aa270a2b1461dc,0x86adac6a3be4eba0,0x79495c4ec93da33a,0xe7175850a43ccaed,0x0b2bc2a163de1bf2]), 
        U384Repr,
        BLS12_381_FIELD
    ),
    U384Repr,
    BLS12_381_EXTENSION_2_FIELD
);

const BLS12_381_G1_GENERATOR: CurvePoint<'static, CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>> = 
    CurvePoint::<'static, CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>> 
    {
        curve: &BLS12_381_G1_CURVE,
        x: BLS12_381_G1_GENERATOR_X,
        y: BLS12_381_G1_GENERATOR_Y,
        z: BLS12_381_FP_ONE,
    };

const BLS12_381_G2_GENERATOR: CurvePoint<'static, CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>> = 
    CurvePoint::<'static, CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>>
    {
        curve: &BLS12_381_G2_CURVE,
        x: BLS12_381_G2_GENERATOR_X,
        y: BLS12_381_G2_GENERATOR_Y,
        z: BLS12_381_FP2_ONE,
    };

const BLS12_381_PAIRING_ENGINE: Bls12Instance<
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
    x: &BLS12_381_X,
    x_is_negative: BLS12_381_X_IS_NEGATIVE,
    twist_type: TwistType::M,
    base_field: &BLS12_381_FIELD,
    curve: &BLS12_381_G1_CURVE,
    curve_twist: &BLS12_381_G2_CURVE,
    fp2_extension: &BLS12_381_EXTENSION_2_FIELD,
    fp6_extension: &BLS12_381_EXTENSION_6_FIELD,
    fp12_extension: &BLS12_381_EXTENSION_12_FIELD,
    prefer_naf: false,
    x_naf: Vec::new()
};

#[cfg(test)]
mod test {
    use super::*;
    use crate::public_interface::decode_g1::serialize_g1_point;

    #[test]
    fn test_engine_bilinearity() {
        use crate::weierstrass::Group;
        use crate::pairings::PairingEngine;

        let p = BLS12_381_G1_GENERATOR.clone();
        let q = BLS12_381_G2_GENERATOR.clone();

        let mut p2 = p.mul(vec![12345678]);
        p2.normalize();

        let mut q2 = q.mul(vec![12345678]);
        q2.normalize();

        let ans1 = BLS12_381_PAIRING_ENGINE.pair(&[p.clone()], &[q2]).unwrap();
        let ans2 = BLS12_381_PAIRING_ENGINE.pair(&[p2], &[q.clone()]).unwrap();
        let ans3 = BLS12_381_PAIRING_ENGINE.pair(&[p], &[q]).unwrap();
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
        let mut input_encoding = serialize_g1_point(48, &BLS12_381_G1_GENERATOR).expect("must serialize a generator");
        let scalar = [0u64];
        let scalar_encoding = serialize_scalar(&scalar, 32);
        input_encoding.extend(scalar_encoding);

        let mut out = BLS12_381_G1_GENERATOR.mul(&scalar);
        out.normalize();

        let output_encoding = serialize_g1_point(48, &out).expect("must serialize a generator");

        output_test_vector(&input_encoding, &output_encoding);
    }
}
