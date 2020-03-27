use crate::field::*;

extern crate once_cell;

use self::once_cell::sync::Lazy;

use super::isogeny::*;
use super::simple_swu::*;
use crate::weierstrass::*;

pub static BLS12_G1_MAPPING_PARAMS: Lazy<
        (
            SwuParameters<CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>>,
            IsogenyParameters<CurveOverFpParameters<'static, U384Repr, PrimeField<U384Repr>>>
        )
        > = Lazy::new(|| {
            super::constants::calculate_bls12_381_g1_mapping_params(&crate::engines::bls12_381::BLS12_381_FIELD)
});

pub static BLS12_G2_MAPPING_PARAMS: Lazy<
        (
            SwuParameters<CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>>,
            IsogenyParameters<CurveOverFp2Parameters<'static, U384Repr, PrimeField<U384Repr>>>
        )
        > = Lazy::new(|| {
            super::constants::calculate_bls12_381_g2_mapping_params(&crate::engines::bls12_381::BLS12_381_EXTENSION_2_FIELD)
});