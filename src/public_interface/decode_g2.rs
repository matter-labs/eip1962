use crate::weierstrass::Group;
use crate::weierstrass::twist;
use crate::weierstrass::cubic_twist;
use crate::field::{SizedPrimeField, field_from_modulus};
use crate::fp::Fp;
use crate::representation::ElementRepr;

use num_bigint::BigUint;
use num_traits::{Zero};

use super::decode_fp::*;

pub(crate) fn decode_g2_point_from_xy_in_fp2<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>,
    GE: ElementRepr,
    G: SizedPrimeField<Repr = GE>
    >
    (
        bytes: &'a [u8], 
        field_byte_len: usize,
        curve: &'a twist::WeierstrassCurveTwist<'a, FE, F, GE, G>
    ) -> Result<(twist::TwistPoint<'a, FE, F, GE, G>, &'a [u8]), ()>
{
    let (x, rest) = decode_fp2(&bytes, field_byte_len, curve.base_field)?;
    let (y, rest) = decode_fp2(&rest, field_byte_len, curve.base_field)?;
    
    let p: twist::TwistPoint<'a, FE, F, GE, G> = twist::TwistPoint::point_from_xy(&curve, x, y);
    
    Ok((p, rest))
}

pub(crate) fn decode_g2_point_from_xy_in_fp3<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>,
    GE: ElementRepr,
    G: SizedPrimeField<Repr = GE>
    >
    (
        bytes: &'a [u8], 
        field_byte_len: usize,
        curve: &'a cubic_twist::WeierstrassCurveTwist<'a, FE, F, GE, G>
    ) -> Result<(cubic_twist::TwistPoint<'a, FE, F, GE, G>, &'a [u8]), ()>
{
    let (x, rest) = decode_fp3(&bytes, field_byte_len, curve.base_field)?;
    let (y, rest) = decode_fp3(&rest, field_byte_len, curve.base_field)?;
    
    let p: cubic_twist::TwistPoint<'a, FE, F, GE, G> = cubic_twist::TwistPoint::point_from_xy(&curve, x, y);
    
    Ok((p, rest))
}