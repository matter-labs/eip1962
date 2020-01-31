use crate::weierstrass::curve::{WeierstrassCurve, CurvePoint};
use crate::field::{SizedPrimeField};
use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::weierstrass::CurveParameters;
use crate::constants::{MaxGroupSizeUint};

use super::decode_fp::*;

use super::decode_utils::{split, decode_group_order_with_length};

use crate::errors::ApiError;

pub(crate) fn parse_group_order_from_encoding<
    'a
    >(encoding: &'a [u8]) -> Result<(usize, MaxGroupSizeUint, &'a [u8]), ApiError>
{
    let ((order_len, order), rest) = decode_group_order_with_length(&encoding)?;
    if order.is_zero() {
        return Err(ApiError::InputError(format!("Group order is zero, file {}, line {}", file!(), line!())))
    }

    Ok((order_len, order, rest))
}

pub(crate) fn parse_ab_in_base_field_from_encoding<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>
    >(
        encoding: &'a [u8], 
        modulus_len: usize,
        base_field: &'a F
    ) -> Result<(Fp<'a, FE, F>, Fp<'a, FE, F>, &'a [u8]), ApiError>
{
    let (a, rest) = decode_fp(&encoding, modulus_len, base_field)?;
    let (b, rest) = decode_fp(&rest, modulus_len, base_field)?;

    Ok((a, b, rest))
}

pub(crate) fn serialize_g1_point<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE> + 'a,
    C: CurveParameters<BaseFieldElement = Fp<'a, FE, F>>
    >
    (
        modulus_len: usize,
        point: &CurvePoint<'a, C>
    ) -> Result<Vec<u8>, ApiError>
{
    let (x, y) = point.into_xy();
    let mut result = Vec::with_capacity(2*modulus_len);
    result.extend(serialize_fp_fixed_len(modulus_len, &x)?);
    result.extend(serialize_fp_fixed_len(modulus_len, &y)?);

    Ok(result)
}

pub(crate) fn decode_g1_point_from_xy<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE> + 'a,
    C: CurveParameters<BaseFieldElement = Fp<'a, FE, F>>
    >
    (
        bytes: &'a [u8], 
        field_byte_len: usize,
        curve: &'a WeierstrassCurve<'a, C>
    ) -> Result<(CurvePoint<'a, C>, &'a [u8]), ApiError>
{
    let (x, rest) = decode_fp(&bytes, field_byte_len, curve.params.params())?;
    let (y, rest) = decode_fp(&rest, field_byte_len, curve.params.params())?;
    
    let p: CurvePoint<'a, C> = CurvePoint::point_from_xy(&curve, x, y);
    
    Ok((p, rest))
}

pub(crate) fn decode_scalar_representation<
    'a
    >
    (
        bytes: &'a [u8], 
        order_byte_len: usize,
        _order: &MaxGroupSizeUint,
    ) -> Result<(MaxGroupSizeUint, &'a [u8]), ApiError>
{
    let (encoding, rest) = split(bytes, order_byte_len, "Input is not long enough to get scalar")?;
    let scalar = MaxGroupSizeUint::from_big_endian(&encoding);
    if &scalar > _order {
        return Err(ApiError::InputError(format!("Scalar is larger than the group order, file {}, line {}", file!(), line!())));
    }

    Ok((scalar, rest))
}

