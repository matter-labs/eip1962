use crate::weierstrass::curve::{WeierstrassCurve, CurvePoint};
use crate::field::{SizedPrimeField, field_from_modulus, PrimeField};
use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::traits::ZeroAndOne;
use crate::weierstrass::CurveParameters;

use super::constants::*;
use super::decode_fp::*;

use num_bigint::BigUint;
use num_traits::{Zero};

use super::decode_utils::split;

use crate::errors::ApiError;

pub(crate) fn parse_base_field_from_encoding<
    'a,
    FE: ElementRepr,
    >(encoding: &'a [u8]) -> Result<(PrimeField<FE>, usize, BigUint, &'a [u8]), ApiError>
{
    let ((modulus, modulus_len), rest) = get_base_field_params(&encoding)?;
    let field = field_from_modulus::<FE>(modulus.clone()).map_err(|_| {
        ApiError::InputError("Failed to create prime field from modulus".to_owned())
    })?;
    if rest.len() < modulus_len {
        return Err(ApiError::InputError("Input is not long enough".to_owned()));
    }

    Ok((field, modulus_len, modulus, rest))
}

pub(crate) fn parse_group_order_from_encoding<
    'a
    >(encoding: &'a [u8]) -> Result<(Vec<u64>, usize, BigUint, &'a [u8]), ApiError>
{
    use crate::field::biguint_to_u64_vec;
    let ((order, order_len), rest) = get_g1_curve_params(&encoding)?;
    let order = BigUint::from_bytes_be(&order);
    if order.is_zero() {
        return Err(ApiError::InputError(format!("Group order is zero, file {}, line {}", file!(), line!())))
    }
    let as_vec = biguint_to_u64_vec(order.clone());

    Ok((as_vec, order_len, order, rest))
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
    let mut result = serialize_fp_fixed_len(modulus_len, &x)?;
    result.extend(serialize_fp_fixed_len(modulus_len, &y)?);

    Ok(result)
}

pub(crate) fn get_base_field_params(bytes: &[u8]) -> Result<((BigUint, usize), &[u8]), ApiError> {
    let (modulus_len, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get modulus length")?;
    let modulus_len = modulus_len[0] as usize;

    let (modulus_encoding, rest) = split(rest, modulus_len, "Input is not long enough to get modulus")?;
    let modulus = BigUint::from_bytes_be(&modulus_encoding);
    if modulus.is_zero() {
        return Err(ApiError::UnexpectedZero("Modulus can not be zero".to_owned()));
    }

    Ok(((modulus, modulus_len), rest))
}

pub(crate) fn get_g1_curve_params(bytes: &[u8]) -> Result<((&[u8], usize), &[u8]), ApiError> {

    let (order_len, rest) = split(bytes, BYTES_FOR_LENGTH_ENCODING, "Input is not long enough to get group size length")?;
    let order_len = order_len[0] as usize;
    let (order_encoding, rest) = split(rest, order_len, "Input is not long enough to get main group order size")?;

    Ok(((order_encoding, order_len), rest))
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
    let (x_encoding, rest) = split(bytes, field_byte_len, "Input is not long enough to get X")?;
    let x = Fp::from_be_bytes(curve.params.params(), x_encoding, true).map_err(|_| {
        ApiError::InputError("Failed to parse X".to_owned())
    })?;
    let (y_encoding, rest) = split(rest, field_byte_len, "Input is not long enough to get Y")?;
    let y = Fp::from_be_bytes(curve.params.params(), y_encoding, true).map_err(|_| {
        ApiError::InputError("Failed to parse Y".to_owned())
    })?;
    
    let p: CurvePoint<'a, C> = CurvePoint::point_from_xy(&curve, x, y);
    
    Ok((p, rest))
}

pub(crate) fn decode_scalar_representation<
    'a
    >
    (
        bytes: &'a [u8], 
        order_byte_len: usize,
        order: &BigUint,
        order_repr: &[u64],
    ) -> Result<(Vec<u64>, &'a [u8]), ApiError>
{
    use crate::field::biguint_to_u64_vec;
    let (encoding, rest) = split(bytes, order_byte_len, "Input is not long enough to get scalar")?;
    let scalar = BigUint::from_bytes_be(&encoding);
    if &scalar >= order {
        return Err(ApiError::InputError(format!("Group order is zero, file {}, line {}", file!(), line!())));
    }
    let mut repr = biguint_to_u64_vec(scalar);
    if repr.len() < order_repr.len() {
        repr.resize(order_repr.len(), 0u64);
    }

    Ok((repr, rest))
}

