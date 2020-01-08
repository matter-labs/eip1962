use crate::field::{SizedPrimeField};
use crate::fp::Fp;
use crate::extension_towers::*;
use crate::extension_towers::fp2;
use crate::extension_towers::fp3;
use crate::representation::{ElementRepr};
use crate::weierstrass::curve::{WeierstrassCurve, CurvePoint};
use crate::traits::FieldElement;
use crate::weierstrass::CurveParameters;
use crate::constants::MaxFieldUint;

use super::decode_fp::*;
use super::constants::*;
use super::decode_utils::split;

use crate::errors::ApiError;

pub(crate) fn create_fp2_extension<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>,
    >
(
    bytes: &'a [u8], 
    modulus: &MaxFieldUint,
    field_byte_len: usize,
    base_field: &'a F,
    need_frobenius: bool
) -> Result<(fp2::Extension2<'a, FE, F>, &'a [u8]), ApiError>
{
    let (extension_degree, rest) = split(bytes, EXTENSION_DEGREE_ENCODING_LENGTH, "Input is not long enough to get extension degree")?;
    if extension_degree[0] != EXTENSION_DEGREE_2 {
        return Err(ApiError::UnknownParameter("Extension degree expected to be 2".to_owned()));
    }

    let (fp_non_residue, rest): (Fp<'a, FE, F>, _) = decode_fp(&rest, field_byte_len, base_field)?;
    if fp_non_residue.is_zero() {
        return Err(ApiError::UnexpectedZero("Fp2 non-residue can not be zero".to_owned()));
    }

    {
        let not_a_square = is_non_nth_root(&fp_non_residue, modulus, 2);
        if !not_a_square {
            if !crate::features::in_fuzzing_or_gas_metering() {
                return Err(ApiError::InputError(format!("Non-residue for Fp2 is actually a residue, file {}, line {}", file!(), line!())));
            }
        }
    }

    let mut extension_2 = fp2::Extension2::new(fp_non_residue);
    if need_frobenius {
        extension_2.calculate_frobenius_coeffs(modulus).map_err(|_| {
            ApiError::UnknownParameter("Failed to calculate Frobenius coeffs for Fp2".to_owned())
        })?;
    }
    
    Ok((extension_2, rest))
}

pub(crate) fn create_fp3_extension<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>,
    >
(
    bytes: &'a [u8], 
    modulus: &MaxFieldUint,
    field_byte_len: usize,
    base_field: &'a F,
    need_frobenius: bool
) -> Result<(fp3::Extension3<'a, FE, F>, &'a [u8]), ApiError>
{
    let (extension_degree, rest) = split(bytes, EXTENSION_DEGREE_ENCODING_LENGTH, "Input is not long enough to get extension degree")?;
    if extension_degree[0] != EXTENSION_DEGREE_3 {
        return Err(ApiError::UnknownParameter("Extension degree expected to be 3".to_owned()));
    }

    let (fp_non_residue, rest): (Fp<'a, FE, F>, _) = decode_fp(&rest, field_byte_len, base_field)?;
    if fp_non_residue.is_zero() {
        return Err(ApiError::UnexpectedZero("Fp3 non-residue can not be zero".to_owned()));
    }

    {
        let not_a_cube = is_non_nth_root(&fp_non_residue, modulus, 3);
        if !not_a_cube {
            if !crate::features::in_fuzzing_or_gas_metering() {
                return Err(ApiError::InputError(format!("Non-residue for Fp3 is actually a residue, file {}, line {}", file!(), line!())));
            }
        }
    }

    let mut extension_3 = fp3::Extension3::new(fp_non_residue);
    if need_frobenius {
        extension_3.calculate_frobenius_coeffs(modulus).map_err(|_| {
            ApiError::UnknownParameter("Failed to calculate Frobenius coeffs for Fp3".to_owned())
        })?;
    }
    
    Ok((extension_3, rest))
}

pub(crate) fn decode_g2_point_from_xy_in_fp2<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE> + 'a,
    C: CurveParameters<BaseFieldElement = fp2::Fp2<'a, FE, F>>
    >
    (
        bytes: &'a [u8], 
        field_byte_len: usize,
        curve: &'a WeierstrassCurve<'a, C>
    ) -> Result<(CurvePoint<'a, C>, &'a [u8]), ApiError>
{
    let (x, rest) = decode_fp2(&bytes, field_byte_len, curve.params.params())?;
    let (y, rest) = decode_fp2(&rest, field_byte_len, curve.params.params())?;
    
    let p: CurvePoint<'a, C> = CurvePoint::point_from_xy(&curve, x, y);
    
    Ok((p, rest))
}

pub(crate) fn decode_g2_point_from_xy_in_fp3<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE> + 'a,
    C: CurveParameters<BaseFieldElement = fp3::Fp3<'a, FE, F>>
    >
    (
        bytes: &'a [u8], 
        field_byte_len: usize,
        curve: &'a WeierstrassCurve<'a, C>
    ) -> Result<(CurvePoint<'a, C>, &'a [u8]), ApiError>
{
    let (x, rest) = decode_fp3(&bytes, field_byte_len, curve.params.params())?;
    let (y, rest) = decode_fp3(&rest, field_byte_len, curve.params.params())?;
    
    let p: CurvePoint<'a, C> = CurvePoint::point_from_xy(&curve, x, y);
    
    Ok((p, rest))
}

pub(crate) fn serialize_g2_point_in_fp2<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE> + 'a,
    C: CurveParameters<BaseFieldElement = fp2::Fp2<'a, FE, F>>
    >
    (
        modulus_len: usize,
        point: &CurvePoint<'a, C>
    ) -> Result<Vec<u8>, ApiError>
{
    let (x, y) = point.into_xy();
    let mut result = Vec::with_capacity(4*modulus_len);
    result.extend(serialize_fp2_fixed_len(modulus_len, &x)?);
    result.extend(serialize_fp2_fixed_len(modulus_len, &y)?);
    
    Ok(result)
}

pub(crate) fn serialize_g2_point_in_fp3<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE> + 'a,
    C: CurveParameters<BaseFieldElement = fp3::Fp3<'a, FE, F>>
    >
    (
        modulus_len: usize,
        point: &CurvePoint<'a, C>
    ) -> Result<Vec<u8>, ApiError>
{
    let (x, y) = point.into_xy();
    let mut result = Vec::with_capacity(6*modulus_len);
    result.extend(serialize_fp3_fixed_len(modulus_len, &x)?);
    result.extend(serialize_fp3_fixed_len(modulus_len, &y)?);
    
    Ok(result)
}

pub(crate) fn parse_ab_in_fp2_from_encoding<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>
    >(
        encoding: &'a [u8], 
        modulus_len: usize,
        field: &'a fp2::Extension2<'a, FE, F>
    ) -> Result<(fp2::Fp2<'a, FE, F>, fp2::Fp2<'a, FE, F>, &'a [u8]), ApiError>
{
    let (a, rest) = decode_fp2(&encoding, modulus_len, field)?;
    let (b, rest) = decode_fp2(&rest, modulus_len, field)?;

    Ok((a, b, rest))
}

pub(crate) fn parse_ab_in_fp3_from_encoding<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>
    >(
        encoding: &'a [u8], 
        modulus_len: usize,
        field: &'a fp3::Extension3<'a, FE, F>
    ) -> Result<(fp3::Fp3<'a, FE, F>, fp3::Fp3<'a, FE, F>, &'a [u8]), ApiError>
{
    let (a, rest) = decode_fp3(&encoding, modulus_len, field)?;
    let (b, rest) = decode_fp3(&rest, modulus_len, field)?;

    Ok((a, b, rest))
}