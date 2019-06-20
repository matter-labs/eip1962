use crate::weierstrass::Group;
use crate::weierstrass::twist;
use crate::weierstrass::cubic_twist;
use crate::field::{SizedPrimeField, field_from_modulus};
use crate::fp::Fp;
use crate::extension_towers::fp2;
use crate::extension_towers::fp3;
use crate::representation::ElementRepr;

use num_bigint::BigUint;
use num_traits::{Zero};

pub(crate) fn decode_fp<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>
    >
    (
        bytes: &'a [u8], 
        field_byte_len: usize,
        base_field: &'a F
    ) -> Result<(Fp<'a, FE, F>, &'a [u8]), ()>
{
    if bytes.len() < field_byte_len {
        return Err(());
    }
    let (x_encoding, rest) = bytes.split_at(field_byte_len);
    let x = Fp::from_be_bytes(base_field, x_encoding, true).map_err(|_| ())?;

    Ok((x, rest))
}

pub(crate) fn decode_fp2<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>
    >
    (
        bytes: &'a [u8], 
        field_byte_len: usize,
        extension_field: &'a fp2::Extension2<'a, FE, F>
    ) -> Result<(fp2::Fp2<'a, FE, F>, &'a [u8]), ()>
{
    if bytes.len() < field_byte_len {
        return Err(());
    }
    let (c0_encoding, rest) = bytes.split_at(field_byte_len);
    let c0 = Fp::from_be_bytes(extension_field.field, c0_encoding, true).map_err(|_| ())?;

    if rest.len() < field_byte_len {
        return Err(());
    }
    let (c1_encoding, rest) = rest.split_at(field_byte_len);
    let c1 = Fp::from_be_bytes(extension_field.field, c1_encoding, true).map_err(|_| ())?;

    let mut x = fp2::Fp2::zero(extension_field);
    x.c0 = c0;
    x.c1 = c1;

    Ok((x, rest))
}

pub(crate) fn decode_fp3<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>
    >
    (
        bytes: &'a [u8], 
        field_byte_len: usize,
        extension_field: &'a fp3::Extension3<'a, FE, F>
    ) -> Result<(fp3::Fp3<'a, FE, F>, &'a [u8]), ()>
{
    if bytes.len() < field_byte_len {
        return Err(());
    }
    let (c0_encoding, rest) = bytes.split_at(field_byte_len);
    let c0 = Fp::from_be_bytes(extension_field.field, c0_encoding, true).map_err(|_| ())?;

    if rest.len() < field_byte_len {
        return Err(());
    }
    let (c1_encoding, rest) = rest.split_at(field_byte_len);
    let c1 = Fp::from_be_bytes(extension_field.field, c1_encoding, true).map_err(|_| ())?;

    if rest.len() < field_byte_len {
        return Err(());
    }
    let (c2_encoding, rest) = rest.split_at(field_byte_len);
    let c2 = Fp::from_be_bytes(extension_field.field, c2_encoding, true).map_err(|_| ())?;

    let mut x = fp3::Fp3::zero(extension_field);
    x.c0 = c0;
    x.c1 = c1;
    x.c2 = c2;

    Ok((x, rest))
}

pub(crate) fn serialize_fp_fixed_len<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>
    >
    (
        field_byte_len: usize,
        element: &'a Fp<'a, FE, F>
    ) -> Result<Vec<u8>, ()>
{
    let mut bytes: Vec<u8> = vec![];
    let element = element.into_repr();
    element.write_be(&mut bytes).map_err(|_| ())?;
    if bytes.len() > field_byte_len {
        bytes.reverse();
        bytes.truncate(field_byte_len);
        bytes.reverse();
    } else if bytes.len() < field_byte_len {
        bytes.reverse();
        bytes.resize(field_byte_len, 0u8);
        bytes.reverse();
    }

    Ok(bytes)
}

pub(crate) fn serialize_fp2_fixed_len<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>
    >
    (
        field_byte_len: usize,
        element: &'a fp2::Fp2<'a, FE, F>
    ) -> Result<Vec<u8>, ()>
{
    let mut bytes = serialize_fp_fixed_len(field_byte_len, &element.c0)?;
    bytes.extend(serialize_fp_fixed_len(field_byte_len, &element.c1)?);

    Ok(bytes)
}

pub(crate) fn serialize_fp3_fixed_len<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>
    >
    (
        field_byte_len: usize,
        element: &'a fp3::Fp3<'a, FE, F>
    ) -> Result<Vec<u8>, ()>
{
    let mut bytes = serialize_fp_fixed_len(field_byte_len, &element.c0)?;
    bytes.extend(serialize_fp_fixed_len(field_byte_len, &element.c1)?);
    bytes.extend(serialize_fp_fixed_len(field_byte_len, &element.c2)?);

    Ok(bytes)
}