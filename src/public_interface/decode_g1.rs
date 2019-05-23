use crate::weierstrass::Group;
use crate::weierstrass::curve::{WeierstrassCurve, CurvePoint};
use crate::field::{SizedPrimeField, field_from_modulus};
use crate::fp::Fp;
use crate::representation::ElementRepr;

use super::constants::*;

use num_bigint::BigUint;
use num_traits::{Zero};

macro_rules! create_base_field {
    ($bytes:expr, $repr:tt) => {
        {
            let ((modulus, modulus_len), rest) = get_base_field_params($bytes)?;
            let field = field_from_modulus::<$repr>(modulus)?;
            if rest.len() < modulus_len {
                return Err(());
            }

            (field, modulus_len, rest)
        }
    }
}

macro_rules! create_base_field_with_modulus {
    ($bytes:expr, $repr:tt) => {
        {
            let ((modulus, modulus_len), rest) = get_base_field_params($bytes)?;
            let field = field_from_modulus::<$repr>(modulus.clone())?;
            if rest.len() < modulus_len {
                return Err(());
            }

            (field, modulus_len, modulus, rest)
        }
    }
}

macro_rules! get_ab_in_base_field {
    ($rest:expr, $field:expr, $modulus_len: expr) => {
        {
            if $rest.len() < $modulus_len {
                return Err(());
            }
            let (a_encoding, rest) = $rest.split_at($modulus_len);
            let a = Fp::from_be_bytes(&$field, a_encoding, true).map_err(|_| ())?;
            if rest.len() < $modulus_len {
                return Err(());
            }
            let (b_encoding, rest) = rest.split_at($modulus_len);
            let b = Fp::from_be_bytes(&$field, b_encoding, true).map_err(|_| ())?;
            
            (a, b, rest)
        }
    }
}

macro_rules! create_group {
    ($bytes:expr, $repr:tt) => {
        {
            let ((order, order_len), rest) = get_g1_curve_params($bytes)?;
            let order = BigUint::from_bytes_be(&order);
            let group = field_from_modulus::<$repr>(order)?;

            (group, order_len, rest)
        }
    }
}

pub(crate) fn serialize_g1_point<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>,
    GE: ElementRepr,
    G: SizedPrimeField<Repr = GE>
    >
    (
        modulus_len: usize,
        point: &CurvePoint<'a, FE, F, GE, G>
    ) -> Result<Vec<u8>, ()>
{
    let (x, y) = point.into_xy();
    let x = x.into_repr();
    let y = y.into_repr();

    let mut x_bytes: Vec<u8> = vec![];
    x.write_be(&mut x_bytes).map_err(|_| ())?;
    if x_bytes.len() > modulus_len {
        x_bytes.reverse();
        x_bytes.truncate(modulus_len);
        x_bytes.reverse();
    } else if x_bytes.len() < modulus_len {
        x_bytes.reverse();
        x_bytes.resize(modulus_len, 0u8);
        x_bytes.reverse();
    }

    let mut y_bytes: Vec<u8> = vec![];
    y.write_be(&mut y_bytes).map_err(|_| ())?;
    if y_bytes.len() > modulus_len {
        y_bytes.reverse();
        y_bytes.truncate(modulus_len);
        y_bytes.reverse();
    } else if y_bytes.len() < modulus_len {
        y_bytes.reverse();
        y_bytes.resize(modulus_len, 0u8);
        y_bytes.reverse();
    }

    let mut result = vec![];
    result.append(&mut x_bytes);
    result.append(&mut y_bytes);

    Ok(result)
}

pub(crate) fn get_base_field_params(bytes: &[u8]) -> Result<((BigUint, usize), &[u8]), ()> {
    if bytes.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(());
    }
    let (modulus_len, rest) = bytes.split_at(BYTES_FOR_LENGTH_ENCODING);
    let modulus_len = modulus_len[0] as usize;

    if rest.len() < modulus_len {
        return Err(());
    }
    let (modulus_encoding, rest) = rest.split_at(modulus_len);
    let modulus = BigUint::from_bytes_be(&modulus_encoding);
    if modulus.is_zero() {
        return Err(());
    }

    Ok(((modulus, modulus_len), rest))
}

pub(crate) fn get_g1_curve_params(bytes: &[u8]) -> Result<((&[u8], usize), &[u8]), ()> {
    if bytes.len() < BYTES_FOR_LENGTH_ENCODING {
        return Err(());
    }

    let (order_len, rest) = bytes.split_at(BYTES_FOR_LENGTH_ENCODING);
    let order_len = order_len[0] as usize;
    if rest.len() < order_len {
        return Err(());
    }
    let (order_encoding, rest) = rest.split_at(order_len);

    Ok(((order_encoding, order_len), rest))
}

pub(crate) fn decode_g1_point_from_xy<
    'a,
    FE: ElementRepr,
    F: SizedPrimeField<Repr = FE>,
    GE: ElementRepr,
    G: SizedPrimeField<Repr = GE>
    >
    (
        bytes: &'a [u8], 
        field_byte_len: usize,
        curve: &'a WeierstrassCurve<'a, FE, F, GE, G>
    ) -> Result<(CurvePoint<'a, FE, F, GE, G>, &'a [u8]), ()>
{
    if bytes.len() < field_byte_len {
        return Err(());
    }
    let (x_encoding, rest) = bytes.split_at(field_byte_len);
    let x = Fp::from_be_bytes(curve.base_field, x_encoding, true).map_err(|_| ())?;
    if rest.len() < field_byte_len {
        return Err(());
    }
    let (y_encoding, rest) = rest.split_at(field_byte_len);
    let y = Fp::from_be_bytes(curve.base_field, y_encoding, true).map_err(|_| ())?;
    
    let p: CurvePoint<'a, FE, F, GE, G> = CurvePoint::point_from_xy(&curve, x, y);
    
    Ok((p, rest))
}

pub(crate) fn decode_scalar_representation<
    'a,
    GE: ElementRepr,
    G: SizedPrimeField<Repr = GE>
    >
    (
        bytes: &'a [u8], 
        order_byte_len: usize,
        group: &G,
    ) -> Result<(GE, &'a [u8]), ()>
{
    if bytes.len() < order_byte_len {
        return Err(());
    }
    let (encoding, rest) = bytes.split_at(order_byte_len);
    let mut repr = GE::default();
    if encoding.len() >= repr.as_ref().len() * 8 {
        repr.read_be(encoding).map_err(|_| ())?;
    } else {
        let mut padded = vec![0u8; repr.as_ref().len() * 8 - encoding.len()];
        padded.extend_from_slice(encoding);
        repr.read_be(&padded[..]).map_err(|_| ())?;
    }

    Ok((repr, rest))
}

