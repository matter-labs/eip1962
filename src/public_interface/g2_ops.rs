use crate::weierstrass::Group;
use crate::weierstrass::twist;
use crate::weierstrass::cubic_twist;
use crate::field::{SizedPrimeField, field_from_modulus};
use crate::fp::Fp;
use crate::representation::ElementRepr;

use num_bigint::BigUint;
use num_traits::{Zero};


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