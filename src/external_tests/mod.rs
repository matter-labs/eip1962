use crate::weierstrass::{Group, CurveOverFpParameters};
use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
use crate::representation::ElementRepr;
use crate::multiexp::peppinger;
use crate::field::*;
use super::constants::*;

use super::decode_g1::*;
use super::decode_utils::*;
use super::decode_fp::*;

use crate::errors::ApiError;

pub fn test_inversion(input: &[u8]) -> Result<Vec<u8>, ApiError> {

}

pub fn test_different_g1_multiplications(bytes: &[u8]) -> Result<(), ApiError>  {
    let (_, modulus, _) = parse_modulus_and_length(&bytes)?;
    let modulus_limbs = num_limbs_for_modulus(&modulus)?;

    let result: Result<(), ApiError> = expand_for_modulus_limbs!(modulus_limbs, Tester, bytes, parse_and_compare_muls); 

    result
}

struct Tester<FE: ElementRepr> {
    _marker:: std::marker::phantom_data<FE>
}

impl<FE: ElementRepr> Tester<FE> {
    fn parse_and_compare_muls(input: &[u8]) -> Result<(), ApiError> {
        let (field, modulus_len, _, rest) = parse_base_field_from_encoding::<FE>(&bytes)?;
        let (a, b, rest) = parse_ab_in_base_field_from_encoding(&rest, modulus_len, &field)?;
        let (order_len, order, rest) = parse_group_order_from_encoding(rest)?;

        let fp_params = CurveOverFpParameters::new(&field);

        let curve = WeierstrassCurve::new(&order.as_ref(), a, b, &fp_params).map_err(|_| {
            ApiError::InputError("Curve shape is not supported".to_owned())
        })?;

        let (p_0, rest) = decode_g1_point_from_xy(rest, modulus_len, &curve)?;
        let (scalar, rest) = decode_scalar_representation(rest, order_len, &order)?;

        if rest.len() != 0 {
            return Err(ApiError::InputError("Input contains garbage at the end".to_owned()));
        }

        if !p_0.is_on_curve() {
            if !crate::features::in_fuzzing_or_gas_metering() {
                return Err(ApiError::InputError(format!("Point is not on curve, file {}, line {}", file!(), line!())));
            }
        }

        let (x_double_and_add, y_double_and_add) = p_0.mul_impl(&scalar).into_xy();
        let (x_wnaf_mul, y_wnaf_mul) = p_0.wnaf_mul_impl(&scalar).into_xy();

        if x_double_and_add != x_wnaf_mul || y_double_and_add != y_wnaf_mul {
            return Err(ApiError::MissingValue);
        }

        return Ok(())
    }
}