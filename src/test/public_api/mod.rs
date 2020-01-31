use crate::public_interface::{OperationType, perform_operation, ApiError};

use crate::weierstrass::{Group, CurveOverFpParameters};
use crate::weierstrass::curve::{CurvePoint, WeierstrassCurve};
use crate::representation::ElementRepr;
use crate::field::*;

use crate::public_interface::decode_g1::*;
use crate::public_interface::decode_utils::*;
use crate::public_interface::decode_fp::*;

#[macro_use]
use crate::{expand_for_modulus_limbs};

fn call_public_api_on_test_vector(data: &[u8]) -> Result<Vec<u8>, ApiError>{
    if data.len() == 0 {
        return Err(ApiError::InputError("input is zero length".to_owned()));
    }
    let op = OperationType::from_u8(data[0]).ok_or(ApiError::MissingValue)?;

    perform_operation(op, &data[0..])
}

#[test]
fn test_api_mul_g1() {
    use hex;

    let hex_string = "02820d8080001026e1318f230000000000000080be6dc885e544bee65620747e191023316695af6989f85b230000000000044eccc6886286fbaee7561d483d50d89e9e9e95af2989f85b230000000000044eccc688628631";

    let data = hex::decode(hex_string).unwrap();

    let result = call_public_api_on_test_vector(&data);

    if let Ok(result) = result {
        println!("Result = {}", hex::encode(&result));
    } else {
        println!("Error = {}", result.err().unwrap());
    }
}

#[test]
fn test_api_different_mul_g1() {
    use hex;

    let hex_string = "02820d8080001026e1318f230000000000000080be6dc885e544bee65620747e191023316695af6989f85b230000000000044eccc6886286fbaee7561d483d50d89e9e9e95af2989f85b230000000000044eccc688628631";

    let data = hex::decode(hex_string).unwrap();

    let result = test_different_g1_multiplications(&data);

    if let Ok(_) = result {
        println!("Result = ok");
    } else {
        println!("Error = {}", result.err().unwrap());
    }
}

pub(crate) fn test_different_g1_multiplications(bytes: &[u8]) -> Result<(), ApiError>  {
    let (_, modulus, _) = parse_modulus_and_length(&bytes)?;
    let modulus_limbs = num_limbs_for_modulus(&modulus)?;

    let result: Result<(), ApiError> = expand_for_modulus_limbs!(modulus_limbs, Tester, bytes, parse_and_compare_muls); 

    result
}

struct Tester<FE: ElementRepr> {
    _marker: std::marker::PhantomData<FE>
}

impl<FE: ElementRepr> Tester<FE> {
    fn parse_and_compare_muls(bytes: &[u8]) -> Result<(), ApiError> {
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

        let (x_double_and_add, y_double_and_add) = p_0.mul_impl(scalar.as_ref()).into_xy();
        let (x_wnaf_mul, y_wnaf_mul) = p_0.wnaf_mul_impl(scalar.as_ref()).into_xy();

        if x_double_and_add != x_wnaf_mul || y_double_and_add != y_wnaf_mul {
            return Err(ApiError::InputError(
                format!("DoubleAndAdd x = {}, y = {}, Wnaf x = {}, y = {}", x_double_and_add, y_double_and_add, x_wnaf_mul, y_wnaf_mul)
            ));
        }

        return Ok(())
    }
}