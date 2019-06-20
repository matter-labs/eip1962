use crate::weierstrass::Group;
use crate::weierstrass::twist;
use crate::weierstrass::cubic_twist;
use crate::field::{SizedPrimeField, field_from_modulus};
use crate::fp::Fp;
use crate::representation::ElementRepr;
use crate::multiexp::peppinger;

use crate::field::{
    U256Repr, 
    U320Repr,
    U384Repr,
    U448Repr,
    U512Repr,
    U576Repr,
    U640Repr,
    U704Repr,
    U768Repr,
    U832Repr,
    U896Repr,
    // U960Repr,
    // U1024Repr
};

use num_bigint::BigUint;

use super::decode_utils::parse_encodings_in_extension;

use super::decode_g2::*;

use super::decode_g1::decode_scalar_representation;
use super::decode_g1::get_base_field_params;
use super::decode_g1::get_g1_curve_params;
use super::decode_fp::{decode_fp2, decode_fp3};
use super::constants::*;

use crate::errors::ApiError;

/// Every call has common parameters (may be redundant):
/// - Lengths of modulus (in bytes)
/// - Field modulus
/// - Extension degree (2/3)
/// - Non-residue
/// - Curve A in Fpk
/// - Curve B in Fpk
/// - Length of a scalar field (curve order) (in bytes)
/// - Curve order

pub trait G2Api {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError>;
}

pub struct G2ApiImplementationFp2<FE: ElementRepr, GE: ElementRepr> {
    _marker_fe: std::marker::PhantomData<FE>,
    _marker_ge: std::marker::PhantomData<GE>
}

macro_rules! get_ab_in_fp2_extension_field {
    ($rest: expr, $field: expr, $modulus_len: expr) => {
        {
            if $rest.len() < 4*$modulus_len {
                return Err(ApiError::InputError("Input length is not enough to get A and B in Fp2".to_owned()));
            }
            let (a, rest) = decode_fp2(& $rest, $modulus_len, & $field)?;
            let (b, rest) = decode_fp2(& rest, $modulus_len, & $field)?;
            
            (a, b, rest)
        }
    }
}

macro_rules! get_ab_in_fp3_extension_field {
    ($rest: expr, $field: expr, $modulus_len: expr) => {
        {
            if $rest.len() < $modulus_len {
                return Err(());
            }
            let (a, rest) = decode_fp3(& $rest, $modulus_len, & $field)?;
            let (b, rest) = decode_fp3(& rest, $modulus_len, & $field)?;
            
            (a, b, rest)
        }
    }
}

impl<FE: ElementRepr, GE: ElementRepr> G2Api for G2ApiImplementationFp2<FE, GE> {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, _, rest) = create_base_field_with_modulus!(bytes, FE);
        let (extension_2, rest) = create_fp2_extension(rest, modulus_len, &field)?;
        let (a, b, rest) = get_ab_in_fp2_extension_field!(rest, extension_2, modulus_len);
        let (group, order_len, rest) = create_group!(rest, GE);

        let curve = twist::WeierstrassCurveTwist::new(&group, &extension_2, a, b);

        let (mut p_0, rest) = decode_g2_point_from_xy_in_fp2(rest, modulus_len, &curve)?;
        let (p_1, _rest) = decode_g2_point_from_xy_in_fp2(rest, modulus_len, &curve)?;

        p_0.add_assign(&p_1);

        serialize_g2_point_in_fp2(modulus_len, &p_0)   
    }

    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (field, modulus_len, _, rest) = create_base_field_with_modulus!(bytes, FE);
        let (extension_2, rest) = create_fp2_extension(rest, modulus_len, &field)?;
        let (a, b, rest) = get_ab_in_fp2_extension_field!(rest, extension_2, modulus_len);
        let (group, order_len, rest) = create_group!(rest, GE);

        let curve = twist::WeierstrassCurveTwist::new(&group, &extension_2, a, b);

        let (p_0, rest) = decode_g2_point_from_xy_in_fp2(rest, modulus_len, &curve)?;
        let (scalar, _rest) = decode_scalar_representation(rest, order_len, &group)?;

        let p = p_0.mul(&scalar);

        serialize_g2_point_in_fp2(modulus_len, &p)   
    }

    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        unimplemented!();
        // let (field, modulus_len, _, rest) = create_base_field_with_modulus!(bytes, FE);
        // let (extension_2, rest) = create_fp2_extension(&rest, modulus_len, &field)?;
        // let (a, b, rest) = get_ab_in_fp2_extension_field!(rest, extension_2, modulus_len);
        // let (group, order_len, rest) = create_group!(rest, GE);

        // let curve = twist::WeierstrassCurveTwist::new(&group, &extension_2, a, b);

        // let expected_pair_len = 4*modulus_len + order_len;
        // if rest.len() % expected_pair_len != 0 {
        //     return Err(());
        // }

        // let expected_pairs = rest.len() / expected_pair_len;
        // if expected_pairs == 0 {
        //     return Err(());
        // }

        // let mut global_rest = rest;
        // let mut pairs = vec![];

        // for _ in 0..expected_pairs {
        //     let (p, local_rest) = decode_g2_point_from_xy_in_fp2(global_rest, modulus_len, &curve)?;
        //     let (scalar, local_rest) = decode_scalar_representation(local_rest, order_len, &group)?;
        //     pairs.push((p, scalar));
        //     global_rest = local_rest;
        // }

        // let result = peppinger(pairs);

        // serialize_g2_point_in_fp2(modulus_len, &result)   
    }
}

pub struct PublicG2Api;

impl G2Api for PublicG2Api {
    fn add_points(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (modulus, _, extension_degree, _, _, order, _, _) = parse_encodings_in_extension(&bytes)?;
        let modulus_limbs = (modulus.bits() / 64) + 1;
        let order_limbs = (order.bits() / 64) + 1;

        let result: Result<Vec<u8>, ApiError> = match extension_degree {
            EXTENSION_DEGREE_2 => {
                let result: Result<Vec<u8>, ApiError> = match (modulus_limbs, order_limbs) {
                    (4, 4) => {
                        G2ApiImplementationFp2::<U256Repr, U256Repr>::add_points(&bytes)
                    },
                    (5, 4) => {
                        G2ApiImplementationFp2::<U320Repr, U256Repr>::add_points(&bytes)
                    },
                    (5, 5) => {
                        G2ApiImplementationFp2::<U320Repr, U320Repr>::add_points(&bytes)
                    },
                    (field_limbs, group_limbs) => {
                        unimplemented!("unimplemented for {} modulus and {} group limbs", field_limbs, group_limbs);
                    }
                };

                result
            },
            _ => {
                unimplemented!("Extension degree other than 2 is not yet implemented");
            }
        };

        result
    }

    fn mul_point(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (modulus, _, extension_degree, _, _, order, _, _) = parse_encodings_in_extension(&bytes)?;
        let modulus_limbs = (modulus.bits() / 64) + 1;
        let order_limbs = (order.bits() / 64) + 1;

        let result: Result<Vec<u8>, ApiError> = match extension_degree {
            EXTENSION_DEGREE_2 => {
                let result: Result<Vec<u8>, ApiError> = match (modulus_limbs, order_limbs) {
                    (4, 4) => {
                        G2ApiImplementationFp2::<U256Repr, U256Repr>::mul_point(&bytes)
                    },
                    (5, 4) => {
                        G2ApiImplementationFp2::<U320Repr, U256Repr>::mul_point(&bytes)
                    },
                    (5, 5) => {
                        G2ApiImplementationFp2::<U320Repr, U320Repr>::mul_point(&bytes)
                    },
                    (field_limbs, group_limbs) => {
                        unimplemented!("unimplemented for {} modulus and {} group limbs", field_limbs, group_limbs);
                    }
                };

                result
            },
            _ => {
                unimplemented!("Extension degree other than 2 is not yet implemented");
            }
        };

        result
    }

    fn multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
        let (modulus, _, extension_degree, _, _, order, _, _) = parse_encodings_in_extension(&bytes)?;
        let modulus_limbs = (modulus.bits() / 64) + 1;
        let order_limbs = (order.bits() / 64) + 1;

        let result: Result<Vec<u8>, ApiError> = match extension_degree {
            EXTENSION_DEGREE_2 => {
                let result: Result<Vec<u8>, ApiError> = match (modulus_limbs, order_limbs) {
                    (4, 4) => {
                        G2ApiImplementationFp2::<U256Repr, U256Repr>::add_points(&bytes)
                    },
                    (5, 4) => {
                        G2ApiImplementationFp2::<U320Repr, U256Repr>::add_points(&bytes)
                    },
                    (5, 5) => {
                        G2ApiImplementationFp2::<U320Repr, U320Repr>::add_points(&bytes)
                    },
                    (field_limbs, group_limbs) => {
                        unimplemented!("unimplemented for {} modulus and {} group limbs", field_limbs, group_limbs);
                    }
                };

                result
            },
            _ => {
                unimplemented!("Extension degree other than 2 is not yet implemented");
            }
        };

        result
    }
}