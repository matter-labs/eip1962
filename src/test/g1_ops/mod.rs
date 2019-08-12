pub(crate) mod bls12;
pub(crate) mod bn;

use crate::public_interface::{G1Api, PublicG1Api};
use crate::errors::ApiError;

pub(crate) fn call_g1_engine_add(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
    PublicG1Api::add_points(&bytes)
}

pub(crate) fn call_g1_engine_mul(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
    PublicG1Api::mul_point(&bytes)
}

pub(crate) fn call_g1_engine_multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
    PublicG1Api::multiexp(&bytes)
}