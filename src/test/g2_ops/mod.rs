pub(crate) mod bls12;
pub(crate) mod bn;
pub(crate) mod mnt4;
pub(crate) mod mnt6;

use crate::public_interface::{G2Api, PublicG2Api};
use crate::errors::ApiError;

pub(crate) fn call_g2_engine_add(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
    PublicG2Api::add_points(&bytes)
}

pub(crate) fn call_g2_engine_mul(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
    PublicG2Api::mul_point(&bytes)
}

pub(crate) fn call_g2_engine_multiexp(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
    PublicG2Api::multiexp(&bytes)
}