pub(crate) mod bls12;
pub(crate) mod bn;
pub(crate) mod mnt4;
pub(crate) mod mnt6;

use crate::public_interface::{PairingApi, PublicPairingApi};
use crate::errors::ApiError;

pub(crate) fn call_pairing_engine(bytes: &[u8]) -> Result<Vec<u8>, ApiError> {
    PublicPairingApi::pair(&bytes)
}