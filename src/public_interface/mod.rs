#[macro_use]
mod decode_g1;

mod decode_g2;
mod decode_fp;
mod decode_utils;

mod g1_ops;
mod pairing_ops;

pub mod constants;

pub use pairing_ops::{PairingApi, PublicPairingApi};