// Copyright 2015-2017 Parity Technologies
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

//! Efficient large, fixed-size big integers and hashes.

#![cfg_attr(not(feature = "std"), no_std)]

#[doc(hidden)]
pub use byteorder;

// Re-export libcore using an alias so that the macros can work without
// requiring `extern crate core` downstream.
#[doc(hidden)]
pub use core as core_;

#[doc(hidden)]
pub use rustc_hex;

#[cfg(feature = "quickcheck")]
#[doc(hidden)]
pub use qc;

#[cfg(feature = "quickcheck")]
#[doc(hidden)]
pub use rand;

#[doc(hidden)]
pub use static_assertions;

#[cfg(feature = "unroll")]
pub use crunchy::unroll;

#[macro_use]
#[rustfmt::skip]
mod uint;
pub use crate::uint::*;

pub use self::arith_impl::*;

mod arith_impl {

    /// Calculate a - b - borrow, returning the result and modifying
    /// the borrow value.
    #[inline(always)]
    pub fn sbb(a: u64, b: u64, borrow: &mut u64) -> u64 {

        let tmp = (1u128 << 64).wrapping_add(u128::from(a)).wrapping_sub(u128::from(b)).wrapping_sub(u128::from(*borrow));

        *borrow = if tmp >> 64 == 0 { 1 } else { 0 };

        tmp as u64
    }

    /// Calculate a + b + carry, returning the sum and modifying the
    /// carry value.
    #[inline(always)]
    pub fn adc(a: u64, b: u64, carry: &mut u64) -> u64 {

        let tmp = u128::from(a).wrapping_add(u128::from(b)).wrapping_add(u128::from(*carry));

        *carry = (tmp >> 64) as u64;

        tmp as u64
    }

    /// Calculate a + carry, returning the sum and modifying the
    /// carry value.
    #[inline(always)]
    pub fn add_carry(a: u64, carry: &mut u64) -> u64 {

        let tmp = u128::from(a).wrapping_add(u128::from(*carry));

        *carry = (tmp >> 64) as u64;

        tmp as u64
    }

    /// Calculate a + (b * c) + carry, returning the least significant digit
    /// and setting carry to the most significant digit.
    #[inline(always)]
    pub fn mac_with_carry(a: u64, b: u64, c: u64, carry: &mut u64) -> u64 {

        let tmp = (u128::from(a)).wrapping_add(u128::from(b).wrapping_mul(u128::from(c))).wrapping_add(u128::from(*carry));

        *carry = (tmp >> 64) as u64;

        tmp as u64
    }
}