#![cfg_attr(not(feature = "std"), no_std)]

#[macro_use]
#[cfg(not(feature = "std"))]
extern crate alloc as std;

extern crate uint;

mod field;
mod field_construction;

pub use self::field::MaxFieldUint;
pub use self::field_construction::MaxFieldSquaredUint;