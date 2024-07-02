#![cfg_attr(not(feature = "std"), no_std)]

#[macro_use]
#[cfg(not(feature = "std"))]
extern crate alloc as std;

extern crate uint;

mod group;
mod loop_param;

pub use self::group::MaxGroupSizeUint;
pub use self::loop_param::MaxLoopParametersUint;