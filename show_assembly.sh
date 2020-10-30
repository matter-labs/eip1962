#!/bin/sh
RUSTFLAGS="-C target-feature=+bmi2" cargo +nightly asm crate::field::U384Repr::mont_mul_assign --rust