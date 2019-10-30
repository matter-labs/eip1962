#!/bin/sh
cargo clean -p eth_pairings
GAS_METERING=1 cargo test --release -- --nocapture run_single_curve_arithmetic_ops