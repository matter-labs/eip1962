#!/bin/sh
cargo clean -p eth_pairings
# GAS_METERING=1 cargo +nightly test --release -- --ignored --nocapture run_single_one_off_test
GAS_METERING=1 cargo test --release -- --ignored --nocapture run_single_one_off_test
# GAS_METERING=1 cargo test --release -- --ignored --nocapture run_single_curve_and_fields_construction
# GAS_METERING=1 cargo test --release -- --nocapture run_single_curve_arithmetic_ops