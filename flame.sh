#!/bin/sh
cargo clean -p eth_pairings
GAS_METERING=1 cargo test --no-run
GAS_METERING=1 flamegraph -o my_flamegraph.svg ./target/debug/eth_pairings-4c90044880fe6f08 --nocapture run_single_curve_arithmetic_ops