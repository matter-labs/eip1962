#!/bin/sh
cargo clean -p eth_pairings
GAS_METERING=1 cargo test --no-run --release
# GAS_METERING=1 flamegraph -o my_flamegraph.svg ./target/debug/eth_pairings-4c90044880fe6f08 --nocapture run_single_curve_arithmetic_ops
# GAS_METERING=1 sudo flamegraph -o my_flamegraph.svg ./target/debug/eth_pairings-4c90044880fe6f08 --nocapture run_flamer_one_off_bls12
GAS_METERING=1 sudo flamegraph -o flamegraph_g2_const.svg ./target/release/eth_pairings-90094f439ad4e815 --nocapture flame_one_off_field_and_curve_constuction_g2_ext3
# GAS_METERING=1 sudo flamegraph -o my_flamegraph.svg ./target/release/eth_pairings-b6574d4c7bd68eef --nocapture run_flamer_one_off_bls12

