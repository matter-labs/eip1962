#!/bin/sh
cargo clean -p eth_pairings
# GAS_METERING=1 cargo test --release -- --nocapture --ignored run_arithmetic_ops_pseudo_curves_monte_carlo
GAS_METERING=1 cargo test --release -- --nocapture --ignored run_deterministic_search_over_parameter_space_for_g1_and_g2