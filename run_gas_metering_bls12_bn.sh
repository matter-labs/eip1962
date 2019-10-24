#!/bin/sh
cargo clean -p eth_pairings
GAS_METERING=1 cargo test --release -- --nocapture --ignored run_bls12_bn_pseudo_curves_monte_carlo