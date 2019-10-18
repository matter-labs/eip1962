#!/bin/sh
cargo clean -p eth_pairings
GAS_METERING=1 cargo test --release -- --nocapture run_pseudo_curves_monte_carlo