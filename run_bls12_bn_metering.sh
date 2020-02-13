#!/bin/sh
cargo clean -p eth_pairings
RAYON_NUM_THREADS=4 GAS_METERING=1 cargo test --release -- --nocapture --ignored parallel_measure_bls12_bn_pairing_costs

