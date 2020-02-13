#!/bin/sh
cargo clean -p eth_pairings
RAYON_NUM_THREADS=4 NUM_SAMPLES=5000 cargo test --release --features=gas_metering_mode -- --nocapture --ignored parallel_measure_bls12_bn_pairing_costs

