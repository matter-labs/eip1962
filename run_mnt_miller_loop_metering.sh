#!/bin/sh
cargo clean -p eth_pairings
RAYON_NUM_THREADS=4 NUM_SAMPLES=3000 cargo test --release --features=gas_metering_mode -- --nocapture --ignored parallel_measure_miller_loop_pairing_costs_mnt

