#!/bin/sh
cargo clean -p eth_pairings
RAYON_NUM_THREADS=4 NUM_SAMPLES=1000 cargo test --release --features=gas_metering_mode -- --nocapture --ignored parallel_measure_one_off_pairing_costs

