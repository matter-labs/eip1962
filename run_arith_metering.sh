#!/bin/sh
cargo clean -p eth_pairings
RAYON_NUM_THREADS=4 NUM_SAMPLES=200 cargo test --release --features=gas_metering_mode -- --nocapture --ignored run_deterministic_parallel_search_no_filtering
