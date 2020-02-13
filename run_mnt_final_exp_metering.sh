#!/bin/sh
cargo clean -p eth_pairings
RAYON_NUM_THREADS=4 GAS_METERING=1 cargo test --release --features=gas_metering_mode -- --nocapture --ignored parallel_measure_final_exp_pairing_costs_mnt

