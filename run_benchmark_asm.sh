#!/bin/sh
# cargo +nightly bench --features=benchmarks bench_bls12_381_engine_
# cargo +nightly bench --features=benchmarks bench_bls12_377_engine_
RUSTFLAGS="-C target-feature=+bmi2" cargo +nightly bench --features=benchmarks,asm bench_eip_2537_