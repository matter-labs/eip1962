#!/bin/sh
cargo +nightly fuzz run fuzz_target_api --jobs=12 -- -max_len=8192