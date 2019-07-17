#!/bin/sh
NUM_JOBS=12 cargo +nightly fuzz run fuzz_target_compare --jobs=12 -- -max_len=8192