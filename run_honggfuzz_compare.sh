#!/bin/sh
cd honggfuzz
NUM_JOBS=12 cargo hfuzz run fuzz_target_compare