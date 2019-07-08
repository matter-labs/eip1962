#![no_main]
#[macro_use] extern crate libfuzzer_sys;
extern crate eth_pairings;

fuzz_target!(|data: &[u8]| {
    let _ = eth_pairings::gas_meter::GasMeter::meter(&data);
});
