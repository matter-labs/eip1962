#![no_main]
#[macro_use] extern crate libfuzzer_sys;
extern crate eth_pairings;
extern crate eth_pairings_cpp;
extern crate hex;

fuzz_target!(|data: &[u8]| {
    let native = eth_pairings::public_interface::API::run(&data);
    let cpp = eth_pairings_cpp::run(&data);
    if native.is_err() {
        if !cpp.is_err() {
            panic!("Native result returned error, while C++ returned {}", hex::encode(&cpp.unwrap()))
        }
        assert!(cpp.is_err());
    } else {
        let n = native.expect("result");
        let c = cpp.expect("cpp result");
        if n != c {
            panic!("Native result = {}, C++ result = {}", hex::encode(&n), hex::encode(&c));
        }
    }
});
