#![allow(dead_code)]

extern crate eth_pairings;
extern crate eth_pairings_cpp;
extern crate hex;

mod run_on_csv;
mod run_on_fuzzer_inputs;

fn run(data: &[u8]) {
    println!("Input = {}", hex::encode(&data));
    let native = eth_pairings::public_interface::API::run(&data);
    let cpp = eth_pairings_cpp::run(&data);
    if native.is_err() {
        if !cpp.is_err() {
            let n = native.err();
            let c = cpp.unwrap();
            println!("Native result returned error {:?}, while C++ returned {}", n, hex::encode(&c));
            panic!("Native result returned error {:?}, while C++ returned {}", n, hex::encode(&c));
        }
    } else {
        let n = native.expect("result");
        if cpp.is_err() {
            let c = cpp.err();
            println!("Native result = {}, while C++ returned error {:?}", hex::encode(&n), c);
            panic!("Native result = {}, while C++ returned error {:?}", hex::encode(&n), c);
        }
        let c = cpp.expect("cpp result");
        if n != c {
            println!("Native result = {}, C++ result = {}", hex::encode(&n), hex::encode(&c));
            panic!("Native result = {}, C++ result = {}", hex::encode(&n), hex::encode(&c));
        } else {
            println!("Native and C++ results coincide on {}", hex::encode(&n));
        }
    }
}

#[test]
fn cross_check_on_input() {
    let filename = "slow-unit-f3a20e0db984a4a6759a80b01af14b7b6ae30367";
    use std::time::Instant;
    use std::io::Read;
    use std::fs::File;
    let mut file = File::open(&format!("../artifacts/fuzz_target_api/{}", filename)).expect("must open");
    let mut input_data = vec![];
    file.read_to_end(&mut input_data).expect("must read");
    assert!(input_data.len() != 0);
    let now = Instant::now();
    self::run(&input_data[..]);
    let elapsed = now.elapsed().as_micros();
    println!("Api call taken in {} micros", elapsed);
}

#[test]
fn cross_check_on_hongg_input() {
    let filename = "SIGABRT.EXC_CRASH.PC.00007fff7543a2c6.STACK.00000003402ccb03.ADDR.0000000000000000.fuzz";
    use std::time::Instant;
    use std::io::Read;
    use std::fs::File;
    let mut file = File::open(&format!("../honggfuzz/hfuzz_workspace/fuzz_target_compare/{}", filename)).expect("must open");
    let mut input_data = vec![];
    file.read_to_end(&mut input_data).expect("must read");
    assert!(input_data.len() != 0);
    let now = Instant::now();
    self::run(&input_data[..]);
    let elapsed = now.elapsed().as_micros();
    println!("Api call taken in {} micros", elapsed);
}