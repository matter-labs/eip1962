#![no_main]
#[macro_use] extern crate libfuzzer_sys;
extern crate eth_pairings;
extern crate eth_pairings_cpp;
extern crate hex;

fuzz_target!(|data: &[u8]| {
    let native = eth_pairings::public_interface::API::run(&data);
    let cpp = eth_pairings_cpp::run(&data);
    match (native, cpp) {
        (Ok(n), Ok(c)) => {
            if n != c {
                // println!("Native result = {}, C++ result = {}", hex::encode(&n), hex::encode(&c));
                panic!("Native result = {}, C++ result = {}", hex::encode(&n), hex::encode(&c));
            } else {
                // println!("Native and C++ results coincide on {}", hex::encode(&n));
            }
        },
        (Err(n), Err(c)) => {
            // println!("Native and C++ results coincide on error: {:?}, {:?}", n, c);
        },
        (Ok(n), Err(c)) => {
            // println!("Input = {}", hex::encode(&data));
            // println!("Native result = {}, while C++ returned error {:?}", hex::encode(&n), c);
            panic!("Native result = {}, while C++ returned error {:?}", hex::encode(&n), c);
        },
        (Err(n), Ok(c)) => {
            // println!("Input = {}", hex::encode(&data));
            // println!("Native result returned error {:?}, while C++ returned {}", n, hex::encode(&c));
            panic!("Native result returned error {:?}, while C++ returned {}", n, hex::encode(&c));
        }
    }
});
