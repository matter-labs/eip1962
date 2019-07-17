#[macro_use] extern crate honggfuzz;
extern crate eth_pairings;
extern crate eth_pairings_cpp;
extern crate hex;

fn main() {
    // Here you can parse `std::env::args and 
    // setup / initialize your project

    // You have full control over the loop but
    // you're supposed to call `fuzz` ad vitam aeternam
    loop {
        // The fuzz macro gives an arbitrary object (see `arbitrary crate`)
        // to a closure-like block of code.
        // For performance reasons, it is recommended that you use the native type
        // `&[u8]` when possible.
        // Here, this slice will contain a "random" quantity of "random" data.
        fuzz!(|data: &[u8]| {
            let native = eth_pairings::public_interface::API::run(&data);
            let cpp = eth_pairings_cpp::run(&data);
            if native.is_err() {
                if !cpp.is_err() {
                    let c = cpp.unwrap();
                    println!("Native result returned error, while C++ returned {}", hex::encode(&c));

                    panic!("Native result returned error, while C++ returned {}", hex::encode(&c));
                }
            } else {
                let n = native.expect("result");
                let c = cpp.expect("cpp result");
                if n != c {
                    println!("Native result = {}, C++ result = {}", hex::encode(&n), hex::encode(&c));

                    panic!("Native result = {}, C++ result = {}", hex::encode(&n), hex::encode(&c));
                }
            }
        });
    }
}