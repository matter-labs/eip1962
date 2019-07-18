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
            match (native, cpp) {
                (Ok(n), Ok(c)) => {
                    if n != c {
                        println!("Native result = {}, C++ result = {}", hex::encode(&n), hex::encode(&c));
                        panic!("Native result = {}, C++ result = {}", hex::encode(&n), hex::encode(&c));
                    } else {
                        println!("Native and C++ results coincide on {}", hex::encode(&n));
                    }
                },
                (Err(n), Err(c)) => {
                    println!("Native and C++ results coincide on error: {:?}, {:?}", n, c);
                },
                (Ok(n), Err(c)) => {
                    println!("Input = {}", hex::encode(&data));
                    println!("Native result = {}, while C++ returned error {:?}", hex::encode(&n), c);
                    panic!("Native result = {}, while C++ returned error {:?}", hex::encode(&n), c);
                },
                (Err(n), Ok(c)) => {
                    println!("Input = {}", hex::encode(&data));
                    println!("Native result returned error {:?}, while C++ returned {}", n, hex::encode(&c));
                    panic!("Native result returned error {:?}, while C++ returned {}", n, hex::encode(&c));
                }
            }
        });
    }
}