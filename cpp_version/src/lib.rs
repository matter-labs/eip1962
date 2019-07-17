mod api;

pub fn run(bytes: &[u8]) -> Result<Vec<u8>, ()> {
    let mut result = vec![0u8; 4096];
    let input = bytes.as_ptr() as *const std::os::raw::c_char;
    let output = result.as_mut_ptr() as *mut std::os::raw::c_char;
    let input_len = bytes.len() as std::os::raw::c_int;

    let result_len = unsafe { self::api::run(input, input_len, output) };
    if result_len == 0 {
        return Err(());
    }

    result.truncate(result_len as usize);

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    extern crate hex;

    #[test]
    fn test() {
        let result = run(&vec![1u8]);
        println!("Result = {:?}", result);
    }

    #[test]
    fn run_on_hongg_input() {
        use hex;
        let filename = "SIGABRT.EXC_CRASH.PC.00007fff7543a2c6.STACK.000000035c4d6d60.ADDR.0000000000000000.fuzz";
        // use std::time::Instant;
        use std::io::Read;
        use std::fs::File;
        let mut file = File::open(&format!("../honggfuzz/hfuzz_workspace/fuzz_target_compare/{}", filename)).expect("must open");
        let mut input_data = vec![];
        file.read_to_end(&mut input_data).expect("must read");
        assert!(input_data.len() != 0);
        println!("Input = {}", hex::encode(&input_data));
        // let now = Instant::now();
        let result = run(&input_data[..]);
        println!("Result = {:?}", result);

        // let elapsed = now.elapsed().as_micros();
        // let gas_estimate = crate::gas_meter::GasMeter::meter(&input_data[..]);
        // if result.is_err() {
        //     println!("Api call failed in {} micros", elapsed);
        //     if gas_estimate.is_ok() {
        //         println!("Gas estimate was {}", gas_estimate.unwrap());
        //     } else {
        //         println!("Gas estimate failed");
        //     }
        // } else {
        //     println!("Api call was ok in {} micros", elapsed);
        //     if gas_estimate.is_ok() {
        //         println!("Gas estimate was {}", gas_estimate.unwrap());
        //     } else {
        //         println!("Gas estimate failed");
        //     }
        // }
    }
}