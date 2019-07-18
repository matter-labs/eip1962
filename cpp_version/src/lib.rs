mod api;

const MAX_OUTPUT_LEN: usize = 256*3*2; // output if Fp = 1023 bits, Fp3 extension, affine point

pub fn run(bytes: &[u8]) -> Result<Vec<u8>, String> {
    let mut result = vec![0u8; MAX_OUTPUT_LEN];
    let mut error_description_buffer = vec![0u8; MAX_OUTPUT_LEN];
    let input = bytes.as_ptr() as *const std::os::raw::c_char;
    let output = result.as_mut_ptr() as *mut std::os::raw::c_char;
    let error_buffer = error_description_buffer.as_mut_ptr() as *mut std::os::raw::c_char;
    let input_len = bytes.len() as u32;
    let mut output_len = 0u32;
    let mut error_description_len = 0u32;

    let success = unsafe { self::api::run(
        input, 
        input_len, 
        output, 
        &mut output_len as *mut u32,
        error_buffer,
        &mut error_description_len as *mut u32
    ) };
    if success == 0 {
        if error_description_len == 0 {
            return Err("C++ api returned empty error description".to_string());
        }
        error_description_buffer.truncate(error_description_len as usize);
        let error_description_string = std::ffi::CString::new(error_description_buffer);
        match error_description_string {
            Ok(c_string) => {
                let string = c_string.into_string();
                match string {
                    Ok(string) => {
                        return Err(string);
                    },
                    Err(err) => {
                        return Err(format!("Error on conversion of string description, {:?}", err));
                    }
                }
            },
            Err(n_error) => {
                return Err(format!("CString containts empty bytes in a middle, {:?}", n_error));
            }
        }
    }

    result.truncate(output_len as usize);

    Ok(result)
}

pub fn meter(bytes: &[u8]) -> Result<u64, String> {
    let input = bytes.as_ptr() as *const std::os::raw::c_char;
    let input_len = bytes.len() as u32;
    let mut gas = 0u64;

    let success = unsafe { self::api::meter_gas(
        input, 
        input_len, 
        &mut gas as *mut u64
    ) };
    if success == 0 {
        return Err("Failed to meter gas".to_string());
    }

    Ok(gas)
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