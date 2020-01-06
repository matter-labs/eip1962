use super::unified_api::{OperationType, PREALLOCATE_FOR_ERROR_BYTES, ;

// this is C interface
#[no_mangle]
pub extern "C" fn c_perform_operation(
    op: ::std::os::raw::c_char,
    i: *const ::std::os::raw::c_char,
    i_len: u32,
    o: *mut ::std::os::raw::c_char,
    o_len: *mut u32,
    err: *mut ::std::os::raw::c_char,
    char_len: *mut u32) -> u32 
{            
    use std::io::Write;

    let op_u8: u8 = unsafe { std::mem::transmute(op) };
    let err_out_i8: &mut [i8] = unsafe { std::slice::from_raw_parts_mut(err, PREALLOCATE_FOR_ERROR_BYTES) };
    let mut err_out: &mut [u8] = unsafe { std::mem::transmute(err_out_i8) };

    let operation = OperationType::from_u8(op_u8);

    if operation.is_none() {
        let written = err_out.write(b"Unknown operation type\0");
        if let Ok(bytes_written) = written {
            unsafe { *char_len = bytes_written as u32 };
        } else {
            unsafe { *char_len = 0u32 };
        }

        return 1u32;
    }

    let operation = operation.expect("is some");
    
    let input_i8: & [i8] = unsafe { std::slice::from_raw_parts(i, i_len as usize) };
    let input: &[u8] = unsafe { std::mem::transmute(input_i8) };

    let raw_out_i8: &mut [i8] = unsafe { std::slice::from_raw_parts_mut(o, PREALLOCATE_FOR_RESULT_BYTES) };
    let mut raw_out: &mut [u8] = unsafe { std::mem::transmute(raw_out_i8) };

    let result = perform_operation(operation, input);

    match result {
        Ok(result) => {
            let written = raw_out.write(result.as_ref());
            if let Ok(bytes_written) = written {
                unsafe { *o_len = bytes_written as u32 };
                return 0u32;
            }

            let written = err_out.write(b"Failed to write the result\0");
            if let Ok(bytes_written) = written {
                unsafe { *char_len = bytes_written as u32 };
            } else {
                unsafe { *char_len = 0u32 };
            }

            return 1u32;
        },
        Err(error) => {
            let err_description = error.to_string();
            let written = err_out.write(err_description.as_bytes());
            if let Ok(bytes_written) = written {
                unsafe { *char_len = bytes_written as u32 };
            } else {
                unsafe { *char_len = 0u32 };
            }

            return 1u32;
        }
    }
} 