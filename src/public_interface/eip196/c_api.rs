// For C style API caller has to preallocate some buffers for results 
pub const EIP196_PREALLOCATE_FOR_ERROR_BYTES: usize = 256;
pub const EIP196_PREALLOCATE_FOR_RESULT_BYTES: usize = 32 * 2; // maximum for G2 point

use static_assertions::const_assert;
const_assert!(EIP196_PREALLOCATE_FOR_RESULT_BYTES == super::SERIALIZED_G1_POINT_BYTE_LENGTH);

#[allow(non_camel_case_types)]
#[repr(u8)]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Eip196OperationType {
    ADD = 1,
    MUL = 2,
    PAIR = 3,
}

impl Eip196OperationType {
    pub fn from_u8(value: u8) -> Option<Self> {
        match value {
            EIP196_ADD_OPERATION_RAW_VALUE => {
                Some(Eip196OperationType::ADD)
            },
            EIP196_MUL_OPERATION_RAW_VALUE => {
                Some(Eip196OperationType::MUL)
            },
            EIP196_PAIR_OPERATION_RAW_VALUE => {
                Some(Eip196OperationType::PAIR)
            },
            _ => {
                None
            }
        }
    }

    pub fn as_u8(&self) -> u8 {
        *self as u8
    }
}

pub const EIP196_ADD_OPERATION_RAW_VALUE: u8 = Eip196OperationType::ADD as u8;
pub const EIP196_MUL_OPERATION_RAW_VALUE: u8 = Eip196OperationType::MUL as u8;
pub const EIP196_PAIR_OPERATION_RAW_VALUE: u8 = Eip196OperationType::PAIR as u8;

// this is C interface
#[no_mangle]
pub extern "C" fn eip196_perform_operation(
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
    let err_out_i8: &mut [libc::c_char] = unsafe { std::slice::from_raw_parts_mut(err, EIP196_PREALLOCATE_FOR_ERROR_BYTES) };
    let mut err_out: &mut [u8] = unsafe { std::mem::transmute(err_out_i8) };

    let operation = Eip196OperationType::from_u8(op_u8);

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
    
    let input_i8: & [libc::c_char] = unsafe { std::slice::from_raw_parts(i, i_len as usize) };
    let input: &[u8] = unsafe { std::mem::transmute(input_i8) };

    let raw_out_i8: &mut [libc::c_char] = unsafe { std::slice::from_raw_parts_mut(o, EIP196_PREALLOCATE_FOR_RESULT_BYTES) };
    let mut raw_out: &mut [u8] = unsafe { std::mem::transmute(raw_out_i8) };

    let result = match operation {
        Eip196OperationType::ADD => super::EIP196Executor::add(&input).map(|r| r[..].to_vec()),
        Eip196OperationType::MUL => super::EIP196Executor::mul(&input).map(|r| r[..].to_vec()),
        Eip196OperationType::PAIR => super::EIP196Executor::pair(&input).map(|r| r[..].to_vec()),
    };

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
