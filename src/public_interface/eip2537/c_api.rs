// For C style API caller has to preallocate some buffers for results 
pub const EIP2537_PREALLOCATE_FOR_ERROR_BYTES: usize = 256;
pub const EIP2537_PREALLOCATE_FOR_RESULT_BYTES: usize = 64 * 2 * 2; // maximum for G2 point

use static_assertions::const_assert;
const_assert!(EIP2537_PREALLOCATE_FOR_RESULT_BYTES == super::SERIALIZED_G2_POINT_BYTE_LENGTH);

#[repr(u8)]
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Eip2537OperationType {
    BLS12_G1ADD = 1,
    BLS12_G1MUL = 2,
    BLS12_G1MULTIEXP = 3,
    BLS12_G2ADD = 4,
    BLS12_G2MUL = 5,
    BLS12_G2MULTIEXP = 6,
    BLS12_PAIR = 7,
    BLS12_MAP = 8,
}

impl Eip2537OperationType {
    pub fn from_u8(value: u8) -> Option<Self> {
        match value {
            BLS12_G1ADD_OPERATION_RAW_VALUE => {
                Some(Eip2537OperationType::BLS12_G1ADD)
            },
            BLS12_G1MUL_OPERATION_RAW_VALUE => {
                Some(Eip2537OperationType::BLS12_G1MUL)
            },
            BLS12_G1MULTIEXP_OPERATION_RAW_VALUE => {
                Some(Eip2537OperationType::BLS12_G1MULTIEXP)
            },
            BLS12_G2ADD_OPERATION_RAW_VALUE => {
                Some(Eip2537OperationType::BLS12_G2ADD)
            },
            BLS12_G2MUL_OPERATION_RAW_VALUE => {
                Some(Eip2537OperationType::BLS12_G2MUL)
            },
            BLS12_G2MULTIEXP_OPERATION_RAW_VALUE => {
                Some(Eip2537OperationType::BLS12_G2MULTIEXP)
            },
            BLS12_PAIR_OPERATION_RAW_VALUE => {
                Some(Eip2537OperationType::BLS12_PAIR)
            },
            BLS12_MAP_OPERATION_RAW_VALUE => {
                Some(Eip2537OperationType::BLS12_MAP)
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

pub const BLS12_G1ADD_OPERATION_RAW_VALUE: u8 = Eip2537OperationType::BLS12_G1ADD as u8;
pub const BLS12_G1MUL_OPERATION_RAW_VALUE: u8 = Eip2537OperationType::BLS12_G1MUL as u8;
pub const BLS12_G1MULTIEXP_OPERATION_RAW_VALUE: u8 = Eip2537OperationType::BLS12_G1MULTIEXP as u8;

pub const BLS12_G2ADD_OPERATION_RAW_VALUE: u8 = Eip2537OperationType::BLS12_G2ADD as u8;
pub const BLS12_G2MUL_OPERATION_RAW_VALUE: u8 = Eip2537OperationType::BLS12_G2MUL as u8;
pub const BLS12_G2MULTIEXP_OPERATION_RAW_VALUE: u8 = Eip2537OperationType::BLS12_G2MULTIEXP as u8;

pub const BLS12_PAIR_OPERATION_RAW_VALUE: u8 = Eip2537OperationType::BLS12_PAIR as u8;
pub const BLS12_MAP_OPERATION_RAW_VALUE: u8 = Eip2537OperationType::BLS12_MAP as u8;

// this is C interface
#[no_mangle]
pub extern "C" fn eip2537_perform_operation(
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
    let err_out_i8: &mut [i8] = unsafe { std::slice::from_raw_parts_mut(err, EIP2537_PREALLOCATE_FOR_ERROR_BYTES) };
    let mut err_out: &mut [u8] = unsafe { std::mem::transmute(err_out_i8) };

    let operation = Eip2537OperationType::from_u8(op_u8);

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

    let raw_out_i8: &mut [i8] = unsafe { std::slice::from_raw_parts_mut(o, EIP2537_PREALLOCATE_FOR_RESULT_BYTES) };
    let mut raw_out: &mut [u8] = unsafe { std::mem::transmute(raw_out_i8) };

    let result = match operation {
        Eip2537OperationType::BLS12_G1ADD => super::EIP2537Executor::g1_add(&input).map(|r| r[..].to_vec()),
        Eip2537OperationType::BLS12_G1MUL => super::EIP2537Executor::g1_mul(&input).map(|r| r[..].to_vec()),
        Eip2537OperationType::BLS12_G1MULTIEXP => super::EIP2537Executor::g1_multiexp(&input).map(|r| r[..].to_vec()),
        Eip2537OperationType::BLS12_G2ADD => super::EIP2537Executor::g2_add(&input).map(|r| r[..].to_vec()),
        Eip2537OperationType::BLS12_G2MUL => super::EIP2537Executor::g2_mul(&input).map(|r| r[..].to_vec()),
        Eip2537OperationType::BLS12_G2MULTIEXP => super::EIP2537Executor::g2_multiexp(&input).map(|r| r[..].to_vec()),
        Eip2537OperationType::BLS12_PAIR => super::EIP2537Executor::pair(&input).map(|r| r[..].to_vec()),
        Eip2537OperationType::BLS12_MAP => super::EIP2537Executor::map_to_curve(&input).map(|r| r[..].to_vec()),
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