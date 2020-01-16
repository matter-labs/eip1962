# Integration

## Exposed interface

Both C++ implementation (from [here](https://github.com/matter-labs/eip1962_cpp), folder `eip1962cpp`) and Rust (with `features = ["c_api"]`) expose the following C interface
```
uint32_t c_perform_operation(char op,
                             const char *i,
                             uint32_t i_len,
                             char *o,
                             uint32_t *o_len,
                             char *err,
                             uint32_t *char_len);


uint32_t c_meter_operation(char op,
                             const char *i,
                             uint32_t i_len,
                             uint64_t *gas,
                             char *err,
                             uint32_t *char_len);
```

(`c_meter_operation` is not yet exposed by Rust);

There are also constants describing the length of preallocated byte arrays.

```
pub const PREALLOCATE_FOR_ERROR_BYTES: usize = 256;
pub const PREALLOCATE_FOR_RESULT_BYTES: usize = 768;
```

Example of integration using Golang (with CGO)

```

package eip1962

/*
#cgo CXXFLAGS: -std=c++17
#cgo CXXFLAGS: -I./include
#include "wrapper.h"
*/
import "C"

import (
	"errors"
	"unsafe"
)

const maxErrLen = 256
const maxOutputLen = 768

var (
	ErrInvalidMsgLen = errors.New("invalid data length, need >= bytes")
	ErrCallFailed    = errors.New("library call returned an error")
)

// Call calls the C++ implementation of the EIP
func Call(operation uint8, data []byte) ([]byte, error) {
	if len(data) == 0 {
		return nil, ErrInvalidMsgLen
	}
	ilen := len(data)
	outputBytes := make([]byte, maxOutputLen)
	olen := uint32(0)
	errStringBytes := make([]byte, maxOutputLen)
	errStringLen := uint32(0)

	var (
        op         = (C.char)(operation)
		inputdata  = (*C.char)(unsafe.Pointer(&data[0]))
		inputlen   = (C.uint32_t)(ilen)
		outputdata = (*C.char)(unsafe.Pointer(&outputBytes[0]))
		outputlen  = (*C.uint32_t)(unsafe.Pointer(&olen))
		errdata    = (*C.char)(unsafe.Pointer(&errStringBytes[0]))
		errlen     = (*C.uint32_t)(unsafe.Pointer(&errStringLen))
	)

	result := C.run(op, inputdata, inputlen, outputdata, outputlen, errdata, errlen)
	if result == 1 {
		// parse error string
		return nil, ErrCallFailed
	}

	return outputBytes[:olen], nil
}

```