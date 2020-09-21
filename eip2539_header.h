#ifndef eip2539_bindings_h
#define eip2539_bindings_h

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define EIP2539_PREALLOCATE_FOR_ERROR_BYTES 256

#define EIP2539_PREALLOCATE_FOR_RESULT_BYTES ((64 * 2) * 2)

#define BLS12_G1ADD_OPERATION_RAW_VALUE 1
#define BLS12_G1MUL_OPERATION_RAW_VALUE 2
#define BLS12_G1MULTIEXP_OPERATION_RAW_VALUE 3
#define BLS12_G2ADD_OPERATION_RAW_VALUE 4
#define BLS12_G2MUL_OPERATION_RAW_VALUE 5
#define BLS12_G2MULTIEXP_OPERATION_RAW_VALUE 6
#define BLS12_PAIR_OPERATION_RAW_VALUE 7

uint32_t eip2539_perform_operation(char op,
                                   const char *i,
                                   uint32_t i_len,
                                   char *o,
                                   uint32_t *o_len,
                                   char *err,
                                   uint32_t *char_len);

#endif /* eip2539_bindings_h */
