#ifndef eip196_bindings_h
#define eip196_bindings_h

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define EIP196_PREALLOCATE_FOR_ERROR_BYTES 256

#define EIP196_PREALLOCATE_FOR_RESULT_BYTES 64

#define EIP196_ADD_OPERATION_RAW_VALUE 1
#define EIP196_MUL_OPERATION_RAW_VALUE 2
#define EIP196_PAIR_OPERATION_RAW_VALUE 3

uint32_t eip196_perform_operation(char op,
                                   const char *i,
                                   uint32_t i_len,
                                   char *o,
                                   uint32_t *o_len,
                                   char *err,
                                   uint32_t *char_len);

#endif /* eip196_bindings_h */
