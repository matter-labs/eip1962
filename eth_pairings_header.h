#ifndef eth_pairings_bindings_h
#define eth_pairings_bindings_h

#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define MAX_GROUP_BYTE_LEN 128

#define MAX_MODULUS_BYTE_LEN 128

#define PREALLOCATE_FOR_ERROR_BYTES 256

#define PREALLOCATE_FOR_RESULT_BYTES 768

#define SIGN_ENCODING_LENGTH 1

#define SIGN_MINUS 1

#define SIGN_PLUS 0

#define TWIST_TYPE_D 2

#define TWIST_TYPE_LENGTH 1

#define TWIST_TYPE_M 1

#define G1ADD 1

#define G1MUL 2

#define G1MULTIEXP 3

#define G2ADD 4

#define G2MUL 5

#define G2MULTIEXP 6

#define BLS12PAIR 7

#define BNPAIR 8

#define MNT4PAIR 9

#define MNT6PAIR 10

uint32_t c_perform_operation(char op,
                             const char *i,
                             uint32_t i_len,
                             char *o,
                             uint32_t *o_len,
                             char *err,
                             uint32_t *char_len);

#endif /* eth_pairings_bindings_h */