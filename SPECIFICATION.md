# Specification

Universality of the precompile requires to think about edge cases. Proper splitting of structure ensures that if initial parameters such as field modulus non-residues for extension, etc. are passed correctly, then all the remaining arithmetic is well-defined

## ABI and parsing checks

One precompile is intended to provide a variety of operations under the same address depending from the passed parameters. Full set of operations is the following:
```
pub const OPERATION_ENCODING_LENGTH: usize = 1;

pub const OPERATION_G1_ADD: u8 = 0x01;
pub const OPERATION_G1_MUL: u8 = 0x02;
pub const OPERATION_G1_MULTIEXP: u8 = 0x03;

pub const OPERATION_G2_ADD: u8 = 0x04;
pub const OPERATION_G2_MUL: u8 = 0x05;
pub const OPERATION_G2_MULTIEXP: u8 = 0x06;

pub const OPERATION_PAIRING: u8 = 0x07;
```

After operation type is defined by the first byte of the input the rest is passed to the corresponding function for a specific operation.

G1 additions, multiplications and multiexponentiations are defined for any curve in the Weierstrass form with `b != 0`. Operations in G2 are operations over the curve defined over some extension field. There are only two such extensions supported: degree 2 and degree 3.

```
pub const EXTENSION_DEGREE_ENCODING_LENGTH: usize = 1;
pub const EXTENSION_DEGREE_2: u8 = 0x02;
pub const EXTENSION_DEGREE_3: u8 = 0x03;
```

Pairing operation is defined only for the following families of curves:
- BN
- BLS12
- MNT4
- MNT6
- Ate pairing for a generic curve in the Weierstrass form with `k = 6` and extension tower `Fp - Fp3 - Fp6` (close sibling of MNT6)

### Common ABI for G1 operations

ABI usually consists of two parts: one defines a base field and a curve, with another is operation dependent and encodes points or scalars for a corresponding operations. 

Signature of all the public functions is just `Operation(&[u8])`, so just a pointer to the array of bytes is passes to the function. Remember that operation type is already stripped from the byte array.

```
pub const BYTES_FOR_LENGTH_ENCODING: usize = 1;
```

Important! `take(N)` operation comsumes(!) first `N` bytes from the byte array

Algorithm:
- ensure that length of the byte array is `> BYTES_FOR_LENGTH_ENCODING`, otherwise return error
- `take(BYTES_FOR_LENGTH_ENCODING)` and parse it as unsigned integer `modulus_length` for encoding of length of the next value. Limit `BYTES_FOR_LENGTH_ENCODING = 1` ensures that byte length of the next parameter that is modulus of the prime field is bounded (not more than 255 bytes), so it's a first sanity check
- ensure that `modulus_length > 0`, otherwise return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is a modulus of a base prime field `base_field_modulus`. Ensure that top (most significant) byte is non-zero. This is not an attack and will only require caller to pay more gas, but it's a trivial check
- ensure that `base_field_modulus >= 3` and `base_field_modulus` is odd. This is the second sanity check and also guarantees that Montgommery form that is used for all the field elements is well-defined (is requires `gcd(modulus, R) == 1` with `R` being power of two in our cases). There is no primarity testing, but arithmetic operations now will not trigger panics
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `A` coefficient for a curve in the Weierstrass form
- ensure that `A < base_field_modulus`, otherwise return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `B` coefficient for a curve in the Weierstrass form
- ensure that `B < base_field_modulus`, otherwise return error

This list is a naive algorithm. Real implementation will merge e.g. reading of `A`, range check and parsing it as a representation of some form of the field element in one operation.

Arithmetic is made using some fixed-length representation of a field element. This implementation follows an approach to represent them as a fixed length array `[u64; NUM_LIMBS`], such that for a modulus `M`: `2^(MODULUS_LIMGS*64) > 2*M` to ensure that one never has to take care about carry flags. In this case a field element with `255` bit modulus would be represented as `[u64; 4]`, but `256` bit modulus will be already represented as `[u64; 5]`

To put some sane limit of the length of the arithmetics modulus `M` must be `< 1024` bits in length, so it's representable as `[u64; 16]` array







