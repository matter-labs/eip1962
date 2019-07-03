# Specification

Universality of the precompile requires to think about edge cases. Proper splitting of structure ensures that if initial parameters such as field modulus non-residues for extension, etc. are passed correctly, then all the remaining arithmetic is well-defined. This document can be viewed as both API spec, implementation guide and use guide for some aspects

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

Signature of all the public functions is just `Operation(byte_array)`, so just a pointer to the array of bytes is passes to the function. Remember that operation type is already stripped from the byte array.

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
- `take(BYTES_FOR_LENGTH_ENCODING)` and parse it as unsigned integer `subgroup_order_length` for encoding of length of the next value
- ensure that `subgroup_order_length > 0`, otherwise return error
- ensure that length of the byte array is `> subgroup_order_length`, otherwise return error
- `take(subgroup_order_length)` and parse it as BigEndian encoding of an unsigned integer `main_subgroup_order`
- ensure that `main_subgroup_order >= 3` and is odd as a sanity check
- there are two purposes to requure `main_subgroup_order` to be included into the ABI:
  - Upfront estimation of the worst-case scenario for a difficulty of the multiplication and multiexponentiation operations
  - One should not enforce point being in a correct subgroup for correctness of operations and if caller expects to process user's inputs he should make such check as a separate call. Nevertheless, this check is MANDATORY for pairing operations
- later should follow operation-specific parameters
  
This list is a naive algorithm. Real implementation will merge e.g. reading of `A`, range check and parsing it as a representation of some form of the field element in one operation.

Arithmetic is made using some fixed-length representation of a field element. This implementation follows an approach to represent them as a fixed length array `[u64; NUM_LIMBS`], such that for a modulus `M`: `2^(MODULUS_LIMGS*64) > 2*M` to ensure that one never has to take care about carry flags. In this case a field element with `255` bit modulus would be represented as `[u64; 4]`, but `256` bit modulus will be already represented as `[u64; 5]`

To put some sane limit of the length of the arithmetics modulus `M` must be `< 1024` bits in length, so it's representable as `[u64; 16]` array

#### Specific ABI part for addition in G1

At this point one knows `modulus_length` and `base_field_modulus`, so now one can parse point for an addition operation

- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `x` of the point `P0`
- ensure that `P0.x < base_field_modulus`, otherwise return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `y` of the point `P0`
- ensure that `P0.y < base_field_modulus`, otherwise return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `x` of the point `P1`
- ensure that `P1.x < base_field_modulus`, otherwise return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `y` of the point `P1`
- ensure that `P1.y < base_field_modulus`, otherwise return error
- at this point curve is well-defined, along with arithmetic on it
- check that `is_on_curve(P0)` and `is_on_curve(P1)`, otherwise return error
- perform an addition depending on a form of the Weierstrass curve (`A == 0` or `A != 0`, `B` is always non-zero) `P_res = P0 + P1`
- operations are most likely to be performed in Jacobial coordinates, so perform normalization into affine coordinates. This require to make inversion of `P_res.z`. If `P_res.z == 0` return point of infinity (this is an expected result of the addition operation in Jacobian coordinates for `P` and `-P`), otherwise inverse `P_res.z` and perform normalization

#### Specific ABI part for multiplication in G1

At this point one knows `modulus_length` and `base_field_modulus`, so now one can parse point for an addition operation

- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `x` of the point `P0`
- ensure that `P0.x < base_field_modulus`, otherwise return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `y` of the point `P0`
- ensure that `P0.y < base_field_modulus`, otherwise return error
- ensure that length of the byte array is `>= subgroup_order_length`, otherwise return error
- `take(subgroup_order_length)` and parse it as BigEndian encoding of an unsigned integer that is an `scalar` for multiplication operation
- ensure that `scalar <= main_subgroup_order`, otherwise return error. Scalar can be equal to the group order for usefull operations for caller
- perform a multiplication depending on a form of the Weierstrass curve (`A == 0` or `A != 0`, `B` is always non-zero) `P_res = scalar*P0`
- operations are most likely to be performed in Jacobial coordinates, so perform normalization into affine coordinates. This require to make inversion of `P_res.z`. If `P_res.z == 0` return point of infinity, otherwise inverse `P_res.z` and perform normalization

Possible variations: explicitly encode length of the BigEndian byte representation of scalar (difficult for gas estimations), or use same level of granularity for any scalar that is an integer number of limbs. In this case one would also want to change common ABI part to allow zero byte padding of a main subgroup order

#### Specific ABI part for multiexponentiation in G1

At this point one knows `modulus_length` and `base_field_modulus`, so now one can parse point for an addition operation

- calculate expected byte length of one `(point, scalar)` pair as `expected_pair_len = 2*modulus_length + subgroup_order_length`
- ensure that length of the byte array is a multiple of `expected_pair_len`, otherwise return error
- calculate number of pairs by dividing the rest of the byte array by `expected_pair_len`
- in a loop:
  - `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `x` of the point `P0`
  - ensure that `P0.x < base_field_modulus`, otherwise return error
  - ensure that length of the byte array is `> modulus_length`, otherwise return error
  - `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `y` of the point `P0`
  - ensure that `P0.y < base_field_modulus`, otherwise return error
  - `take(subgroup_order_length)` and parse it as BigEndian encoding of an unsigned integer that is an `scalar` for multiplication operation
  - ensure that `scalar <= main_subgroup_order`, otherwise return error
- perform a miltiexponentiation depending on a form of the Weierstrass curve (`A == 0` or `A != 0`, `B` is always non-zero), output `P_res`
- operations are most likely to be performed in Jacobial coordinates, so perform normalization into affine coordinates. This require to make inversion of `P_res.z`. If `P_res.z == 0` return point of infinity, otherwise inverse `P_res.z` and perform normalization

#### Notes for G1 ABI:

For purposes of caller's convenience it may be reasonable to have some granularity for lengths of encodings of various element. For example, make encodings multiples of 8 bytes (64 bits) to roughtly correspond to limb bits on x64 machine

### Common ABI for G2 operations

Operations on a "twist" are defined and expected to be used for pairing friendly curves. E.G. original protocol of BLS aggregated signatures requires multiplication in G2, as well as some SNARK verification equations.

To save space only common ABI part for G2 is described, with specific part being similar to G1 part.
```
pub const EXTENSION_DEGREE_ENCODING_LENGTH: usize = 1;
```

Algorithm:
- ensure that length of the byte array is `> BYTES_FOR_LENGTH_ENCODING`, otherwise return error
- `take(BYTES_FOR_LENGTH_ENCODING)` and parse it as unsigned integer `modulus_length` for encoding of length of the next value. Limit `BYTES_FOR_LENGTH_ENCODING = 1` ensures that byte length of the next parameter that is modulus of the prime field is bounded (not more than 255 bytes), so it's a first sanity check
- ensure that `modulus_length > 0`, otherwise return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is a modulus of a base prime field `base_field_modulus`. Ensure that top (most significant) byte is non-zero. This is not an attack and will only require caller to pay more gas, but it's a trivial check
- ensure that `base_field_modulus >= 3` and `base_field_modulus` is odd. This is the second sanity check and also guarantees that Montgommery form that is used for all the field elements is well-defined (is requires `gcd(modulus, R) == 1` with `R` being power of two in our cases). There is no primarity testing, but arithmetic operations now will not trigger panics
- ensure that length of the byte array is `> EXTENSION_DEGREE_ENCODING_LENGTH`, otherwise return error
- `take(EXTENSION_DEGREE_ENCODING_LENGTH)` and parse it as unsigned integer `extension_degree` for encoding of extension degree for a twist. Only `extension_degree == 2` or `extension_degree == 3` are supported, other values should return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `non_residue` - a quadratic or cubic non-residue to make an extension
- check that `non_residue` is non-square or non-cube depending of `extension_degree` to have extension well-formed. If it's not - return error
- - ensure that length of the byte array is `> modulus_length*extension_degree`, otherwise return error
- `take(modulus_length*extension_degree)` and parse it as `extension_degree` densely packed BigEndian encodings of unsigned integers that are coefficient of an element in extension field. Coefficients follow from smallest degree: if element is represented as a polynomial `c0 + c1*x + c2*x^2` then coefficients are parsed as `c0`, `c1`, `c2` one after another. That is an `A` coefficient for a curve twist in the Weierstrass form
- ensure that each of `c*` coefficientsi `< base_field_modulus`, otherwise return error
- ensure that length of the byte array is `> modulus_length*extension_degree` and perform similar checks for `B` coefficient
- `take(BYTES_FOR_LENGTH_ENCODING)` and parse it as unsigned integer `subgroup_order_length` for encoding of length of the next value
- ensure that `subgroup_order_length > 0`, otherwise return error
- ensure that length of the byte array is `> subgroup_order_length`, otherwise return error
- `take(subgroup_order_length)` and parse it as BigEndian encoding of an unsigned integer `main_subgroup_order`
- ensure that `main_subgroup_order >= 3` and is odd as a sanity check

Operations in G2 are the same as for G1 with a difference only in encoding of point coordinates now being in the extension of the base field.

### Common ABI for pairing operations

Due to difference in properties and required parameters for different families of the curves first one has to parse curve type:

```
pub const CURVE_TYPE_LENGTH: usize = 1;
pub const BLS12: u8 = 0x01;
pub const BN: u8 = 0x02;
pub const MNT4: u8 = 0x03;
pub const MNT6: u8 = 0x04;
```

For now encoding for generic Ate pairing for `k=6` extension is not assigned.

- ensure that length of the byte array is `> CURVE_TYPE_LENGTH`, otherwise return error
- `take(CURVE_TYPE_LENGTH)` and parse it as unsigned integer `curve_type` for encoding of the curve type. Check that curve type is known an follow the part that is curve specific

### ABI for pairing operations on BLS12 curves

Important constants:
```
pub const TWIST_TYPE_LENGTH: usize = 1;
pub const TWIST_TYPE_M: u8 = 0x01;
pub const TWIST_TYPE_D: u8 = 0x02;

pub const SIGN_ENCODING_LENGTH: usize = 1;
pub const SIGN_PLUS: u8 = 0x00;
pub const SIGN_MINUS: u8 = 0x01;
```

Note that BLS12 is a family of curves that are parametrized by a single scalar `x`, twist type that is either `M` (multiplication) or `D` (division), and structure of the extension tower (non-residues)

- parse `modulus_length`, `base_field_modulus`, `main_subgroup_order` and `A` and `B` coefficients following the G1 ABI
- check that `A==0` (true for BLS12 curves)
- parse `non_residue_for_fp2` (that is an element of the base field) that is used to construct `Fp2` extension following the logic described for G2 ABI. Check that it's non-square, otherwise return error
- parse `non_residue_for_fp6` (that is an element of `Fp2`) that is used to construct `Fp6` extension following the logic described in G2 ABI for parsing `Fp2` elements. Check that it's non-cube in `Fp2`, otherwise return error
- ensure that length of the byte array is `> TWIST_TYPE_LENGTH`, otherwise return error
- `take(TWIST_TYPE_LENGTH)` and parse it as unsigned integer `twist_type` for encoding of the curve twist type. Check that twist type is either `M` or `D`, otherwise return error
- ensure that length of the byte array is `> BYTES_FOR_LENGTH_ENCODING`, otherwise return error
- `take(BYTES_FOR_LENGTH_ENCODING)` and parse it as unsigned integer `x_length` for encoding of length of the next value
- ensure that length of the byte array is `> x_length`, otherwise return error
- `take(x_length)` and parse it as unsigned integer `x`. Empirical testing will put a bound on the sane limits for `x`, that is to be determined
- ensure that length of the byte array is `> SIGN_ENCODING_LENGTH`, otherwise return error
- `take(SIGN_ENCODING_LENGTH)` and parse it as unsigned integer `x_sign` for encoding of the sign of `x`. Check that it's either `SIGN_PLUS` or `SIGN_MINUS`
- ensure that length of the byte array is `> BYTES_FOR_LENGTH_ENCODING`, otherwise return error
- `take(BYTES_FOR_LENGTH_ENCODING)` and parse it as unsigned integer `num_pairs` for encoding number of point pairing for an operation. This is handy for gas estimation
- calculate expected byte length of encoding of the points pair as `6*modulus_length` and ensure that length of rest of the byte array is equal to `6*modulus_length*num_pairs`
- parse `num_pairs` point pairs as `(G1_point, G2_point)` tuples encoded following rules described in G1 and G2 ABI
- for each point perform on-curve check, otherwise return an error
- for each point perform subgroup check (by multiplying by main subgroup order and comparing to point of infinity). IMPORTANT: while this check is expensive, it's mandatory to ensure correctness of the pairing operation
- perform pairing following formulas for BLS12 curves. Such operation requires inversions that may not exist, in this case return as error
- if result of a pairing (element of `Fp12`) is equal to identity - return single byte `0x01`, otherwise return `0x00` following the existing ABI for BN254 precompile

While BLS12 curve has 

Reader may note that we never got coefficients `A_fp2` and `B_fp2` for curve twist. Depending on a twist type (`M` or `D`) `B_fp2` = `B*non_residue_for_fp6` or `B_fp2` = `B/non_residue_for_fp6` respectively

#### ABI for other curve families

To be extended

### General notices for implementation

While the specification described above guarantees that all field operations are well-defined, one still has two important points to discuss:
- Implementation of gas estimation based on that parsed parameters. One should grab only the minimal number of parameters required for estimation and perform an estimation as soon as possible. For BLS12 pairing (that has the most number of parameters in this spec) one should parform parsing for encoding correctness and delay computations (even such as on-curve checks and that non-residues are valid) until gas is estimated based on the following parameters:
  -  `base_field_modulus` - will give number of limbs required for implementation of field operations, with most of operations having quadratic complexity in number of limbs. Will also determine difficulty of the final exponentiation and some precomputations giving higher-order terms in gas schedule
  -  `main_subgroup_order` - will determine difficulty of subgroup checks
  -  `x` - specifically hamming weight and bit count will determine that difficulty of the Miller loop
- "Garbage in - garbase out". While any valid implementation will give the same results on the same valid input parameters (that define a legitimate curve), caller can supply random values that would pass the ABI checks and gas estimations, but be meaningless in a sense of operations performed. Such input is called "garbage". It's not immediately clean whether or not different implementations will return the same output in this case. This is important for consensus between different implementations and there are two ways to resolve this issue:
  - Use single reference maintained implementation between clients. As history has shown there were bugs in major Ethereum clients in terms of BN254 precompile implementations: Parity didn't properly called the Rust code from Zcash and had excessive exponentiations and insane running time, Geth had no subgroup checks initially, and Pantheon still missed subgroup checks. Implementation of cryptographic primitives is difficult and it's reasonable to focus on universal implementation for all the major clients, for example as a Rust implementation with C ABI. This automatically eliminates problem of output on "garbage" input for consensus purporses. Other implementations are still necessary for cross-check in "valid" scenarios and should be as independent from the reference implementation as possible
  - Have few independent implementation following the same formulas for different operations. Such formulas are well-known and (if validity checks from ABI sections are performed) will give the same results for the same input. This approach requires extra implementation effort, but consistency for "valid" cases again should be cross-checked against independent "distant" implementations

