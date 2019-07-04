# ABI speficiation

## Supported operations

The precompile provides multiple elliptic curve operations. The full set of operations is defined as follows:

|Operation            |Code|
|---------------------|----|
|OPERATION_G1_ADD     |0x01|
|OPERATION_G1_MUL     |0x02|
|OPERATION_G1_MULTIEXP|0x03|
|OPERATION_G2_ADD     |0x04|
|OPERATION_G2_MUL     |0x05|
|OPERATION_G2_MULTIEXP|0x06|
|OPERATION_PAIRING    |0x07|

`OPERATION_G1_ADD`, `OPERATION_G1_MUL` and `OPERATION_G1_MULTIEXP` are operations of additon, multiplication and multiexponentiation for G1 elements of any curve in the Weierstrass form with `b != 0`.

`OPERATION_G2_ADD`, `OPERATION_G2_MUL` and `OPERATION_G2_MULTIEXP` are operations for G2 elements for a curve defined over some extension field.

`OPERATION_PAIRING` is the pairing operation. The following curve families are supported:

- BN
- BLS12
- MNT4
- MNT6
- Ate pairing for a generic curve in the Weierstrass form with `k = 6` and extension tower `Fp - Fp3 - Fp6` (close sibling of MNT6)

## Precompile input (call data)

Call data must be a correctly encoded ABI data string of two elements:

|Value  |Type       |Length  |
|-------|-----------|--------|
|op_code|uint8      |1 byte  |
|op_data|bytes_array|variable|

The first byte of the input specifies the type of the operation. The remaining data is passed to the corresponding operation handler (see details below).

All numbers are passed in **big endian** encoding.

Incorrect data input is always handled and returns an error.

## op_data for G1 operations

`op_data` for all G1 operations consists of a common prefix followed by the operands.

The common prefix must have the following form:

|Value              |Length                    |Comment                    |
|-------------------|--------------------------|---------------------------|
|field_length       |1 byte                    |                           |
|modulus            |`field_length` bytes      |Fq modulus                 |
|a                  |`field_length` bytes      |Curve's a coefficient      |
|b                  |`field_length` bytes      |Curve's b coefficient      |
|group_order_length |1 bytes                   |                           |                    
|group_order        |`group_order_length` bytes|Group order                |

The operands are described below for each operation.

### OPERATION_G1_ADD operands

|Value              |Length                    |                                  |
|-------------------|--------------------------|----------------------------------|
|lhs                |`2*field_length` bytes    |First point's X and Y coordinates |
|rhs                |`2*field_length` bytes    |Second point's X and Y coordinates|

### OPERATION_G1_MUL operands

|Value              |Length                    |                                  |
|-------------------|--------------------------|----------------------------------|
|lhs                |`2*field_length` bytes    |First point's X and Y coordinates |
|rhs                |`group_order_length` bytes|Sсalar multiplication factor      |

### OPERATION_G1_MULTIEXP operands

The multiexponentiation operation can take arbitrary number of operands. Each of the operands must be encoded in the following form:

|Value              |Length                    |                                  |
|-------------------|--------------------------|----------------------------------|
|point              |`2*field_length` bytes    |Point's X and Y coordinates       |
|scalar             |`group_order_length` bytes|Sсalar order of exponentiation    |


## op_data for G2 operations

`op_data` for all G2 operations consists of a common prefix followed by the operands.

The common prefix must have the following form:

|Value              |Length                    |Comment                       |
|-------------------|--------------------------|------------------------------|
|field_length       |1 byte                    |                              |
|modulus            |`field_length` bytes      |Fq modulus                    |
|extension_degree   |1 bytes                   |Only values 2 or 3 are allowed|
|fp_non_residue     |`field_length` bytes      |Fp2 non-residue               |
|a                  |`field_length` bytes      |Curve's a coefficient         |
|b                  |`field_length` bytes      |Curve's b coefficient         |
|group_order_length |1 bytes                   |                              |                    
|group_order        |`group_order_length` bytes|Group order                   |

The operands are described below for each operation. They follow the same schema as for G1 operations, except that all points are encoded in the required extension degree.

### OPERATION_G2_ADD operands

|Value              |Length                                   |                                                          |
|-------------------|-----------------------------------------|----------------------------------------------------------|
|lhs                |`extension_degree*field_length` bytes    |First point's coordinates in the extension field          |
|rhs                |`extension_degree*field_length` bytes    |Second point's X and Y coordinates in the extension field |

### OPERATION_G2_MUL operands

|Value              |Length                                   |                                                         |
|-------------------|-----------------------------------------|---------------------------------------------------------|
|lhs                |`extension_degree*field_length` bytes    |First point's coordinates in the extension field         |
|rhs                |`group_order_length` bytes|Sсalar multiplication factor                                            |

### OPERATION_G1_MULTIEXP operands

The multiexponentiation operation can take arbitrary number of operands. Each of the operands must be encoded in the following form:

|Value              |Length                                   |                                                         |
|-------------------|-----------------------------------------|---------------------------------------------------------|
|point              |`extension_degree*field_length` bytes    |Point's coordinates in the extension field               |
|scalar             |`group_order_length` bytes|Sсalar order of exponentiation                                          |

## op_data for Pairing operations

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

Note that BLS12 is a family of curves that are parametrized by a single scalar `x`, twist type that is either `M` (multiplication) or `D` (division), and structure of the extension tower (non-residues). Nevertheless this ABI required caller to submit `base_field_modulus` and `main_subgroup_order` explicitly. It's also much more convenient for any observer to check validity of parameters for a given known BLS12 curve (e.g. `BLS12-381`). 

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

