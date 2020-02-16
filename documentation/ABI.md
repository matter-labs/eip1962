# ABI interface

## Breaking changes history

- Encoding of group or order is now NOT required to be dense (so first byte is allowed to be zero)
- This also means that group order is calculated from full encoded byte length. In principle one can prepend 8 zero bytes and pay higher price for operations.
- Scalars for multiplication are now NOT required to be less or equal than the group order. This allows caller to have modular reduction "for free" and is already accounted in our pricing model
- There is now an optional byte BEFORE G1 or G2 point encoding in pairing calls indicating whether this point must be subgroup checked or not

## Supported operations

The precompile provides multiple elliptic curve operations. The full set of operations is defined as follows:

|Operation                  |Code|
|---------------------------|----|
|OPERATION_G1_ADD           |0x01|
|OPERATION_G1_MUL           |0x02|
|OPERATION_G1_MULTIEXP      |0x03|
|OPERATION_G2_ADD           |0x04|
|OPERATION_G2_MUL           |0x05|
|OPERATION_G2_MULTIEXP      |0x06|
|OPERATION_PAIRING_BLS12    |0x07|
|OPERATION_PAIRING_BN       |0x08|
|OPERATION_PAIRING_MNT4     |0x09|
|OPERATION_PAIRING_MNT6     |0x0a|

These operations perform internal addressing of what should be done with provided encoded input and do NOT correspond to the set of addresses that would be assigned to the precompile.

`OPERATION_G1_ADD`, `OPERATION_G1_MUL` and `OPERATION_G1_MULTIEXP` are operations of additon, multiplication and multiexponentiation for elements on any curve in the Weierstrass form with `b != 0` defined over base field.

`OPERATION_G2_ADD`, `OPERATION_G2_MUL` and `OPERATION_G2_MULTIEXP` are operations for elements on any curve in the Weierstrass form with `b != 0` defined over field extension of degree `2` or `3`

Following curve families are supported for pairing operations:
- BN
- BLS12
- MNT4
- MNT6

## Constants

These constants limit a scope of the curves that are supported by the precompile. E.g. `NUM_LIMBS_MAX` equal to `16` limits the number of bits of the modulus to `1023`. As `16` limbs of `64` bits can be encoded densely with maximum of `128` bytes it's also constrained here.

- NUM_LIMBS_MIN = 4;
- NUM_LIMBS_MAX = 16;
- NUM_GROUP_LIMBS_MIN = 1;
- NUM_GROUP_LIMBS_MAX = 16;
- MAX_MODULUS_BYTE_LEN = 128;
- MAX_GROUP_BYTE_LEN = 128;

## Sane limits

These limits are just simple upper bound to prevent long computations. Those do not in principle affect correctness, but allow to reject *unreasonably* heavy computational calls early.

NOTE: These limits may change

- MAX_BLS12_X_BIT_LENGTH = 128;
- MAX_BN_U_BIT_LENGTH = 128;
- MAX_BLS12_X_HAMMING = 128;
- MAX_BN_SIX_U_PLUS_TWO_HAMMING = 128;
- MAX_ATE_PAIRING_ATE_LOOP_COUNT = 2032;
- MAX_ATE_PAIRING_ATE_LOOP_COUNT_HAMMING = 2032;
- MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH = 2032;
- MAX_ATE_PAIRING_FINAL_EXP_W1_BIT_LENGTH = 2032;

## Zero point (point of infinity) encoding convension

Points of infinity are encoded as points with zero `X` and `Y` coordinates for both inputs and outputs. Precompile only works with curves in short Weierstrass form with `b != 0` thus point `(0,0)` is not on curve.

## Encoding of elements in extension fields

For an element in extension field `c0 + c1*v + c2*v^3` where `v` is non-residue corresponding encoding is `(c0, c1, c2)`.
This means that for e.g. when field modulus is encoded with 32 bytes (and length of the modulus defines length of the elements of the *base field* as explained also below) when one decodes an Fp2 element then in total one expects to encounter 64 bytes where first 32 bytes encode `c0` and next 32 bytes encode `c1`.

## Junk at the end of byte string

If after parsing of all the parameters the end of the byte string is not encountered then error must be returned.

## Internal representation of field elements

For reference implementations field elements were represented as a set of `N` 64-bit words (limbs) where highest word (most significant) is underfilled (the highest bit of the highest word is unused) that allowed one to perform additions without extra overflow checks. Thus for `B` bits modulus one has to use value of `modulus_libs = B/64 + 1` limbs (`/` is floor division) and such number of limbs is used for gas schedule that is described in the separate document.

## Encoding of the modulus

Encoding of the modulus (our main parameter) is explicitly unrolled in the specs below, but we present it here explicitly

- `field_length` is encoded as a single byte
- `field_length <= MAX_MODULUS_BYTE_LEN`
- `field_length > 0`
- next `field_length` bytes are interpreted as the BE encoding of the modulus `base_field_modulus`
- top byte of `base_field_modulus` is non-zero (byte encoding is dense)
- number of bits of `base_field_modulus` is less than `1024` (not less than or equal!). This is a limit for the maximum number of limbs supported to `16` as was described in the section above.

## Encoding of the `boolean` parameters

- boolean `false` is encoded as a single byte `0x00`
- boolean `true` is encoded as a single byte `0x01`

## Encoding of the `sign` parameters

- sign `+` is encoded as a single byte `0x00`
- sign `-` is encoded as a single byte `0x01`

## Encoding of the `twist type` parameters

- twist type `M` is encoded as a single byte `0x01`
- twist type `D` is encoded as a single byte `0x02`

## Encoding of the field elements and extension field elements

This encoding applied to e.g. `a` or `b` coefficients of the curve being used.

- length of the element is `field_length` (also explicitly unrolled in particular ABI structures later)
- encoded field element is strictly less than the modulus
- for field extensions (`Fp2`, `Fp3`) each of the coefficients is encoded properly

## Encoding of the group order

Encoding of this parameter is explicitly unrolled in the specs below, but we present it here explicitly.

- `group_order_length` is encoded as a single byte
- `group_order_length > 0`
- `group_order_length < MAX_GROUP_BYTE_LEN`
- next `group_order_length` bytes are interpreted as the BE encoding of the `group_order`
- `group_order != 0`

We mention again that `group_order` encoding is NOT dense.

## Encoding of loop parameters

Loop parameters are used in pairing operations and all have "sane" limits. Such limits also impose intrinsic restrictions on how they can be encoded. Let's take an example of the `x` parameter for `BLS12` curve to demonstrate all the checks. There `x` parameter has a limit `MAX_BLS12_X_BIT_LENGTH` for it's bit length and `MAX_BLS12_X_HAMMING` for it's hamming weight (hamming weight check is not a part of the *decoding routine*, but it's a part of *validation*).

- `x_length` is encoded as a single byte
- `x_length > 0`
- calculate `max_x_length` as `(MAX_BLS12_X_BIT_LENGTH + 7) / 8` (e.g. one can densely encode 16 bits with 2 bytes maximum)
- `x_length <= max_x_length`
- next `x_length` bytes are interpreted as the BE encoding of the `x`
- encoding of `x` is dense, same as for a modulus
- bit length of `x` is smaller or equal than `MAX_BLS12_X_BIT_LENGTH` (explicit check)

In this manner all the loop parameters can be decoded.

### Reason why group order encoding is not dense

Here we present an example why restriction of dense group order encoding was removed as well as the check that scalar value for multiplications is no longer required to be less than or equal to the group order.

In many cryptographic operations it's possible that a scalar parameter for multiplication or multiexponentiation call is generated roughly as the following:

- `scalar = Hash(...)` where `Hash` function returns a byte array of some length `hash_len`
- group order has length `group_order_length < hash_len` or `group_order < scalar`
- in the old version of the ABI caller would have to do modular reduction to bring scalar in range (make `scalar <= group order`). While such reduction is simple in Ethereum if `scalar` fits in 256 bits, it's not so trivial for larger scalars. Now this reduction is performed by the precompiled invocation itself and usually it's also free for a user: gas pricing only depended on how long the `group_order` in terms of 64 bit limbs, so if `group_order_length` was not divisible by 8 it still was rounded up even while arithmetic operations were not performed.

## Precompile input (call data)

It is expected that various precompile addresses will be exposed for different particular operations. For any operation an an input byte string will be interpreted as described below in the corresponding sections.

All numbers are passed in **big endian** encoding.

Incorrect data input is **always** handled and returns an error (gas estimation function **also** returns an error if input is malformed. It does NOT perform expensive arithmetic checks that involve e.g. divisions, but checks that input has proper length, encoding of specific elements follows the conventions, etc).

Validation steps for properly encoded and sane data are described along with the data encoding format. These steps are largely an input validation and happen before main arithmetic is performed. 

## Shared input data for all G1 operations

Input data for all G1 operations consists of a common prefix followed by the operands.

The common prefix must have the following form:

|Value              |Length                    |Comment                    |
|-------------------|--------------------------|---------------------------|
|field_length       |1 byte                    |                           |
|base_field_modulus |`field_length` bytes      |Field modulus              |
|a                  |`field_length` bytes      |Curve's a coefficient      |
|b                  |`field_length` bytes      |Curve's b coefficient      |
|group_order_length |1 byte                    |                           |                    
|group_order        |`group_order_length` bytes|Group order                |

The operands are described below for each operation.

The following list of validation steps is performed during parsing:
- `field length <= MAX_MODULUS_BYTE_LEN`
- `field length > 0`
- top byte of `base_field_modulus` is non-zero (byte encoding is dense)
- number of bits of `base_field_modulus` is less than `1024` (not less than or equal!)
- `base_field_modulus` is odd
- `base_field_modulus > 3`
- `a` is properly encoded (*not performed during gas estimation*)
- `b` is properly encoded (*not performed during gas estimation*)
- `b != 0` (*not performed during gas estimation*)
- `group_order_length > 0`
- `group_order_length < MAX_GROUP_BYTE_LEN`
- ~~top byte of `group_order` is non-zero (byte encoding is dense)~~
- `group_order != 0` (*not performed during gas estimation*)

### OPERATION_G1_ADD operands

|Value              |Length                    |                                  |
|-------------------|--------------------------|----------------------------------|
|lhs                |`2*field_length` bytes    |First point's X and Y coordinates |
|rhs                |`2*field_length` bytes    |Second point's X and Y coordinates|

Validations:
- all coordinates encodings are `<base_field_modulus`
- points are on curve

Return value:

`2*field_length` bytes - encoded X and Y coordinates of the result point

### OPERATION_G1_MUL operands

|Value              |Length                    |                                  |
|-------------------|--------------------------|----------------------------------| 
|lhs                |`2*field_length` bytes    |First point's X and Y coordinates |
|rhs                |`group_order_length` bytes|Sсalar multiplication factor      |

Validations:
- all coordinates encodings are valid (*not performed during gas estimation*)
- point is on curve (*not performed during gas estimation*)
- ~~`rhs <= group_order` (note `<=` check)~~

Return value:

`2*field_length` bytes - encoded X and Y coordinates of the result point

### OPERATION_G1_MULTIEXP operands

The multiexponentiation operation can take arbitrary number of operands. Each of the operands must be encoded in the following form:

|Value              |Length                    |                                  |
|-------------------|--------------------------|----------------------------------|
|num_pairs          |1 byte                    | number of (point, scalar) pairs for multiexponentiation  |
|-------------------|--------------------------|----------------------------------|
|point              |`2*field_length` bytes    |Point's X and Y coordinates       |
|scalar             |`group_order_length` bytes|Sсalar order of exponentiation    |
|-------------------|--------------------------|----------------------------------|
|...           |...|...|
|-------------------|--------------------------|----------------------------------|
|point              |`2*field_length` bytes    |Point's X and Y coordinates       |
|scalar             |`group_order_length` bytes|Sсalar order of exponentiation    |
|-------------------|--------------------------|----------------------------------|


here sequence of `(point, scalar)` are must repeat `num_pairs` times

Validations:
- all coordinates encodings are valid (*not performed during gas estimation*)
- `num_pairs > 0`
- all points are on curve (*not performed during gas estimation*)
- ~~for each scalar: `scalar <= group_order` (note `<=` check)~~

Return value:

`2*field_length` bytes - encoded X and Y coordinates of the result point


## op_data for G2 operations

Input data for all G2 operations consists of a common prefix followed by the operands.

The common prefix must have the following form:

|Value              |Length                    |Comment                       |
|-------------------|--------------------------|------------------------------|
|field_length       |1 byte                    |                              |
|base_field_modulus |`field_length` bytes      |Field modulus                 |
|extension_degree   |1 byte                    |Only values 2 or 3 are allowed|
|fp_non_residue     |`field_length` bytes      |Non-residue for Fp2 or Fp3           |
|a                  |`field_length*extension_degree` bytes      |Curve's a coefficient (in Fp2 or Fp3)        |
|b                  |`field_length*extension_degree` bytes      |Curve's b coefficient (in Fp2 or Fp3)        |
|group_order_length |1 byte                    |                              |                    
|group_order        |`group_order_length` bytes|Group order                   |

The operands are described below for each operation. They follow the same schema as for G1 operations, except that all points are encoded in the required extension degree.

Validations:
- All the validations from G1 operations
- Additionally:
  - check that `fp_non_residue` is indeed not a square (for extension degree 2) or cube (for degree 3) root by performing exponention (*not performed during gas estimation*):
    - For extension degree 2:  `fp_non_residue ^ ((base_field_modulus - 1) / 2) != 1`
    - For extension degree 3:  `fp_non_residue ^ ((base_field_modulus - 1) / 3) != 1`
    - For both cases check that division `(base_field_modulus - 1) / 2)` or `(base_field_modulus - 1) / 3)` is exact (otherwise `base_field_modulus == 1 mod 2` or `base_field_modulus == 1 mod 3` respectively)

### OPERATION_G2_ADD operands

|Value              |Length                                   |                                                          |
|-------------------|-----------------------------------------|----------------------------------------------------------|
|lhs                |`extension_degree*field_length` bytes    |First point's coordinates in the extension field          |
|rhs                |`extension_degree*field_length` bytes    |Second point's X and Y coordinates in the extension field |

Validations:
- All the validations from G1 operations

Return value:

`2*field_length*extension_degree` bytes - encoded X and Y coordinates of the result point

### OPERATION_G2_MUL operands

|Value              |Length                                   |                                                         |
|-------------------|-----------------------------------------|---------------------------------------------------------|
|lhs                |`extension_degree*field_length` bytes    |First point's coordinates in the extension field         |
|rhs                |`group_order_length` bytes|Sсalar multiplication factor                                            |

Validations:
- All the validations from G1 operations

Return value:

`2*field_length*extension_degree` bytes - encoded X and Y coordinates of the result point

### OPERATION_G2_MULTIEXP operands

The multiexponentiation operation can take arbitrary number of operands. Each of the operands must be encoded in the following form:

|Value              |Length                                   |                                                         |
|-------------------|-----------------------------------------|---------------------------------------------------------|
|num_pairs          |1 byte                    | number of (point, scalar) pairs for multiexponentiation |
|-------------------|--------------------------|----------------------------------|
|point              |`extension_degree*field_length` bytes    |Point's coordinates in the extension field               |
|scalar             |`group_order_length` bytes|Sсalar order of exponentiation                                          |
|-------------------|--------------------------|----------------------------------|
|...           |...|...|
|-------------------|--------------------------|----------------------------------|
|point              |`extension_degree*field_length` bytes    |Point's coordinates in the extension field               |
|scalar             |`group_order_length` bytes|Sсalar order of exponentiation                                          |
|-------------------|--------------------------|----------------------------------|

Validations:
- All the validations from G1 operations

Return value:

`2*field_length*extension_degree` bytes - encoded X and Y coordinates of the result point

## Pairing operations

Pairing operations require much more steps in validation that is performed during parsing, as well as for different curve types ABI formats differ a lot.

### ABI for pairing operations on BLS12 curves

Note that BLS12 is a family of curves that are parametrized by a single scalar `x`, twist type that is either `M` (multiplication) or `D` (division), and structure of the extension tower (non-residues). Nevertheless this ABI required caller to submit `base_field_modulus` and `main_subgroup_order` explicitly. It's also much more convenient for any observer to check validity of parameters for a given known BLS12 curve (e.g. `BLS12-381`).

|Value              |Length                    |Comment                                      |
|-------------------|--------------------------|---------------------------------------------|
|field_length       |1 byte                    |                                             |
|base_field_modulus |`field_length` bytes      |Fq modulus                                   |
|a                  |`field_length` bytes      |Curve's a coefficient                        |
|b                  |`field_length` bytes      |Curve's b coefficient                        |
|group_order_length |1 bytes                   |                                             |                 
|main_subgroup_order|`group_order_length` bytes|Main subgroup order                          |
|fp2_non_residue    |`field_length` bytes      |Non-residue for Fp 2                         |
|fp6_non_residue    |`2*field_length` bytes    |Non-residue for Fp 6                         |
|twist_type         |1 bytes                   |Can be either 0x01 for M or 0x02 for D       |
|x_length           |1 bytes                   |                                             |
|x                  |`x_length` bytes          |                                             |
|sign               |1 bytes                   |0 for plus, 1 for minus, sign of `x`         |
|num_pairs          |1 bytes                   |Number of point pairs                        |
|pairs              |`2 + 6*field_length*num_pairs`|Point pairs encoded as `(check_g1_boolean, G1_point, check_g2_boolean, G2_point)`|

Validations:
- All validations from G1 common prefix section
- `fp2_non_residue` is not a square root (*not performed during gas estimation*)
- `fp6_non_residue` is not a 6-th root (*not performed during gas estimation*)
- during computations of Frobenius endomorphism coefficients for all the field extensions (Fp2, Fp6 and Fp12) perform the following checks (*not performed during gas estimation*):
  - `base_field_modulus == 1 mod 2` 
  - `base_field_modulus == 1 mod 3` 
  - `base_field_modulus == 1 mod 6` 
- `x_length` > 0
- `x != 0`
- encoding of `x` is dense(!)
- bit length of `x` is smaller or equal than `MAX_BLS12_X_BIT_LENGTH`
- hamming weight of `x` is smaller or equalt than `MAX_BLS12_X_HAMMING`
- `num_pairs > 0`
- all points are on the corresponding curves (*not performed during gas estimation*)
- ~~all points are in the claimed subgroups (!)~~
- for G1 or G2 points where the corresponding `check_g1_boolean` or `check_g2_boolean` is `true` points are checked to be in the correct subgroup (*not performed during gas estimation*)
- calculate a total number of `check_g1_boolean == true` and `check_g2_boolean == true` into the separate variables `num_g1_checks` and `num_g2_checks` (used for gas estimation only)
- filter out pairs where there are zero-points (so those do not contribute to result). If no points left return single byte `0x01`.  

Return value:

If result of a pairing (element of `Fp12`) is equal to identity - return single byte `0x01`, otherwise return `0x00` following the existing ABI for BN254 precompile.

### ABI for pairing operations on BN curves

|Value              |Length                    |Comment                                      |
|-------------------|--------------------------|---------------------------------------------|
|field_length       |1 byte                    |                                             |
|base_field_modulus |`field_length` bytes      |Fq modulus                                   |
|a                  |`field_length` bytes      |Curve's a coefficient                        |
|b                  |`field_length` bytes      |Curve's b coefficient                        |
|group_order_length |1 bytes                   |                                             |                 
|main_subgroup_order|`group_order_length` bytes|Main subgroup order                          |
|fp2_non_residue    |`field_length` bytes      |Non-residue for Fp 2                         |
|fp6_non_residue    |`2*field_length` bytes    |Non-residue for Fp 6                         |
|twist_type         |1 bytes                   |Can be either 0x01 for M or 0x02 for D       |
|u_length           |1 bytes                   |                                             |
|u                  |`u_length` bytes          |                                             |
|sign               |1 bytes                   |0 for plus, 1 for minus, sign of `u`         |
|num_pairs          |1 bytes                   |Number of point pairs                        |
|pairs              |`2 + 6*field_length*num_pairs`|Point pairs encoded as `(check_g1_boolean, G1_point, check_g2_boolean, G2_point)`|

Validations:
- All validations from G1 common prefix section
- `fp2_non_residue` is not a square root (*not performed during gas estimation*)
- `fp6_non_residue` is not a 6-th root (*not performed during gas estimation*)
- during computations of Frobenius endomorphism coefficients for all the field extensions (Fp2, Fp6 and Fp12) perform the following checks (*not performed during gas estimation*):
  - `base_field_modulus == 1 mod 2` 
  - `base_field_modulus == 1 mod 3` 
  - `base_field_modulus == 1 mod 6` 
- `u_length` > 0
- `u != 0`
- encoding of `u` is dense(!)
- bit length of `u` is smaller or equal than `MAX_BN_U_BIT_LENGTH`
- hamming weight of `|6u + 2|` is smaller or equal than `MAX_BN_SIX_U_PLUS_TWO_HAMMING`
- `num_pairs > 0`
- all points are on the corresponding curves (*not performed during gas estimation*)
- ~~all points are in the claimed subgroups (!)~~
- for G1 or G2 points where the corresponding `check_g1_boolean` or `check_g2_boolean` is `true` points are checked to be in the correct subgroup (*not performed during gas estimation*)
- calculate a total number of `check_g1_boolean == true` and `check_g2_boolean == true` into the separate variables `num_g1_checks` and `num_g2_checks` (used for gas estimation only)
- filter out pairs where there are zero-points (so those do not contribute to result). If no points left return single byte `0x01`.

Return value:

If result of a pairing (element of `Fp12`) is equal to identity - return single byte `0x01`, otherwise return `0x00` following the existing ABI for BN254 precompile.

### ABI for pairing operations on MNT4 curves

|Value              |Length                    |Comment                                      |
|-------------------|--------------------------|---------------------------------------------|
|field_length       |1 byte                    |                                             |
|base_field_modulus |`field_length` bytes      |Fq modulus                                   |
|a                  |`field_length` bytes      |Curve's a coefficient                        |
|b                  |`field_length` bytes      |Curve's b coefficient                        |
|group_order_length |1 bytes                   |                                             |                 
|main_subgroup_order|`group_order_length` bytes|Main subgroup order                          |
|fp2_non_residue    |`field_length` bytes      |Non-residue for Fp 2                         |
|loop_byte_length   |1 bytes                   |                                             |
|ate_loop_parameter                   |`loop_byte_length` bytes          |                                             |
|ate_loop_sign               |1 bytes                   |0 for plus, 1 for minus, sign of `ate_loop_parameter`         |
|exp_w0_byte_length   |1 bytes                   |                                             |
|exp_w0                   |`exp_w0_byte_length` bytes          |                                             |
|exp_w1_byte_length   |1 bytes                   |                                             |
|exp_w1                   |`exp_w1_byte_length` bytes          |                                             |
|exp_w0_sign               |1 bytes                   |0 for plus, 1 for minus, sign of `exp_w0`         |
|num_pairs          |1 bytes                   |Number of point pairs                        |
|pairs              |`2 + 6*field_length*num_pairs`|Point pairs encoded as `(check_g1_boolean, G1_point, check_g2_boolean, G2_point)`|

Validations:
- All validations from G1 common prefix section
- `fp2_non_residue` is not a 4-th root (*not performed during gas estimation*)
- during computations of Frobenius endomorphism coefficients for all the field extensions (Fp2 and Fp4) perform the following checks (*not performed during gas estimation*):
  - `base_field_modulus == 1 mod 2` 
  - `base_field_modulus == 1 mod 4` 
- `loop_byte_length > 0`
- `ate_loop_parameter != 0`
- encoding of `ate_loop_parameter` is dense(!)
- bit length of `ate_loop_parameter` is smaller or equal than `MAX_ATE_PAIRING_ATE_LOOP_COUNT`
- hamming weight of `ate_loop_parameter` is smaller or equalt than `MAX_ATE_PAIRING_ATE_LOOP_COUNT_HAMMING`
- `exp_w0_byte_length > 0`
- `exp_w0 != 0`
- encoding of `exp_w0` is dense(!)
- `exp_w1_byte_length > 0`
- `exp_w1 != 0`
- encoding of `exp_w1` is dense(!)
- hamming weight of `exp_w0` is smaller or equalt than `MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH`
- hamming weight of `exp_w1` is smaller or equalt than `MAX_ATE_PAIRING_FINAL_EXP_W1_BIT_LENGTH`
- `num_pairs > 0`
- all points are on the corresponding curves (*not performed during gas estimation*)
- ~~all points are in the claimed subgroups (!)~~
- for G1 or G2 points where the corresponding `check_g1_boolean` or `check_g2_boolean` is `true` points are checked to be in the correct subgroup (*not performed during gas estimation*)
- calculate a total number of `check_g1_boolean == true` and `check_g2_boolean == true` into the separate variables `num_g1_checks` and `num_g2_checks` (used for gas estimation only)
- filter out pairs where there are zero-points (so those do not contribute to result). If no points left return single byte `0x01`.

Return value:

If result of a pairing (element of `Fp4`) is equal to identity - return single byte `0x01`, otherwise return `0x00` following the existing ABI for BN254 precompile.

### ABI for pairing operations on MNT6 curves

|Value              |Length                    |Comment                                      |
|-------------------|--------------------------|---------------------------------------------|
|field_length       |1 byte                    |                                             |
|base_field_modulus |`field_length` bytes      |Fq modulus                                   |
|a                  |`field_length` bytes      |Curve's a coefficient                        |
|b                  |`field_length` bytes      |Curve's b coefficient                        |
|group_order_length |1 bytes                   |                                             |                 
|main_subgroup_order|`group_order_length` bytes|Main subgroup order                          |
|fp2_non_residue    |`field_length` bytes      |Non-residue for Fp 3                         |
|loop_byte_length   |1 bytes                   |                                             |
|ate_loop_parameter                   |`loop_byte_length` bytes          |                                             |
|ate_loop_sign               |1 bytes                   |0 for plus, 1 for minus, sign of `ate_loop_parameter`         |
|exp_w0_byte_length   |1 bytes                   |                                             |
|exp_w0                   |`exp_w0_byte_length` bytes          |                                             |
|exp_w1_byte_length   |1 bytes                   |                                             |
|exp_w1                   |`exp_w1_byte_length` bytes          |                                             |
|exp_w0_sign               |1 bytes                   |0 for plus, 1 for minus, sign of `exp_w0`         |
|num_pairs          |1 bytes                   |Number of point pairs                        |
|pairs              |`2 + 8*field_length*num_pairs`|Point pairs encoded as `(check_g1_boolean, G1_point, check_g2_boolean, G2_point)`|

Validations:
- All validations from G1 common prefix section
- `fp2_non_residue` is not a 6-th root (*not performed during gas estimation*)
- during computations of Frobenius endomorphism coefficients for all the field extensions (Fp2 and Fp4) perform the following checks (*not performed during gas estimation*):
  - `base_field_modulus == 1 mod 2` 
  - `base_field_modulus == 1 mod 3` 
- `loop_byte_length > 0`
- `ate_loop_parameter != 0`
- encoding of `ate_loop_parameter` is dense(!)
- bit length of `ate_loop_parameter` is smaller or equal than `MAX_ATE_PAIRING_ATE_LOOP_COUNT`
- hamming weight of `ate_loop_parameter` is smaller or equalt than `MAX_ATE_PAIRING_ATE_LOOP_COUNT_HAMMING`
- `exp_w0_byte_length > 0`
- `exp_w0 != 0`
- encoding of `exp_w0` is dense(!)
- `exp_w1_byte_length > 0`
- `exp_w1 != 0`
- encoding of `exp_w1` is dense(!)
- hamming weight of `exp_w0` is smaller or equalt than `MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH`
- hamming weight of `exp_w1` is smaller or equalt than `MAX_ATE_PAIRING_FINAL_EXP_W1_BIT_LENGTH`
- `num_pairs > 0`
- all points are on the corresponding curves (*not performed during gas estimation*)
- ~~all points are in the claimed subgroups (!)~~
- for G1 or G2 points where the corresponding `check_g1_boolean` or `check_g2_boolean` is `true` points are checked to be in the correct subgroup (*not performed during gas estimation*)
- calculate a total number of `check_g1_boolean == true` and `check_g2_boolean == true` into the separate variables `num_g1_checks` and `num_g2_checks` (used for gas estimation only)
- filter out pairs where there are zero-points (so those do not contribute to result). If no points left return single byte `0x01`.

Return value:

If result of a pairing (element of `Fp6`) is equal to identity - return single byte `0x01`, otherwise return `0x00` following the existing ABI for BN254 precompile.


## Example of the input parsing

The following byte string (hex encoded) represents a call data to the BLS12 pairing function to perform a pairing for one pair of points:

```
301a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000042073eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff000000011a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaaa0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010108d20100000001000001010117f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e101024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb813e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b828010606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be
```

Let's parse it:
- `30` single byte - `field_length`, 48 bytes (`0x30` in hex)
- `1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab` - 48 bytes that encode the modulus in BE encoding. This is a prime field of BLS12-381 curve
- `000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000` - 48 bytes that encode the `a` coefficient for the short Weierstrass curve. Indeed the BLS12-381 curve has `a = 0`
- `000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000004` - 48 bytes that encode the `b` coefficient for the short Weierstrass curve. Indeed the BLS12-381 curve has `b = 4`
- `20` single byte - `group_order_length`, 32 bytes
- `73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001` - 32 bytes that encode the main subgroup order
- `1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaaa` - 48 bytes that encode an `Fp` element that is quadratic non-residue and is used to construct `Fp2` extension. It's `-1` (`modulus - 1`)
- Fp2 element that is not 6th root and forms `Fp6` and `Fp12` extensions. It is encoded as two `Fp` elemements
  - `000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001` - 48 bytes
  - `000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001` - 48 bytes
  - Fp2 element is `1 + 1*v` for BLS12-381
- `01` single byte - encodes twist type. For BLS12-381 twist type is `M` - multiplication.
- `08` single byte - encodes `x_length`
- `d201000000010000` - 8 bytes, modulus of `x` parameter of BLS12 curves
- `01` single byte - encodes that `x` parameter is negative
- `01` single byte - encodes `num_pairs`. We call it for pairing of one pair
  - `01` single byte - encodes `check_g1_boolean = true` and we want to perform a subgroup check for the G1 point
  - G1 point
    - `17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb` - 48 bytes, x coordiante
    - `08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1` - 48 bytes, y coordiante
  - `01` single byte - encodes `check_g2_boolean = true` and we want to perform a subgroup check for the G2 point
  - G2 point
    - X coordinate in Fp2
      - `024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8` - 48 bytes, first coefficient
      - `13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e` - 48 bytes, second coefficient
    - Y coordinate in Fp2
      - `0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801` - 48 bytes, first coefficient
      - `0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be` - 48 bytes, second coefficient

