# ABI interface

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

- NUM_LIMBS_MIN = 4;
- NUM_LIMBS_MAX = 16;
- NUM_GROUP_LIMBS_MIN = 1;
- NUM_GROUP_LIMBS_MAX = 16;
- MAX_MODULUS_BYTE_LEN = 128;
- MAX_GROUP_BYTE_LEN = 128;

## Sane limits

These limits are just simple upper bound to prevent long computations. Those do not in principle affect correctness, but allow to reject heavy computational calls early

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

Points of infinity are encoded as points with zero `X` and `Y` coordinates for both inputs and outputs. Precompile only works with curves in Weierstrass form with `b != 0` thus point `(0,0)` is not on curve.

## Encoding of elements in extension fields

For an element in extension field `c0 + c1*v + c2*v^3` where `v` is non-residue corresponding encoding is `(c0, c1, c2)`.

## Junk at the end of byte string

If after parsing the of all the parameters end of the byte string is not encountered then error must be returned.

## Representation of field elements

For reference implementations field elements were represented as a set of `N` 64-bit words where highest word (most significant) is underfilled (the highest bit of the highest word is unused) that allowed one to perform additions without extra overflow checks. Thus for `B` bits modulus one has to use `B/64 + 1` limbs (`/` is floor division) and such number of limbs is used for gas schedule that is described in the separate document.

## Precompile input (call data)

It is expected that various precompile addresses will be exposed for different particular operations. For any operation an an input byte string will be interpreted as described below in the corresponding sections.

All numbers are passed in **big endian** encoding.

Incorrect data input is always handled and returns an error.

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
- `a < base_field_modulus`
- `b < base_field_modulus`
- `b != 0`
- `group_order_length > 0`
- `group_order_length < MAX_GROUP_BYTE_LEN`
- top byte of `group_order` is non-zero (byte encoding is dense)
- `group_order != 0`

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
- all coordinates encodings are `<base_field_modulus`
- point is on curve
- `rhs <= group_order` (note `<=` check)

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
- all coordinates encodings are `<base_field_modulus`
- `num_pairs > 0`
- all points are on curve
- for each scalar: `scalar <= group_order` (note `<=` check)

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
  - check that `fp_non_residue` is indeed not a square (for extension degree 2) or cube (for degree 3) root by performing exponention:
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
|pairs              |`6*field_length*num_pairs`|Point pairs encoded as `(G1_point, G2_point)`|

Validations:
- All validations from G1 common prefix section
- `fp2_non_residue` is not a square root
- `fp6_non_residue` is not a 6-th root
- during computations of Frobenius endomorphism coefficients for all the field extensions (Fp2, Fp6 and Fp12) perform the following checks:
  - `base_field_modulus == 1 mod 2` 
  - `base_field_modulus == 1 mod 3` 
  - `base_field_modulus == 1 mod 6` 
- `x != 0`
- bit length of `x` is smaller or equal than `MAX_BLS12_X_BIT_LENGTH`
- hamming weight of `x` is smaller or equalt than `MAX_BLS12_X_HAMMING`
- `num_pairs > 0`
- all points are on the corresponding curves
- all points are in the claimed subgroups (!)
- filter out pairs where there are zero-points (so those do not contribute to result). If no points left return single byte `0x00`.  

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
|u                  |`x_length` bytes          |                                             |
|sign               |1 bytes                   |0 for plus, 1 for minus, sign of `u`         |
|num_pairs          |1 bytes                   |Number of point pairs                        |
|pairs              |`6*field_length*num_pairs`|Point pairs encoded as `(G1_point, G2_point)`|

Validations:
- All validations from G1 common prefix section
- `fp2_non_residue` is not a square root
- `fp6_non_residue` is not a 6-th root
- during computations of Frobenius endomorphism coefficients for all the field extensions (Fp2, Fp6 and Fp12) perform the following checks:
  - `base_field_modulus == 1 mod 2` 
  - `base_field_modulus == 1 mod 3` 
  - `base_field_modulus == 1 mod 6` 
- `u != 0`
- bit length of `u` is smaller or equal than `MAX_BN_U_BIT_LENGTH`
- hamming weight of `|6u + 2|` is smaller or equalt than `MAX_BN_SIX_U_PLUS_TWO_HAMMING`
- `num_pairs > 0`
- all points are on the corresponding curves
- all points are in the claimed subgroups (!)
- filter out pairs where there are zero-points (so those do not contribute to result). If no points left return single byte `0x00`.

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
|twist_type         |1 bytes                   |Can be either 0x01 for M or 0x02 for D       |
|loop_byte_length   |1 bytes                   |                                             |
|loop                   |`loop_byte_length` bytes          |                                             |
|loop_sign               |1 bytes                   |0 for plus, 1 for minus, sign of `loop`         |
|exp_w0_byte_length   |1 bytes                   |                                             |
|exp_w0                   |`exp_w0_byte_length` bytes          |                                             |
|exp_w1_byte_length   |1 bytes                   |                                             |
|exp_w1                   |`exp_w1_byte_length` bytes          |                                             |
|exp_w0_sign               |1 bytes                   |0 for plus, 1 for minus, sign of `exp_w0`         |
|num_pairs          |1 bytes                   |Number of point pairs                        |
|pairs              |`6*field_length*num_pairs`|Point pairs encoded as `(G1_point, G2_point)`|

Validations:
- All validations from G1 common prefix section
- `fp2_non_residue` is not a 4-th root
- during computations of Frobenius endomorphism coefficients for all the field extensions (Fp2 and Fp4) perform the following checks:
  - `base_field_modulus == 1 mod 2` 
  - `base_field_modulus == 1 mod 4` 
- `loop != 0`
- bit length of `loop` is smaller or equal than `MAX_ATE_PAIRING_ATE_LOOP_COUNT`
- hamming weight of `loop` is smaller or equalt than `MAX_ATE_PAIRING_ATE_LOOP_COUNT_HAMMING`
- `exp_w0 != 0`
- `exp_w1 != 0`
- hamming weight of `exp_w0` is smaller or equalt than `MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH`
- hamming weight of `exp_w1` is smaller or equalt than `MAX_ATE_PAIRING_FINAL_EXP_W1_BIT_LENGTH`
- `num_pairs > 0`
- all points are on the corresponding curves
- all points are in the claimed subgroups (!)
- filter out pairs where there are zero-points (so those do not contribute to result). If no points left return single byte `0x00`.

Return value:

If result of a pairing (element of `Fp12`) is equal to identity - return single byte `0x01`, otherwise return `0x00` following the existing ABI for BN254 precompile.

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
|twist_type         |1 bytes                   |Can be either 0x01 for M or 0x02 for D       |
|loop_byte_length   |1 bytes                   |                                             |
|loop                   |`loop_byte_length` bytes          |                                             |
|loop_sign               |1 bytes                   |0 for plus, 1 for minus, sign of `loop`         |
|exp_w0_byte_length   |1 bytes                   |                                             |
|exp_w0                   |`exp_w0_byte_length` bytes          |                                             |
|exp_w1_byte_length   |1 bytes                   |                                             |
|exp_w1                   |`exp_w1_byte_length` bytes          |                                             |
|exp_w0_sign               |1 bytes                   |0 for plus, 1 for minus, sign of `exp_w0`         |
|num_pairs          |1 bytes                   |Number of point pairs                        |
|pairs              |`8*field_length*num_pairs`|Point pairs encoded as `(G1_point, G2_point)`|

Validations:
- All validations from G1 common prefix section
- `fp2_non_residue` is not a 6-th root
- during computations of Frobenius endomorphism coefficients for all the field extensions (Fp2 and Fp4) perform the following checks:
  - `base_field_modulus == 1 mod 2` 
  - `base_field_modulus == 1 mod 3` 
- `loop != 0`
- bit length of `loop` is smaller or equal than `MAX_ATE_PAIRING_ATE_LOOP_COUNT`
- hamming weight of `loop` is smaller or equalt than `MAX_ATE_PAIRING_ATE_LOOP_COUNT_HAMMING`
- `exp_w0 != 0`
- `exp_w1 != 0`
- hamming weight of `exp_w0` is smaller or equalt than `MAX_ATE_PAIRING_FINAL_EXP_W0_BIT_LENGTH`
- hamming weight of `exp_w1` is smaller or equalt than `MAX_ATE_PAIRING_FINAL_EXP_W1_BIT_LENGTH`
- `num_pairs > 0`
- all points are on the corresponding curves
- all points are in the claimed subgroups (!)
- filter out pairs where there are zero-points (so those do not contribute to result). If no points left return single byte `0x00`.

Return value:

If result of a pairing (element of `Fp12`) is equal to identity - return single byte `0x01`, otherwise return `0x00` following the existing ABI for BN254 precompile.