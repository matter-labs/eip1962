# ABI interface

## Supported operations

EIP-1962 includes multiple elliptic curve operations. All operations will be performed on separate precompiled contracts. The mapping precompiles (listed on the left) to the addresses at which the precompile resides is defined as follows:

|Operation            |Code|
|---------------------|----|
|OPERATION_G1_ADD     |0x09|
|OPERATION_G1_MUL     |0x0A|
|OPERATION_G1_MULTIEXP|0x0B|
|OPERATION_G2_ADD     |0x0C|
|OPERATION_G2_MUL     |0x0D|
|OPERATION_G2_MULTIEXP|0x0E|
|OPERATION_PAIRING    |0x0F|

`OPERATION_G1_ADD`, `OPERATION_G1_MUL` and `OPERATION_G1_MULTIEXP` are operations of additon, multiplication and multiexponentiation for G1 elements of any curve in the Weierstrass form with `b != 0`.

`OPERATION_G2_ADD`, `OPERATION_G2_MUL` and `OPERATION_G2_MULTIEXP` are operations for G2 elements for a curve defined over some extension field.

`OPERATION_PAIRING` is the pairing operation. The following curve families are supported:

- BN
- BLS12
- MNT4
- MNT6

## Precompile input (call data)

Call data must be a correctly encoded ABI data string of two elements:

|Value  |Type       |Length  |
|-------|-----------|--------|
|op_data|bytes_array|variable|
|operands|bytes_array|variable|

The data is passed to the corresponding operation handler (see details below).

All numbers are passed in **big endian** encoding.

Incorrect data input is always handled and returns an error.

## op_data for G1 operations

`op_data` for all G1 operations consists of a common prefix followed by the operands.

The common prefix must have the following form:

|Value              |Length                    |Comment                    |
|-------------------|--------------------------|---------------------------|
|field_length       |1 byte                    |                           |
|base_field_modulus |`field_length` bytes      |Fq modulus                 |
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
|num_pairs          |1 byte                    | number of (point, scalar) pairs for multiexponentiation  
|point              |`2*field_length` bytes    |Point's X and Y coordinates       |
|scalar             |`group_order_length` bytes|Sсalar order of exponentiation    |


## op_data for G2 operations

`op_data` for all G2 operations consists of a common prefix followed by the operands.

The common prefix must have the following form:

|Value              |Length                    |Comment                       |
|-------------------|--------------------------|------------------------------|
|field_length       |1 byte                    |                              |
|base_field_modulus |`field_length` bytes      |Fq modulus                    |
|extension_degree   |1 bytes                   |Only values 2 or 3 are allowed|
|fp_non_residue     |`field_length` bytes      |Non-residue for Fp 2          |
|a                  |`field_length` bytes      |Curve's a coefficient         |
|b                  |`field_length` bytes      |Curve's b coefficient         |
|group_order_length |1 bytes                   |                              |                    
|group_order        |`group_order_length` bytes|Group order                   |

The operands are described below for each operation. They follow the same schema as for G1 operations, except that all points are encoded in the required extension degree.

### OPERATION_G2_ADD operands

|Value              |Length                                   |                                                          |
|-------------------|-----------------------------------------|----------------------------------------------------------|
|lhs                |`2*extension_degree*field_length` bytes    |First point's coordinates in the extension field          |
|rhs                |`2*extension_degree*field_length` bytes    |Second point's X and Y coordinates in the extension field |

### OPERATION_G2_MUL operands

|Value              |Length                                   |                                                         |
|-------------------|-----------------------------------------|---------------------------------------------------------|
|lhs                |`2*extension_degree*field_length` bytes    |First point's coordinates in the extension field         |
|rhs                |`group_order_length` bytes|Sсalar multiplication factor                                            |

### OPERATION_G2_MULTIEXP operands

The multiexponentiation operation can take arbitrary number of operands. Each of the operands must be encoded in the following form:

|Value              |Length                                   |                                                         |
|-------------------|-----------------------------------------|---------------------------------------------------------|
|num_pairs          |1 byte                    | number of (point, scalar) pairs for multiexponentiation 
|point              |`2*extension_degree*field_length` bytes    |Point's coordinates in the extension field               |
|scalar             |`group_order_length` bytes|Sсalar order of exponentiation                                          |

## op_data for Pairing operations

The first byte of `op_data` for every Pairing operation is the curve type, as defined below:

|Curve | Type |
|------|------|
|BLS12 | 0x01 |
|BN    | 0x02 |
|MNT4  | 0x03 |
|MNT6  | 0x04 |

### ABI for pairing operations on BLS12 curves

Note that BLS12 is a family of curves that are parametrized by a single scalar `x`, twist type that is either `M` (multiplication) or `D` (division), and structure of the extension tower (non-residues). Nevertheless this ABI required caller to submit `base_field_modulus` and `main_subgroup_order` explicitly. It's also much more convenient for any observer to check validity of parameters for a given known BLS12 curve (e.g. `BLS12-381`).

|Value              |Length                    |Comment                                      |
|-------------------|--------------------------|---------------------------------------------|
|curve_type         |1 byte                    |See table below                              |
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
|sign               |1 bytes                   |0 for plus, 1 for minus                      |
|num_pairs          |1 bytes                   |Number of point pairs                        |
|pairs              |`6*field_length*num_pairs`|Point pairs encoded as `(G1_point, G2_point)`|

Return value:

If result of a pairing (element of `Fp12`) is equal to identity - return single byte `0x01`, otherwise return `0x00` following the existing ABI for BN254 precompile.
