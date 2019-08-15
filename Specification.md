# Specification

Universality of the precompile requires to think about edge cases. Proper splitting of the structure ensures that if initial parameters (such as field modulus non-residues for extension, etc.) are passed correctly, then all the remaining arithmetic is well-defined. This document can be viewed as API spec, implementation guide, and user guide for some aspects

## Supported operations ABI

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

`OPERATION_G2_ADD`, `OPERATION_G2_MUL` and `OPERATION_G2_MULTIEXP` are operations for G2 elements for a curve defined over some extension field. There are only two such extensions supported: degree 2 and degree 3.

`OPERATION_PAIRING` is the pairing operation. The following curve families are supported:

- BN
- BLS12
- MNT4
- MNT6

### Precompile input (call data)

Call data must be a correctly encoded ABI data string of two elements:

|Value  |Type       |Length  |
|-------|-----------|--------|
|op_data|bytes_array|variable|
|operands|bytes_array|variable|

The data is passed to the corresponding operation handler (see details below).

All numbers are passed in **big endian** encoding.

Incorrect data input is always handled and returns an error.

### op_data for G1 operations

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

#### OPERATION_G1_ADD operands

|Value              |Length                    |                                  |
|-------------------|--------------------------|----------------------------------|
|lhs                |`2*field_length` bytes    |First point's X and Y coordinates |
|rhs                |`2*field_length` bytes    |Second point's X and Y coordinates|

#### OPERATION_G1_MUL operands

|Value              |Length                    |                                  |
|-------------------|--------------------------|----------------------------------| 
|lhs                |`2*field_length` bytes    |First point's X and Y coordinates |
|rhs                |`group_order_length` bytes|Sсalar multiplication factor      |

#### OPERATION_G1_MULTIEXP operands

The multiexponentiation operation can take arbitrary number of operands. Each of the operands must be encoded in the following form:

|Value              |Length                    |                                  |
|-------------------|--------------------------|----------------------------------|
|num_pairs          |1 byte                    | number of (point, scalar) pairs for multiexponentiation  
|point              |`2*field_length` bytes    |Point's X and Y coordinates       |
|scalar             |`group_order_length` bytes|Sсalar order of exponentiation    |


### op_data for G2 operations

`op_data` for all G2 operations consists of a common prefix followed by the operands.

The common prefix must have the following form:

|Value              |Length                    |Comment                       |
|-------------------|--------------------------|------------------------------|
|field_length       |1 byte                    |                              |
|base_field_modulus |`field_length` bytes      |Fq modulus                    |
|extension_degree   |1 bytes                   |Only values 2 or 3 are allowed|
|fp_non_residue     |`field_length` bytes      |Non-residue for Fp 2          |
|a                  |`extension_degree*field_length` bytes      |Curve's a coefficient         |
|b                  |`extension_degree*field_length` bytes      |Curve's b coefficient         |
|group_order_length |1 bytes                   |                              |                    
|group_order        |`group_order_length` bytes|Group order                   |

The operands are described below for each operation. They follow the same schema as for G1 operations, except that all points are encoded in the required extension degree.

#### OPERATION_G2_ADD operands

|Value              |Length                                   |                                                          |
|-------------------|-----------------------------------------|----------------------------------------------------------|
|lhs                |`2*extension_degree*field_length` bytes    |First point's coordinates in the extension field          |
|rhs                |`2*extension_degree*field_length` bytes    |Second point's X and Y coordinates in the extension field |

#### OPERATION_G2_MUL operands

|Value              |Length                                   |                                                         |
|-------------------|-----------------------------------------|---------------------------------------------------------|
|lhs                |`2*extension_degree*field_length` bytes    |First point's coordinates in the extension field         |
|rhs                |`group_order_length` bytes|Sсalar multiplication factor                                            |

#### OPERATION_G2_MULTIEXP operands

The multiexponentiation operation can take arbitrary number of operands. Each of the operands must be encoded in the following form:

|Value              |Length                                   |                                                         |
|-------------------|-----------------------------------------|---------------------------------------------------------|
|num_pairs          |1 byte                    | number of (point, scalar) pairs for multiexponentiation 
|point              |`2*extension_degree*field_length` bytes    |Point's coordinates in the extension field               |
|scalar             |`group_order_length` bytes|Sсalar order of exponentiation                                          |

### op_data for Pairing operations

The first byte of `op_data` for every Pairing operation is the curve type, as defined below:

|Curve | Type |
|------|------|
|BLS12 | 0x01 |
|BN    | 0x02 |
|MNT4  | 0x03 |
|MNT6  | 0x04 |

#### ABI for pairing operations on BLS12 curves

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

## Implementation ABI and parsing checks

The input is passed to the corresponding function of the operation called.

G1 additions, multiplications and multiexponentiations are defined for any curve in the Weierstrass form with `b != 0`. Operations in G2 are performed over the curve defined over some extension field. There are only two such extensions supported: degree 2 and degree 3.

Important constants:

|Name                            |Value |
|--------------------------------|------|
|EXTENSION_DEGREE_ENCODING_LENGTH|  1   |
|EXTENSION_DEGREE_2              | 0x02 |
|EXTENSION_DEGREE_3              | 0x03 |

Pairing operation is defined only for the following families of curves:
- BN
- BLS12
- MNT4
- MNT6
- Ate pairing for a generic curve in the Weierstrass form with `k = 6` and extension tower `Fp - Fp3 - Fp6` (sibling of MNT6)
- Ate pairing for a generic curve in the Weierstrass form with `k = 4` and extension tower `Fp - Fp2 - Fp4` (sibling of MNT4)

### General notice about ABI

It's mandatory that after all required parameters for a curve were read by the parser, there should be no data left to read -> no garbage at the end of input!

### ABI and parsing checks for G1 operations

ABI consists of two parts: one defines a base field and a curve, with another is operation dependent and encodes points or scalars for a corresponding operations. 

Signature of public function is just `Operation(byte_array)` (one public function), so just a pointer to the array of bytes is passes to the function.

#### Common ABI and parsing checks for G1 operations

Important! `take(N)` operation comsumes(!) first `N` bytes from the byte array

Algorithm:
- ensure that length of the byte array is `> BYTES_FOR_LENGTH_ENCODING`, otherwise return error
- `take(BYTES_FOR_LENGTH_ENCODING)` and parse it as unsigned integer `modulus_length` for encoding of length of the next value. Limit `BYTES_FOR_LENGTH_ENCODING = 1` ensures that byte length of the next parameter that is modulus of the prime field is bounded (not more than 255 bytes), so it's a first sanity check
- ensure that `modulus_length > 0`, otherwise return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is a modulus of a base prime field `base_field_modulus`
- ensure that top (most significant) byte of `base_field_modulus` is non-zero. This is not an attack and will only require caller to pay more gas, but it's a trivial check
- ensure that `base_field_modulus >= 3` and `base_field_modulus` is odd. This is the second sanity check and also guarantees that Montgommery form that is used for all the field elements is well-defined (is requires `gcd(modulus, R) == 1` with `R` being power of two in our cases). There is no primarity testing, but arithmetic operations now will not trigger panics
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `A` coefficient for a curve in the Weierstrass form
- ensure that `A < base_field_modulus`, otherwise return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `B` coefficient for a curve in the Weierstrass form
- ensure that `B < base_field_modulus`, otherwise return error
- ensure that `B > 0`, otherwise return error
- ensure that length of the byte array is `> BYTES_FOR_LENGTH_ENCODING`, otherwise return error
- `take(BYTES_FOR_LENGTH_ENCODING)` and parse it as unsigned integer `subgroup_order_length` for encoding of length of the next value
- ensure that `subgroup_order_length > 0`, otherwise return error
- ensure that length of the byte array is `> subgroup_order_length`, otherwise return error
- `take(subgroup_order_length)` and parse it as BigEndian encoding of an unsigned integer that is a main subgroup order `main_subgroup_order`
- ensure that `main_subgroup_order > 0`
- there are two purposes to require `main_subgroup_order` to be included into the ABI:
  - Upfront estimation of the worst-case scenario for a difficulty of the multiplication and multiexponentiation operations
  - One should not enforce point being in a correct subgroup for correctness of operations and if caller expects to process user's inputs he should make such check as a separate call. Nevertheless, this check is MANDATORY for pairing operations
- later should follow operation-specific parameters, depending on curve
  
This list is a naive algorithm. Real implementation will merge e.g. reading of `A`, range check and parsing it as a representation of some form of the field element in one operation.

Arithmetic is made using some fixed-length representation of a field element. This implementation follows an approach to represent them as a fixed length array `[u64; NUM_LIMBS`], such that for a modulus `M`: `2^(MODULUS_LIMBS*64) > 2*M` to ensure that one never has to take care about carry flags. In this case a field element with `255` bit modulus would be represented as `[u64; 4]`, but `256` bit modulus will be already represented as `[u64; 5]`

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

- `take(BYTES_FOR_LENGTH_ENCODING)` to get `expected_number_pairs`
- calculate expected byte length of one `(point, scalar)` pair as `expected_pair_len = 2*modulus_length + subgroup_order_length`
- ensure that length of the byte array is `== expected_number_pairs*expected_pair_len`, otherwise return error
- in a loop:
  - `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `x` of the point `P0`
  - ensure that `P0.x < base_field_modulus`, otherwise return error
  - `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `y` of the point `P0`
  - ensure that `P0.y < base_field_modulus`, otherwise return error
  - `take(subgroup_order_length)` and parse it as BigEndian encoding of an unsigned integer that is an `scalar` for multiplication operation
  - ensure that `scalar <= main_subgroup_order`, otherwise return error
- perform a miltiexponentiation depending on a form of the Weierstrass curve (`A == 0` or `A != 0`, `B` is always non-zero), output `P_res`
- operations are most likely to be performed in Jacobial coordinates, so perform normalization into affine coordinates. This require to make inversion of `P_res.z`. If `P_res.z == 0` return point of infinity, otherwise inverse `P_res.z` and perform normalization

#### Notes for G1 ABI:

For purposes of caller's convenience it may be reasonable to have some granularity for lengths of encodings of various element. For example, make encodings multiples of 8 bytes (64 bits) to roughtly correspond to limb bits on x64 machine

### ABI and parsing checks for G2 operations

ABI consists of two parts: one defines a base field and a curve, with another is operation dependent and encodes points or scalars for a corresponding operations. 

Signature of public function is just `Operation(byte_array)` (one public function), so just a pointer to the array of bytes is passes to the function.

#### Common ABI and parsing checks for G2 operations

Operations on a "twist" are defined and expected to be used for pairing friendly curves. E.G. original protocol of BLS aggregated signatures requires multiplication in G2, as well as some SNARK verification equations.

To save space only common ABI part for G2 is described, with specific part being similar to G1 part.

Important! `take(N)` operation comsumes(!) first `N` bytes from the byte array

Algorithm:
- ensure that length of the byte array is `> BYTES_FOR_LENGTH_ENCODING`, otherwise return error
- `take(BYTES_FOR_LENGTH_ENCODING)` and parse it as unsigned integer `modulus_length` for encoding of length of the next value. Limit `BYTES_FOR_LENGTH_ENCODING = 1` ensures that byte length of the next parameter that is modulus of the prime field is bounded (not more than 255 bytes), so it's a first sanity check
- ensure that `modulus_length > 0`, otherwise return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is a modulus of a base prime field `base_field_modulus` 
- ensure that top (most significant) byte of `base_field_modulus` is non-zero. This is not an attack and will only require caller to pay more gas, but it's a trivial check
- ensure that `base_field_modulus >= 3` and `base_field_modulus` is odd. This is the second sanity check and also guarantees that Montgommery form that is used for all the field elements is well-defined (is requires `gcd(modulus, R) == 1` with `R` being power of two in our cases). There is no primarity testing, but arithmetic operations now will not trigger panics
- ensure that length of the byte array is `> EXTENSION_DEGREE_ENCODING_LENGTH`, otherwise return error
- `take(EXTENSION_DEGREE_ENCODING_LENGTH)` and parse it as unsigned integer `extension_degree` for encoding of extension degree for a twist. Only `extension_degree == 2` or `extension_degree == 3` are supported, on other values should return error
- ensure that length of the byte array is `> modulus_length`, otherwise return error
- `take(modulus_length)` and parse it as BigEndian encoding of an unsigned integer that is an `non_residue` - a quadratic or cubic non-residue to make an extension
- ensure `non_residue > 0`, otherwise return error
- check that `non_residue` is non-square or non-cube depending of `extension_degree` to have extension well-formed. If it's not - return error
- ensure that length of the byte array is `> modulus_length*extension_degree`, otherwise return error
- `take(modulus_length*extension_degree)` and parse it as `extension_degree` densely packed BigEndian encodings of unsigned integers that are coefficient of an element in extension field. Coefficients follow from smallest degree: if element is represented as a polynomial `c0 + c1*x + c2*x^2` then coefficients are parsed as `c0`, `c1`, `c2` one after another. That is an `A` coefficient for a curve twist in the Weierstrass form
- ensure that each of `c*` coefficients is `< base_field_modulus`, otherwise return error
- ensure that length of the byte array is `> modulus_length*extension_degree`
- perform similar checks for `B` coefficient
- ensure that `B > 0`, otherwise return error
- ensure that length of the byte array is `> BYTES_FOR_LENGTH_ENCODING`, otherwise return error
- `take(BYTES_FOR_LENGTH_ENCODING)` and parse it as unsigned integer `subgroup_order_length` for encoding of length of the next value
- ensure that `subgroup_order_length > 0`, otherwise return error
- ensure that length of the byte array is `> subgroup_order_length`, otherwise return error
- `take(subgroup_order_length)` and parse it as BigEndian encoding of an unsigned integer `main_subgroup_order`
- ensure that `main_subgroup_order > 0`
- there are two purposes to require `main_subgroup_order` to be included into the ABI:
  - Upfront estimation of the worst-case scenario for a difficulty of the multiplication and multiexponentiation operations
  - One should not enforce point being in a correct subgroup for correctness of operations and if caller expects to process user's inputs he should make such check as a separate call. Nevertheless, this check is MANDATORY for pairing operations
- later should follow operation-specific parameters, depending on curve

Operations in G2 are the same as for G1 with a difference only in encoding of point coordinates now being in the extension of the base field (`modulus_length*extension_degree`).

### ABI and parsing checks for pairing operation

ABI consists of two parts: one defines a base field and a curve, with another that encodes points for a pairing operation.

Signature of public function is just `Operation(byte_array)` (one public function), so just a pointer to the array of bytes is passes to the function.

#### Common ABI and parsing checks for pairing operation

Due to difference in properties and required parameters for different families of the curves first one has to parse curve type:

`CURVE_TYPE_LENGTH = 1`

|Curve | Type |
|------|------|
|BLS12 | 0x01 |
|BN    | 0x02 |
|MNT4  | 0x03 |
|MNT6  | 0x04 |

For generic curves (e.g. generated by Cocks-Pinch method) that support Ate pairing one should use operation codes for MNT4/6 curves depending on `k` parameter of the curve.

- ensure that length of the byte array is `> CURVE_TYPE_LENGTH`, otherwise return error
- `take(CURVE_TYPE_LENGTH)` and parse it as unsigned integer `curve_type` for encoding of the curve type. Check that curve type is known an follow the part that is curve specific

#### ABI for pairing operations on BLS12 curves

Important constants:
|Name | Value |
|------|------|
|TWIST_TYPE_LENGTH| 1 |
|TWIST_TYPE_M| 0x01 |
|TWIST_TYPE_D| 0x02 |
|SIGN_ENCODING_LENGTH| 1 |
|SIGN_PLUS| 0x00 |
|SIGN_MINUS| 0x01 |

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
- `take(x_length)` and parse it as unsigned integer `x`. Empirical testing will put a bound on the sane limits for `x`, that is to be determined. It is unsigned `x` so it should be `> 0`, otherwise return error
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
