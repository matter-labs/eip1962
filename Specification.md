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

### ABI and parsing checks for G1 operations

ABI consists of two parts: one defines a base field and a curve, with another is operation dependent and encodes points or scalars for a corresponding operations. 

Signature of public function is just `Operation(byte_array)` (one public function), so just a pointer to the array of bytes is passes to the function.

#### **Common** ABI and parsing checks for G1 operations

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
- create weierstrass curve instance based on `main_subgroup_order`, `A`, `B` and `base_field_modulus` parameters
- later should follow operation-specific parameters, depending on curve
  
This list is a naive algorithm. Real implementation will merge e.g. reading of `A`, range check and parsing it as a representation of some form of the field element in one operation.

Arithmetic is made using some fixed-length representation of a field element. This implementation follows an approach to represent them as a fixed length array `[u64; NUM_LIMBS`], such that for a modulus `M`: `2^(MODULUS_LIMBS*64) > 2*M` to ensure that one never has to take care about carry flags. In this case a field element with `255` bit modulus would be represented as `[u64; 4]`, but `256` bit modulus will be already represented as `[u64; 5]`

To put some sane limit of the length of the arithmetics modulus `M` must be `< 1024` bits in length, so it's representable as `[u64; 16]` array

#### Notes for G1 ABI:

For purposes of caller's convenience it may be reasonable to have some granularity for lengths of encodings of various element. For example, make encodings multiples of 8 bytes (64 bits) to roughtly correspond to limb bits on x64 machine

### ABI and parsing checks for G2 operations

ABI consists of two parts: one defines a base field and a curve, with another is operation dependent and encodes points or scalars for a corresponding operations. 

Signature of public function is just `Operation(byte_array)` (one public function), so just a pointer to the array of bytes is passes to the function.

#### **Common** ABI and parsing checks for G2 operations

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
- `take(BYTES_FOR_LENGTH_ENCODING)` and parse it as unsigned integer `subgroup_order_length` for encoding of length of the next value
- ensure that `subgroup_order_length > 0`, otherwise return error
- ensure that length of the byte array is `> subgroup_order_length`, otherwise return error
- `take(subgroup_order_length)` and parse it as BigEndian encoding of an unsigned integer `main_subgroup_order`
- ensure that `main_subgroup_order > 0`
- there are two purposes to require `main_subgroup_order` to be included into the ABI:
  - Upfront estimation of the worst-case scenario for a difficulty of the multiplication and multiexponentiation operations
  - One should not enforce point being in a correct subgroup for correctness of operations and if caller expects to process user's inputs he should make such check as a separate call. Nevertheless, this check is MANDATORY for pairing operations
- create weierstrass curve instance based on `main_subgroup_order`, `A`, `B`, `base_field_modulus` and `extension_degree` parameters 
- later should follow operation-specific parameters, depending on curve

Operations in G2 are the same as for G1 with a difference only in encoding of point coordinates now being in the extension of the base field.