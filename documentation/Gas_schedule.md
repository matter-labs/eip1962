# Gas schedule for EIP 1962

Gas schedule is actually the most challenging part of this EIP.

Metering was performed for various operation using either worst case variants (e.g. for point multiplication operation where for a certain number of 64 bit "words" (limbs) in the group order representation all those "words" are filled with maximum values that gives maximum bit length and maximum hamming weight) or Monte-Carlo method (for pairing operations) where some validity checks were disabled (performed, but computation is not stopped if it fails).

Metering was performed on the machine that has `35 MGas/second` performance based on `ecrecover` precompile or `23 MGas/second` based on the current `PAIRING` precompile. Decision was made to use an average of `30 MGas/second` for concrete model estimations. Models for different operations when evaluated give results in `gas` units.

## Note on gas estimation arithmetic

For consistency it's REQUIRED that unsigned 64 bit integer is used for gas cost calculations and that all the arithmetic operations (actually only multiplications and additions) are checked for an overflow. If an overflow is encountered than gas estimator MUST return an error. Also note that gas estimation function MAY return values larger than Ethereum block gas limit. It's an expected behavior and it's not an error, even though such transaction will never get into the block.

## Calculation of various parameters used in formulas below

Gas estimation begins with ABI parsing as described in the corresponding [document](https://github.com/matter-labs/eip1962/blob/master/documentation/ABI.md). During parsing and validation performed by the gas estimator (not all the checks are performed. Checks skipped by the gas estimator are labeled) following common (for many operations) parameters will be parsed:

- `base_field_modulus`
- `group_order_length` (if present)
- `num_pairs` (if present)

Based on these common parameters the following values are calculated to form set of parameters required for metering (not all of them may be present in all the operations, but these are the most common ones):

- `modulus_limbs` is a number of limbs required to represent a modulus (affects base arithmetic performance). If `bits` is a bit width of the `base_field_modulus`, than it is calculated as `bits / 64 + 1` (so we do not use top bit and use 5 limbs for 256 bit modulus for ease of implementation)
- `group_order_limbs` is a number of limbs required to represent a group order (affects multiplication and subgroup checks performance). New ABI does not use dense encoding of the group order encoding, so number of words is just `(group_order_length + 7) / 8)`. It in principle allows one to not even parse the encoding as an integer, but to use the length only. ~~If `bits` is a bit width of the group order, than it is calculated as `(bits + 63)/ 64` (different from modulus!)~~
- `num_pairs` parameter represents a either a number of `(point, scalar)` pairs for multiexponentiation or number of `(G1, G2)` pairs for pairing operation

For pairing operations various loop parameters and parameters used in final exponentiations are parsed from the string. Those are covered in the corresponding sections below.

## G1 and G2 operations (additions, multiplications, multiexp)

For these types models are simple lookup tables and trivial formulas:

- For `addition` operations in base field, Fp2 or Fp3 there are lookup tables based on `modulus_limbs` parameter
- For `multiplication` there are two ways to call the metering:
  - For multiplication and multiexponentiation operation (label it as `multiplication(include_base = true)`) formula is `base + per_limb * group_order_limbs` where for `base` and `per_limb` parameters there are lookup tables based on `modulus_limbs` parameter
  - For subgroup check that is a part of the pairing operations (covered below) one does not need `base` cost (that is one-shot cost not paid during such call) that is labeled as `multiplication(include_base = false)` formula is even simpler: `per_limb * group_order_limbs`
- For `multiexp` formula is `num_pairs * multiplication(include_base = true) * discount / discount_multiplier` where `multiplication(include_base = true)` is a corresponding price of multiplication including(!) one-shot costs, `discount` is given as a lookup table over `num_pairs` parameter for `num_pairs <= 128` and for larger values there is a `max_discount` cap parameter. `discount_multiplier` is a simple scalar to implement multiplicaiton by parameter `<1` using integer multiplication and integer division. Detailed information is given along the model file description at the end of the document

Models are stored in `src/gas_meter/*.json`. Description of the model files is given at the end of this document.

## Pairings

Implementation has clear separatation of Miller loop and final exponentiation, so for all the curves final cost of the pairing operation can be represented as `cost = subgroup_checks + (one_off + final_exp_cost + num_pairs * miller_loop_cost) / multiplier`.
- `subgroup_checks` denoted computations required to perform subgroup checks if requested by caller
- `one_off` parameter may be equal to zero in some cases and denoted costs required to perform all the validations and precomputations (expensive operations like divisions, large exponent powerings, etc). 
- `num_pairs` parameter is factored based on apriori assumptions (more on this below).
- `miller_loop_cost` is a computational cost of running the Miller loop (part of the pairing operation) per single pair of points
- `final_exp_cost` is a computational cost of the final exponentiation. This is one-off operation after the Miller loop
- `multiplier` is present to have good range of parameters for integer operations instead of floating point calculation. 

For each of the curve families there will be given explicit sets of parameters used to calculate corresponding contributions `one_off`, `miller_loop_cost` and `final_exp_cost`.

## Some rationale about the approach to the gas estimation functions

There are various parameters (curve and operation dependent) that apriori define how expensive is an operation. E.g. it's expected that larger `num_pairs` parameter in pairings makes a computation longer and this dependence is *linear* because we just repeat the same computation over larger number of pairs. In a similar manner we can analyze every parameter and end up with a following conclusion that for our range of input parameters for formulas that are NOT the lookup tables (so it's pairing operations over various curves) we can separate the parameters in the following groups

- Affect a computation time in a *linear* manner (but we need to determine a prefactor anyway!):
  - `group_order_limbs` 
  - `num_pairs` 
  - `x_bit_length`
  - `x_hamming_weight`
  - `u_bit_length`
  - `u_hamming_weight`
  - `six_u_plus_two_bit_length`
  - `six_u_plus_two_hamming`
  - `ate_loop_bits`
  - `ate_loop_hamming`
  - `w0_bits`
  - `w0_hamming`
  - `w1_bits`
  - `w1_hamming`
- Affect a computation in non-linear manner:
  - `modulus_limbs`

Thus we choose an approach to use a polynomial fitting to calculate the pairing call gas costs. 

### Polynomial evaluation

Let's use a simple example of the polynomial of the *total* degree 2 over *two* variables `x` and `y`. It has the following form: `p(x, y) = c0 + c01 * x + c10 * y + c11 * x * y + c02 * x^2 + c20 * y^2`. Such polynomial can be sparse if a lot of coefficients `c_` are equal to zero. To evaluate this polynomial we can precompute power series of `x` and `y` up to the 2nd power inclusive and perform a substitution.

## Gas cost analysis

Having apriori identified influence of all the parameters on the time of the computation we perform a polynomial fitting using a polynomial of the *total* degree 2 over *set of variables* `v_0, ..., v_n` and chose a set of such variables as listed in each of the corresponding section for the corresponding curve below. Such lists are given in a form `(variable, power)` that means that if `power > 1` we make a new variable `v_k = variable ^ k` for `1 <= k <= power` and use them in the set of variables. Also during polynomial fitting we put an additional restriction that we expect our polynomial to be sparse (for example, there is no computation that's execution time would scale as e.g. `ate_loop_bits * w1_bits`, so we try to enforce this). We also split the full cost of the pairing call to the `base` that performs all the checks and precomputations **before** the actual pairing call. Then we also identify the `miller_loop_cost` cost that we apriori know scales linearly with `num_pairs` (thus we multiply `miller_loop_cost * num_pairs` in the formulas above) and `final_exp_cost` that is also one-shot computation. Such separation allows one to better perform the metering cause e.g. `base` cost depends from the length of the modulus (`modulus_limbs` parameter) only and it's cost can be measured by the greedy algorithm.

## Lookup tables encoding in JSON files

As many parameters in gas estimation formulas are taken from the lookup tables, here we give a description of the encoding.

JSON field `json_field` contatining a lookup table for a variable `var` over parameter `param` is encoded as a vector (`[...]`) of variables of `data_point` type. `data_point` type is a vector (`[...]`) of two integer(!) elements: [`param_value`, `var_value`] where `param_value` indicates for what value of `param` this `var_value` is returned for a variable `var`.

## Polynomial model encoding in JSON files

Full description of the model files for each of the operations is given at the end of the document. Here we only cover in details how encoding of the sparse polynomial is performed cause it's the only non-trivial part.

We need a convention to encode such a polynomial. As was mentioned, it's sparse, and of the total degree 2, so encoding is quite simple:
- we encode a vector (`[...]` in JSON) of coefficients
- each coefficient is a vector (`[...]` again) of two elements: `coefficient_value` and `variables`
- `coefficient_value` is an (unsigned) integer
- `variables` is a vector (`[...]` again) of pairs `(var_num, power)`. Pairs are encoded as vectors of length 2.
- `var_num` is zero-enumerated index of the variable in the corresponding list for each curve. E.g. `var_num = 1` for BLS12 curve for `miller_loop_cost` means `group_order_limbs`
- `power` is a power of the value of the corresponding variable

These information is enough to evaluate every subterm of the polynomial.

Example: for BLS12 curve we have a term `[6309, [[1, 1], [2, 2]]]` in the `final_exp_cost`
- First we parse it at `(coefficient_value, variables)`: 
  - `coefficient_value = 6309`
  - `variables = [[1, 1], [2, 2]]`
- Now for a list of `variables`
  - First `(var_num, power) = [1, 1]` means `x_hamming_weight` and power of 1
  - Second `(var_num, power) = [2, 2]` means `modulus_limbs` and power of 2
- Subterm will look as `6309 * x_hamming_weight^1 * modulus_limbs^2`

NOTE: don't forget to divide by `multiplier` in the full evaluation!

### Common parameters for all the pairing operations

As a result of ABI parsing one gets values of two variables:
- `num_g1_checks`
- `num_g2_checks`

that indicated for how many points in G1/G2 we need to perform the subgroup check. Subgroup check is a trivial multiplication by the group order, so for all the pairing operations `subgroup_checks = multiplication_in_g1(include_base = false) * num_g1_checks + multiplication_in_g2(include_base = false) * num_g2_checks` where `multiplication_in_g1` and `multiplication_in_g2` are `multiplication` costs based on the corresponding models (namely multiplication in Fp, Fp2 or Fp3). In both cases `base` cost of the multiplication is NOT included

### BSL12

`BLS12` curves are generated from the single parameter `x`. 

Additional variables for estimation of pairing operation cost are:
- number of bits in the `x` parameter (from ABI decoding), called `x_bit_length` later
- hamming weight of the `x` parameter (from ABI decoding), called `x_hamming_weight` later

Model for BLS12 curve pairings is located in the JSON file named `bls12_model.json`.

- `one_off` cost is not used and not present in the model
- for `miller_loop_cost` input parameters are `(x_bit_length, 1), (x_hamming_weight, 1), (group_order_limbs, 1), (modulus_limbs, 6)`, 
- for `final_exp_cost` input parameters are `(x_bit_length, 1), (x_hamming_weight, 1), (modulus_limbs, 6)`
- `multiplication_in_g1` is based on the model file `g1_multiplication.json`
- `multiplication_in_g2` is based on the model file `g2_multiplication_ext2.json`

Model file (JSON) contain the following fields:
- `multiplier` - single integer encoding `multiplier`
- `miller` - encoding of the `miller_loop_cost` polynomial model
- `final_exp` - encoding of the `final_ext_cost` polynomial model

### BN

`BN` curves are generated from the single parameter `x` (or `u` in different literature). 

Base variables for estimation of pairing operation cost are:
- number of bits in `u` (from ABI decoding), called `u_bit_length` later
- hamming weight of `u` (from ABI decoding), called `u_hamming_weight` later
- number of bits in `|6u + 2|` (calculated during ABI decoding), called `six_u_plus_two_bit_length` later
- hamming weight of `|6u + 2|` (calcualted during ABI decoding), called `six_u_plus_two_hamming` later

Model for BLS12 curve pairings is located in the JSON file named `bn_model.json`.

- `one_off` cost is not used and not present in the model
- for `miller_loop_cost` input parameters are `(six_u_plus_two_bit_length, 1), (six_u_plus_two_hamming, 1), (group_order_limbs, 1), (modulus_limbs, 6)`, 
- for `final_exp_cost` input parameters are `(u_bit_length, 1), (u_hamming_weight, 1), (modulus_limbs, 6)`, 
- `multiplication_in_g1` is based on the model file `g1_multiplication.json`
- `multiplication_in_g2` is based on the model file `g2_multiplication_ext2.json`

Model file (JSON) contain the following fields:
- `multiplier` - single integer encoding `multiplier`
- `miller` - encoding of the `miller_loop_cost` polynomial model
- `final_exp` - encoding of the `final_ext_cost` polynomial model

### MNT4/MNT6

Those Ate pairings are parametrized by the the miller loop scalar labeled `ate_loop_parameter`, and two contributions into final exponentiation labeled as `w0` and `w1`.

Additional parameters:

- number of bits in `ate_loop_parameter` (from ABI decoding), called `ate_loop_bits` later
- hamming weight of `ate_loop_parameter` (from ABI decoding), called `ate_loop_hamming` later
- number of bits in `w0` (from ABI decoding), called `w0_bits` later
- hamming weight of `w0` (from ABI decoding), called `w0_hamming` later
- number of bits in `w1` (from ABI decoding), called `w1_bits` later
- hamming weight of `w1` (from ABI decoding), called `w1_hamming` later


Models for MNT4/MNT6 curve pairings is located in the JSON files named `mnt4_model.json` and `mnt6_model.json` respectively.

Constant `power` is equal to `4` for MNT4 and to `6` for MNT6

- `one_off` is a simple lookup table based on `modulus_limbs`
- for `miller_loop_cost` input parameters are `(group_order_limbs, 1), (ate_loop_bits, 1), (ate_loop_hamming, 1), (modulus_limbs, power)`,  (NOTE reordering from BN/BLS12 models)
- for `final_exp_cost` input parameters are `(w0_bits, 1), (w0_hamming, 1), (w1_bits, 1), (w1_hamming, 1), (modulus_limbs, power)`, 
- `multiplication_in_g1` is based on the model file `g1_multiplication.json`
- `multiplication_in_g2` is based on the model file `g2_multiplication_ext2.json` for MNT4 and on the model file `g2_multiplication_ext3.json` for MNT6.

Model files themselfves (JSONs) contain the following fields:
- `one_off` - lookup table for `one_off` based on `modulus_limbs`
- `multiplier` - single integer encoding `multiplier`
- `miller` - encoding of the `miller_loop_cost` polynomial model
- `final_exp` - encoding of the `final_ext_cost` polynomial model

## Monte-Carlo simulation rationale

Even some "sane" parameter space is too large to perform full greedy evaluation for a further fitting. For pairing-friendlt curves some parameters were drawn from the space and then deterministically test vectors with `2`, `4` and `6` pairs were generated. Simple linear fit on a final execution time immediately gives final exponentiation and Miller loop (per pair) costs using apriory formula from above.

Even small test set of `4498` data points for model fitting and `500` data point for prediction testing is much larger then the final number of non-sparse coefficients in a fitting polynomial (for `3` input parameters in BLS12 and maximum total(!) degree of `6` in every multivariate term) and, later, sparse coefficients. Normalized `R^2` parameter was `>0.93` for all the fits (training and testing). 

## DDoS prevention

Our polynomial model contains only positive coefficients and returns good statistical fitting results (our `R^2` value is `>0.9` and number of data points is much smaller than the number of coefficients that we need to determine). Simultaneously constant for MGas/second is most likely overestimated compared to the current BN curve precompile that gives us a multiplicative execution time protection. In addition we may also set a bottom limit to call pairing operations to be e.g. `10_000` (arbitrary number) gas no matter what was a set of supplied parameters. Any sane call (using even the most efficient pairing curve with a smallest number of limbs) would never give an execution time below `1ms` that is roughly `30_000` gas`


## Model files description

### Addition model files

These model files are
- `g1_addition.json` (for curve defined over the base field Fp)
- `g2_addition_ext2.json` (for curve defined over the quadratic extension field Fp2)
- `g2_addition_ext3.json` (for curve defined over the cubic extension field Fp3)

and encode simple lookup tables based on the `modulus_limbs` parameter.

JSON contains a single field `price` that a lookup table for variable `price` over the parameter `modulus_limbs`.

Addition cost is simply `cost = price` from the corresponding lookup table.

Example: for a BN254 curve and operations in G1 (that has a modulus 254 bits long and thus uses `modulus_limbs = 4`) corresponding `price = 195`

### Multiplication model files

These model files are
- `g1_multiplication.json` (for curve defined over the base field Fp)
- `g2_multiplication_ext2.json` (for curve defined over the quadratic extension field Fp2)
- `g2_multiplication_ext3.json` (for curve defined over the cubic extension field Fp3)

and encode simple lookup tables based on the `modulus_limbs` parameter.

JSON contains a two fields: `base` that a lookup table for variable `base` over the parameter `modulus_limbs` and `per_limb` that a lookup table for variable `per_limb` over the parameter `modulus_limbs`.

Full cost is calcualted as `cost = base + per_limb * group_order_limbs` where `base` and `per_limb` values are taken from the corresponding lookup tables described above.

NOTE: please refer to the beginning of the document for multiplication cost routine description with and without base value)

Example: for a BN254 curve and multiplication in G1 (that has a modulus 254 bits long and thus uses `modulus_limbs = 4`) corresponding `base = 180` and `per_limb = 870`. Caller has encoded scalar using `32` bytes that results in `group_order_limbs = 4`. Thus multiplication cost is `180 + 4*870 = 3660`

### Multiexponentiations discounts file

There is a single file describing discount when performing multiexponentiation compared to individual multiplication.
As described above first gas estimator calculated the corresponding single multiplication price `multiplication(include_base = true)` based on whether multiplication is over Fp/Fp2/Fp3. Then based on the single model file `multiexp_discounts.json` one lookups a discount based on the `num_pairs` parameter from the ABI.

There is only a single model file for all the Fp/Fp2/Fp3 for simplicity. Discounts include worst case discounts between Fp, Fp2 and Fp3 operations (that nevertheless had deviations below `10%`).

JSON files contatains four field
- `discount_multiplier` - integer that encodes a `multiplier` parameter in the multiexponentiation formula
- `max_pairs` - integer that describes what is a largest value of `num_pairs` for which the discount is given as the lookup instead of the fixed cap `max_discount`
- `max_discount` - integer that encodes the discount cap (in case of `num_pairs > max_pairs`)
- `discounts` that a lookup table for variable `discount` over the parameter `num_pairs` when `num_pairs <= max_pairs`

After ABI parsing and obtaining `num_pairs` gas estimator first calculated a `discount` by either lookup table if `num_pairs <= max_pairs` or using `max_discount` if `num_pairs > max_pairs`. Please remember that `num_pairs = 0` is forbidden by the ABI and lookup table gives a discount starting from `num_pairs = 1`, so this function is well-defined.

With `multiplication(include_base = true)` calculated as described at the beginning of this section one calculated `multiexp_cost = num_pairs * multiplication(include_base = true) * discount / discount_multiplier`.

Example: for a BN254 curve and multiexponentiation of `10` scalar-point pairs in G1 corresponding `multiplication(include_base = true)` cost is equal to `3660` (in a section above). `discount_multiplier = 1000` from the model file and lookup over `num_pairs = 10` provides `discount = 380`. Thus final price is `10 * 3660 * 380 / 1000 = 13908` gas (please note the floor integer division).