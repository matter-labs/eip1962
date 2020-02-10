# Gas schedule for EIP 1962

Gas schedule is actually the most challenging part of this EIP.

Metering was performed for various operation using either worst case variants (e.g. for point multiplication operation where for a certain number of `uint64` "words" (limbs) in the group order representation all those "words" are filled with maximum values that gives maximum bit length and maximum hamming weight) or Monte-Carlo method (for pairing operations) where some validity checks were disabled (performed, but computation is not stopped if it fails).

Metering was performed on the machine that has `35 MGas/second` performance based on `ecrecover` precompile or `23 MGas/second` based on the current `PAIRING` precompile. Decision was made to use an average of `30 MGas/second` for concrete model estimations. Models for different operations when evaluated give results in `gas` units.

## Calculation of various parameters used in formulas below

- `num_limbs` is a number of limbs required to represent a modulus (affects base arithmetic performance). If `bits` is a bit width of the modulus, than it is calculated as `bits / 64 + 1` (so we do not use top bit and use 5 limbs for 256 bit modulus for ease of implementation)
- `group_order_limbs` is a number of limbs required to represent a group order (affects multiplication and subgroup checks performance). New ABI does not use dense encoding of the group order encoding, so number of words is just `(group_order_length + 7) / 8)`. It in principle allows one to not even parse the encoding as an integer, but to use the length only. ~~If `bits` is a bit width of the group order, than it is calculated as `(bits + 63)/ 64` (different from modulus!)~~
- `num_pairs` parameter represents a either a number of `(point, scalar)` pairs for multiexponentiation or number of `(G1, G2)` pairs for pairing operation
- For loop parameters and parameters used in final exponentiations for pairings one uses bit width and hamming weight of the corresponding parameters. They are covered in the corresponding sections.

## G1 and G2 operations (additions, multiplications, multiexp)

For these types models are simple lookup tables and trivial formulas:

- For `addition` operations in base field, Fp2 or Fp3 there are lookup tables based on `num_limbs` parameter
- For `multiplication` formula is `base + per_limb * group_order_limbs` where for `base` and `per_limb` parameters there are lookup tables based on `num_limbs` parameter
- For `multiexp` formula is `num_pairs * multiplication * discount / discount_multiplier` where `multiplication` is a corresponding price of multiplication, `discount` is gived as a lookup table over `num_pairs` parameter for `num_pairs <= 128` and for larger values there is a `max_discount` cap parameter. `discount_multiplier` is a simple scalar to implement multiplicaiton by parameter `<1` using integer multiplication and integer division

Models are stored in `src/gas_meter/*.json`

## Pairings

Implementation has clear separatation of Miller loop and final exponentiation, so for all the curves final cost is `cost = (one_off + final_exp_cost + num_pairs * miller_loop_cost) / multiplier` where in some cases `one_off` parameter will be equal to zero. `multiplier` is present to have good range of parameters for integer operations instead of floating point calculation.

## Some rationale about the approach to the gas estimation functions

There are various parameters (curve and operation dependent) that apriori define how expensive is an operation. E.g. it's expected that larger `num_pairs` parameter in pairings makes a computation longer and this dependence is *linear* because we just repeat the same computation over larger number of pairs. In a similar manner we can analyze every parameter and end up with a following conclusion that for our range of input parameters for formulas that are NOT the lookup tables (so it's pairing operations over various curves) we can separate the parameters in the following groups

- Affect a computation time in a *linear* manner:
  - `group_order_limbs` 
  - `num_pairs` 
  - `x_bit_length`
  - `x_hamming_weight`
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

Having apriori identified influence of all the parameters on the time of the computation we perform a polynomial fitting using a polynomial of the *total* degree 2 over *set of variables* `v_0, ..., v_n` and chose a set of such variables as listed in each of the corresponding section for the corresponding curve below. Such lists are given in a form `(variable, power)` that means that if `power > 1` we make a new variable `v_k = variable ^ k` for `1 <= k <= power` and use them in the set of variables. Also during polynomial fitting we put an additional restriction that we expect our polynomial to be sparse (for example, there is no computation that's execution time would scale as e.g. `ate_loop_bits * w1_bits`, so we try to enforce this). We also split the full cost of the pairing call to the `base` that performs all the checks and precomputations **before** the actual pairing call. Then we also identify the `miller_loop_cost` cost that we apriori know scales linearly with `num_pairs` (thus we multiply `miller_loop_cost * num_pairs` in the formulas above) and `final_exp_cost` that is also one-shot computation. Such separation allows one to better perform the metering cause e.g. `base` cost depends from the length of the modulus and it's cost can be measured by the greedy algorithm.

## Polynomial model encoding in JSON files

We also need a convention to encode such a polynomial. As was mentioned, it's sparse, and of the total degree 2, so encoding is quite simple:
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

### BSL12

`BLS12` curves are generated from the single parameter `x`. Additional variables for estimation of pairing operation cost are:
- number of bits in `x`, called `x_bit_length` later
- hamming weight of `x`, called `x_hamming_weight` later

There exist models for `miller_loop_cost` and `final_exp_cost` parameters. `one_off` is not used.

- for `miller_loop_cost` input parameters are `(x_bit_length, 1), (x_hamming_weight, 1), (group_order_limbs, 1), (modulus_limbs, 6)`, 
- for `final_exp_cost` input parameters are `(x_bit_length, 1), (x_hamming_weight, 1), (modulus_limbs, 6)`, 

### BN

`BN` curves are generated from the single parameter `x` (or `u` in different literature). Base variables for estimation of pairing operation cost are:
- number of bits in `x`, called `x_bit_length` later
- hamming weight of `x`, called `x_hamming_weight` later
- number of bits in `|6u + 2|`, called `six_u_plus_two_bit_length` later
- hamming weight of `|6u + 2|`, called `six_u_plus_two_hamming` later

There exist models for `miller_loop_cost` and `final_exp_cost` parameters. `one_off` is not used.

- for `miller_loop_cost` input parameters are `(six_u_plus_two_bit_length, 1), (six_u_plus_two_hamming, 1), (group_order_limbs, 1), (modulus_limbs, 6)`, 
- for `final_exp_cost` input parameters are `(x_bit_length, 1), (x_hamming_weight, 1), (modulus_limbs, 6)`, 

### MNT4/MNT6

Those Ate pairings are parametrized by the the miller loop scalar labeled `ate_loop_parameter`, and two contributions into final exponentiation labeled as `w0` and `w1`.

Additional parameters:

- number of bits in `ate_loop_parameter`, called `ate_loop_bits` later
- hamming weight of `ate_loop_parameter`, called `ate_loop_hamming` later
- number of bits in `w0`, called `w0_bits` later
- hamming weight of `w0`, called `w0_hamming` later
- number of bits in `w1`, called `w1_bits` later
- hamming weight of `w1`, called `w1_hamming` later


There exist models for `miller_loop_cost`, `final_exp_cost` and `one_off` parameters

Constant `power` is equal to `4` for MNT4 and to `6` for MNT6

- `one_off` is a simple lookup table based on `modulus_limbs`
- for `miller_loop_cost` input parameters are `(group_order_limbs, 1), (ate_loop_bits, 1), (ate_loop_hamming, 1), (modulus_limbs, power)`,  (NOTE reordering from BN/BLS12 models)
- for `final_exp_cost` input parameters are `(w0_bits, 1), (w0_hamming, 1), (w1_bits, 1), (w1_hamming, 1), (modulus_limbs, power)`, 

## Monte-Carlo simulation rationale

Even some "sane" parameter space is too large to perform full greedy evaluation for a further fitting. For pairing-friendlt curves some parameters were drawn from the space and then deterministically test vectors with `2`, `4` and `6` pairs were generated. Simple linear fit on a final execution time immediately gives final exponentiation and Miller loop (per pair) costs using apriory formula from above.

Even small test set of `4498` data points for model fitting and `500` data point for prediction testing is much larger then the final number of non-sparse coefficients in a fitting polynomial (for `3` input parameters in BLS12 and maximum total(!) degree of `6` in every multivariate term) and, later, sparse coefficients. Normalized `R^2` parameter was `>0.93` for all the fits (training and testing). 

## DDoS prevention

Our polynomial model contains only positive coefficients and returns good statistical fitting results (our `R^2` value is `>0.9` and number of data points is much smaller than the number of coefficients that we need to determine). Simultaneously constant for MGas/second is most likely overestimated compared to the current BN curve precompile that gives us a multiplicative execution time protection. In addition we may also set a bottom limit to call pairing operations to be e.g. `10_000` (arbitrary number) gas no matter what was a set of supplied parameters. Any sane call (using even the most efficient pairing curve with a smallest number of limbs) would never give an execution time below `1ms` that is roughly `30_000` gas`

