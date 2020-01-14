# Gas schedule for EIP 1962

Gas schedule is actually the most challenging part of this EIP.

Metering was performed for various operation using either worst case variants (e.g. for point multiplication operation where for a certain number of `uint64` "words" in the group order representation all those "words" are filled with maximum values that gives maximum bit length and maximum hamming weight) or Monte-Carlo method (for pairing operations) where some validity checks were disabled (performed, but computation is not stopped if it fails).

Metering was performed on the machine that has `35 MGas/second` performance based on `ecrecover` precompile or `23 MGas/second` based on the current `PAIRING` precompile. Decision was made to use an average of `29 MGas/second` for concrete model estimations. Models for different operations when evaluated give results in `gas` units.

## Calculation of various parameters used in formulas below

- `num_limbs` is a number of limbs required to represent a modulus (affects base arithmetic performance). If `bits` is a bit width of the modulus, than it is calculated as `bits / 64 + 1` (so we do not use top bit and use 5 limbs for 256 bit modulus for ease of implementation)
- `group_order_limbs` is a number of limbs required to represent a group order (affects multiplication and subgroup checks performance). If `bits` is a bit width of the group order, than it is calculated as `(bits + 63)/ 64` (different from modulus!)
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

## Complex model evaluation rules

TODO

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

