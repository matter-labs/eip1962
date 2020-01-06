# Gas schedule for EIP 1962

Gas schedule is actually the most challenging part of this EIP.

All further cost parameters are expressed in microseconds of execution time. It's expected to use a coefficient `15 Mgas/second` for conversion into gas costs unless chosen otherwise.

Also, coefficients of polynomials are floats, but with maximum of `3` digits after decimal point, so all calculations can later be transformed into nanoseconds and integer arithmetic.

## G1 and G2 operations

Published as lookup tables (being updated):
- [Addition](https://github.com/matter-labs/eip1962/G1_and_G2_addition_costs.md)
- [Multiplication](https://github.com/matter-labs/eip1962/G1_and_G2_multiplication_costs.md)
- Multiexp - WIP

## Pairings

High-level input parameters:
- number of point pairs for pairing, called `num_pairs`

Implementation has clear separatation of Miller loop and final exponentiation, so for all the curves final cost is `cost = final_exp_cost + num_pairs * miller_loop_cost`. 

### BSL12

`BLS12` curves are generated from the single parameter `x`. Base variables for estimation of pairing operation cost are:
- number of bits in `x`, called `x_bit_length` later
- hamming weight of `x`, called `x_hamming_weight` later
- number of limbs in representation of the arithmetic, called `modulus_limbs` later

Coefficients for formulas below were determined by measuing execution time on Core i9, 2.9 GHz using Monte-Carlo method to determine `x_bit_length` and `x_hamming_weight`, so performed operation is NOT a correct pairing, but we only care about execution time. Modulus (so, `modulus_limbs`) was used from properly generated curves and was selected randomly for each simulation round to cover all the parameter space.

Total time to API call including deserialization was measured. Fitting was performed using polynomial of maximal degree 6 for each of the variables without constant term. Additional L1 regularization was performed to eliminate most of the coefficients and then coefficients smaller then `0.001` were zeroed with fit accuracy recheck.

`final_exp_cost = 3.942075 * x_bit_length^1 + 1.505175 * x_hamming_weight^1 + 15.113984 * modulus_limbs^1 + 1.865667 * x_bit_length^1 * modulus_limbs^1 + 4.120774 * x_hamming_weight^1 * modulus_limbs^1 + 59.540759 * modulus_limbs^2 + 0.284245 * x_bit_length^1 * modulus_limbs^2 + 0.496913 * x_hamming_weight^1 * modulus_limbs^2 + 5.009175 * modulus_limbs^3 + 0.009092 * x_bit_length^1 * modulus_limbs^3 + 0.016185 * modulus_limbs^4`

`miller_loop_cost = 0.417709 * x_bit_length^1 * modulus_limbs^1 + 0.017594 * x_hamming_weight^1 * modulus_limbs^1 + 18.306288 * modulus_limbs^1 * group_limbs^1 + 0.002166 * x_bit_length^1 * x_hamming_weight^1 * modulus_limbs^1 + 0.144132 * x_bit_length^1 * modulus_limbs^2 + 0.004641 * x_bit_length^1 * modulus_limbs^1 * group_limbs^1 + 0.174109 * x_hamming_weight^1 * modulus_limbs^2 + 4.846864 * modulus_limbs^2 * group_limbs^1 + 0.010811 * modulus_limbs^1 * group_limbs^2 + 0.001364 * modulus_limbs^2 * group_limbs^2`

### BN

`BN` curves are generated from the single parameter `x` (or `u` in different literature). Base variables for estimation of pairing operation cost are:
- number of bits in `x`, called `x_bit_length` later
- hamming weight of `x`, called `x_hamming_weight` later
- number of bits in `|6u + 2|`, called `six_u_plus_two_bit_length` later
- hamming weight of `|6u + 2|`, called `six_u_plus_two_hamming` later
- number of limbs in representation of the arithmetic, called `modulus_limbs` later

Fitting procedure is the same as for `BLS12` curve.

`final_exp_cost =  2.523658 * x_bit_length^1 + 2.862805 * x_hamming_weight^1 + 64.942954 * modulus_limbs^1 + 1.658588 * x_bit_length^1 * modulus_limbs^1 + 1.196372 * x_hamming_weight^1 * modulus_limbs^1 + 58.574598 * modulus_limbs^2 + 0.213338 * x_bit_length^1 * modulus_limbs^2 + 0.328796 * x_hamming_weight^1 * modulus_limbs^2 + 5.034872 * modulus_limbs^3 + 0.002038 * x_bit_length^1 * modulus_limbs^3 + 0.002797 * x_hamming_weight^1 * modulus_limbs^3 + 0.05634 * modulus_limbs^4`

`miller_loop_cost = 0.194238 * six_u_plus_two_bit_length^1 * modulus_limbs^1 + 0.221893 * six_u_plus_two_hamming^1 * modulus_limbs^1 + 18.859488 * modulus_limbs^1 * group_limbs^1 + 0.00177 * six_u_plus_two_bit_length^1 * six_u_plus_two_hamming^1 * modulus_limbs^1 + 0.160538 * six_u_plus_two_bit_length^1 * modulus_limbs^2 + 0.154685 * six_u_plus_two_hamming^1 * modulus_limbs^2 + 4.87253 * modulus_limbs^2 * group_limbs^1`

### MNT4

Not yet measured

### MNT6

Not yet measured

## Monte-Carlo simulation rationale

Even some "sane" parameter space is too large to perform full greedy evaluation for a further fitting. For pairing-friendlt curves some parameters were drawn from the space and then deterministically test vectors with `2`, `4` and `6` pairs were generated. Simple linear fit on a final execution time immediately gives final exponentiation and Miller loop (per pair) costs using apriory formula from above.

Even small test set of `4498` data points for model fitting and `500` data point for prediction testing is much larger then the final number of non-sparse coefficients in a fitting polynomial (for `3` input parameters in BLS12 and maximum total(!) degree of `6` in every multivariate term) and, later, sparse coefficients. Normalized `R^2` parameter was `>0.93` for all the fits (training and testing). 

