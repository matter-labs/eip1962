# Gas schedule for EIP 1962

Gas schedule is actually the most challenging part of this EIP.

All further cost parameters are expressed in microseconds of execution time. It's expected to use a coefficient `15 Mgas/second` for conversion into gas costs unless chosen otherwise.

Also, coefficients of polynomials are floats, but with maximum of `3` digits after decimal point, so all calculations can later be transformed into nanoseconds and integer arithmetic.

## G1 operations

Not published until recalculated with new method

## G2 operations

Not published until recalculated with new method

## Pairings

High-level input parameters:
- number of point pairs for pairing, called `num_pairs`

Implementation has clear separatation of Miller loop and final exponentiation, so for all the curves final cost is `cost = final_exp_cost + num_pairs * miller_loop_cost`. 

### BSL12

`BLS12` curves are generated from the single parameter `x`. Base variables for estimation of pairing operation cost are:
- number of bits in `x`, called `x_bits` later
- hamming weight of `x`, called `x_hamming` later
- number of limbs in representation of the arithmetic, called `limbs` later

Coefficients for formulas below were determined by measuing execution time on Core i9, 2.9 GHz using Monte-Carlo method to determine `x_bits` and `x_hamming`, so performed operation is NOT a correct pairing, but we only care about execution time. Modulus (so, `limbs`) was used from properly generated curves and was selected randomly for each simulation round to cover all the parameter space.

Total time to API call including deserialization was measured. Fitting was performed using polynomial of maximal degree 6 for each of the variables without constant term. Additional L1 regularization was performed to eliminate most of the coefficients and then coefficients smaller then `0.001` were zeroed with fit accuracy recheck.

`final_exp_cost = 0.779 * x_bits^1 * limbs^1 + 2.45 * x_hamming^1 * limbs^1 + 53.569 * limbs^2 + 0.366 * x_bits^1 * limbs^2 + 0.617 * x_hamming^1 * limbs^2 + 4.893 * limbs^3 + 0.009 * x_bits^1 * limbs^3 + 0.042 * limbs^4`

`miller_loop_cost = 0.237 * x_bits^1 * limbs^1 + 0.143 * x_hamming^1 * limbs^1 + 22.11 * limbs^2 + 0.162 * x_bits^1 * limbs^2 + 0.002 * x_hamming^2 * limbs^1 + 0.169 * x_hamming^1 * limbs^2 + 2.074 * limbs^3 + 0.023 * limbs^4`

### BN

Not published until recalculated with new method

### MNT4

Not yet measured

### MNT6

Not yet measured

## Monte-Carlo simulation rationale

Even some "sane" parameter space is too large to perform full greedy evaluation for a further fitting. For pairing-friendlt curves some parameters were drawn from the space and then deterministically test vectors with `2`, `4` and `6` pairs were generated. Simple linear fit on a final execution time immediately gives final exponentiation and Miller loop (per pair) costs using apriory formula from above.

Even small test set of `4498` data points for model fitting and `500` data point for prediction testing is much larger then the final number of non-sparse coefficients in a fitting polynomial (for `3` input parameters in BLS12 and maximum total(!) degree of `6` in every multivariate term) and, later, sparse coefficients. Normalized `R^2` parameter was `>0.93` for all the fits (training and testing). 

