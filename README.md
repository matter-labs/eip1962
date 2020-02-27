# Status

This Rust implementation of EIP1962 is complete to the large extend. If course it's possible to polish further (e.g. make it `no_std` compatible), but largest part is done:

Features:
- [x] Fields implementation
- [x] Weierstrass curves implementation
  - [x] a = 0
  - [x] generic case (a != 0, b != 0)
- [x] Extension towers
  - [x] Fp2
  - [x] Fp3
  - [x] Fp4 as 2 over 2
  - [x] Fp6 as 2 over 3
  - [x] Fp6 as 3 over 2
  - [x] Fp12 as 2 over 3 over 2
- [x] Pairings
  - [x] BLS12 curves family
  - [x] BN family
  - [x] MNT6 family
  - [x] MNT4 family
  - [x] Cocks-Pinch method generated curves in Weierstrass form (Ate pairing) with k=6

Testing:

- Basic properties are tested during development (whitebox testing) in a form of e.g. bilinearity checks for pairings
- Fuzzy testing in cross-checks mode with C++ and Go implementations that catches both crashes in any of the libraries and tests for a consistent output (for consensus purposes) 
  - During such testing most of the checks are disabled, e.g. points are allowed to be not on the curve cause it would be difficult for a fuzzer to find a proper test vector. So such testing covers more edge cases then would be possible in production

# Documentation about EIP1962

See [documentation](https://github.com/matter-labs/eip1962/tree/master/documentation) folder for a complete description and the single source of truth about EIP.

## Original proposal

Original EIP is [here](https://eips.ethereum.org/EIPS/eip-1962)

# Contributors

- Kobi Gurkan, [kobigurk@gmail.com](mailto://kobigurk@gmail.com)

# Resources to consult and use 

- https://eprint.iacr.org/2012/072.pdf
- https://eprint.iacr.org/2013/722.pdf
- https://eprint.iacr.org/2016/130.pdf