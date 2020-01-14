# Purpose

Rust library for EC arithmetics and pairing calculations over various curves with parameters defined at runtime.

Features (WIP):
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
  - [x] Cocks-Pinch method generated curves in Weierstrass form (Ate pairing)
    - [x] Test over a single k=6 curve from Zexe 

# ABI interface

See [ABI.md](ABI.md).

# Performance testing

- [x] Find a way to save on precomputations
- [x] Benchmark Peppinger

# Contributors

- Kobi Gurkan, [kobigurk@gmail.com](mailto://kobigurk@gmail.com)

# Resources to consult and use 

So I do not forget it

- https://eprint.iacr.org/2012/072.pdf
- https://eprint.iacr.org/2013/722.pdf
- https://eprint.iacr.org/2016/130.pdf