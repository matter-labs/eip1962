# Purpose

Rust library for EC arithmetics and pairing calculations over various curves with parameters defined at runtime.

Features (WIP):
- [x] Fields implementation
- [ ] Weierstrass curves implementation
  - [x] a = 0
  - [x] generic case (a != 0, b != 0)
  - [ ] b = 0 
  - [ ] a = -3
- [ ] Extension towers
  - [x] Fp2
  - [x] Fp3
  - [ ] Fp4 as 2 over 2
  - [x] Fp6 as 2 over 3
  - [x] Fp6 as 3 over 2
  - [x] Fp12 as 2 over 3 over 2
- [ ] Pairings
  - [x] BLS12 curves family
  - [ ] BN family
  - [ ] MNT6 family
  - [ ] MNT4 family
  - [ ] Cocks-Pinch method generated curves in Weierstrass form
 
