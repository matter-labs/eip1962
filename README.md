# Purpose

Rust library for EC arithmetics and pairing calculations over various curves with parameters defined at runtime.

Features (WIP):
- [x] Fields implementation
- [x] Weierstrass curves implementation
  - [x] a = 0
  - [x] generic case (a != 0, b != 0)
  - [ ] b = 0 (most likely will not be implemented to avoid point (0,0) being on curve)
  - [ ] a = -3 (not a priority, can be covered by generic case w/o much performance hit and with simpler gas cost schedule)
- [x] Extension towers
  - [x] Fp2
  - [x] Fp3
  - [x] Fp4 as 2 over 2
  - [x] Fp6 as 2 over 3
  - [x] Fp6 as 3 over 2
  - [x] Fp12 as 2 over 3 over 2
- [ ] Pairings
  - [x] BLS12 curves family
  - [x] BN family
  - [x] MNT6 family
  - [x] MNT4 family
  - [ ] Cocks-Pinch method generated curves in Weierstrass form (Ate pairing)
    - [x] Test over a single k=6 curve from Zexe 
 
# Performance testing

- [ ] Find a way to save on precomputations
- [ ] Find more test vectors to fit quadratic gas schedules
- [ ] Benchmark Peppinger


# Resources to consult and use 

So I do not forget it

- https://eprint.iacr.org/2012/072.pdf
- https://eprint.iacr.org/2013/722.pdf
- https://eprint.iacr.org/2016/130.pdf