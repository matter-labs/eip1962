# Current state of affairs

## Done

- [x] Implementation and toolsets for testing (Rust)
- [ ] Implementation in C++
  - [x] Initial + fuzzy testing on it
  - [ ] Update and check for latest changes from Rust one
  - [ ] Port gas metering
  - [ ] Another run of fuzzy testing + fuzzy test multiexponentiation
- [ ] Gas metering for add/mul/multiexp operations in G1/G2 (WIP)
  - [x] Addition (paitially)
  - [x] Multiplication (partially)
  - [ ] Multiexponentiation
- [x] Gas metering for BLS12 curves pairing
- [x] Gas metering for BN curves pairing
- [ ] Gas metering for MNT4 curves pairing
- [ ] Gas metering for MNT6 curves pairing

## Open questions

- Should we separate cases for `A = 0` and `A != 0` in gas metering for add/mul/multiexp? For now I take the worst case between those two
- Make gas metering multicore. For now it's single core and it takes days to complete
- MNT4/6 pairing operations are diffucult to gas meter due to weak dependencies on some parameters. Now I'm splitting it into three stages and will have to measure them independently:
  - Parse/prepare
  - Miller loop 
  - Final exponentiation
- Gas metering tables show "unintuitive" behavior. I think it's ok and just reflects empirical data. We can extend the data set, but it's already quite redundant and pessimistic: every operation is measured more than once and only the worst result is taken amont them. Then if those results produce overlapping data between each other such data is used for averaging and fitting procedures

