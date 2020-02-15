# Document overview

This document gives a high-level overview what is a functionality of EIP 1962 precompile and gives brief examples of it's practical applications, as well as explanation of the main challenges to design and implement this precompile. 

IMPORTANT NOTE: We try to provide this set of documents in a form of both rationale (why and how we made certain design decision) and a spec (how to execute this decision in code).

### External source/tutorial

There are terms that are specific to cryptography and we invite the reader to refer to the [Handbook of Applied Cryptography](http://cacr.uwaterloo.ca/hac/) to find expanations of the terms the he/she does not know or understand well. This book is a huge amount of work and we don't think it's a good idea to copy/paste parts of it to form a short glossary for this spec.

## Benefits of this precompile and what is included

- Support of basic arithmetic (additions and multiplications) for elliptic curves over prime field and it's quadratic and cubic extensions.
- Multiexponentiation operation that allows one to save execution time (so gas) by both the algorithm used and (usually forgotten) by the fact that `CALL` operation in Ethereum is expensive (at the time of writing), so one would have to pay non-negigible overhead if e.g. for multiexponentiation of `100` points would have to call the multipication precompile `100` times and addition for `99` times (roughly `138600` would be saved).
- Pairing operations over wide set of curves
  - BLS12 family. For now it's de-facto standard for all cryptographic implementations cause it provides a required security levels with high efficiency
  - BN family. It is provided both for a legacy reasons (BN254 curve precompile exists in Ethereum) and that such curves can be used to make half-pairing-friendly curve cycles
  - Weierstrass curves with Ate pairing and embedding degree k=4 or k=6. Later in the precompile it's called MNT4/6 curves. Such curves are used as wrapper curves for one-step-recursive constructions (e.g. Zexe) and for full recursive cycles (recursive SNARKs)
- Three different implementations
  - Rust (this repo) that is written from scratch including all the arithmetic.
  - C++ ([repo](https://github.com/matter-labs/eip1962_cpp)) that uses external implementation of field operations from [ctbignum](https://github.com/niekbouman/ctbignum) modern arithmetic library
  - Go ([repo](https://github.com/saitima/eip1962)) that is a completely separate and independent work
- Implementations were fuzzy-tested for months on 32 core server both individually for crashes and pairwise for consistency
- Gas metering scripts - in a contrast with previous EIPs where for gas metering one would have to write his own scripts or benchmarks in this repo there is a full set of benchmarking tests and data processing scripts in Python that perform an analysis and output complete model files that can be used by gas estimator

## Examples of what can be done using Ethereum and this precompile

- BLS signatures can not be used with aggregation of public keys in both G1 or G2 (don't be confused with BLS12 curve!) and can be made with 128 bits of security by utilization of new elliptic curves (existing BN254 precompile provides 80 bits by latest estimates). Applications of this can be e.g. mutlisignatures, sidechains or DAOs
- Feature reach privacy solutions as described in Zexe and was demonstrated by recent work of EY can be implemented using BLS12-377 curve and it's embeddings that use Ate pairing with k=6 embedding degree
- Cheaper verification of Bulletproofs due to multiexponentiation operation both for range proofs and privacy, as well as for arithmetic circuits
- Now it's possible to implement Schnorr signatures verification and various modifications including off-chain aggregations/multisignatures like MuSig and numerous variants
- Arithmetic over Secp256k1 can now be performed in full without tricks with `ecrecover` precompile
- Pedersen arithmetic hashing (usually used inside of the SNARKs, also refered as using jubjub/babyjubjub curves) can be very efficiently implemented using multiexponentiation routine and curve equivalence between Twisted Edwards form and Weierstrass form used by this precompile
- Recursive SNARKs verification using both half-pairing-friendly cycles and full cycles

## What is different from precompile that focuses on a single curve

To have more capabilities we also have to pay a price. For this precompile this price is actually very low and mainly results in longer execution time cause some values can not be precomputed

- During construction of the finite field from the modulus one has to perform few long divisions that results in tens of microseconds of time to construct a field. But this is done only once and results in a fixed gas overhead
- For pairing operations one also has to precompile Frobenius endomorphism coefficients. This overhead applies only to the pairing operation and involves few quite long exponentiations. But again, this is also one-off cost
- Potential lack of multiplicative inverse that is required for some operations. In the section below one can find a document with all the explicit arithmetic formulas and conventions how this case is covered in two places where it can happen: during simple arithmetic (addition/multiplication) and during pairing operations.

## Consistency and validity

For Ethereum precompile it's important (especially if there is more than one implementation) to be crash-resistant and input-output consistent to preserve the network consensus. Separate important aspect is a gas schedule that should not allow DDoS attacks and it's covered in the separate gas schedule document (all the links are at the end of this document).

Here we focus on input-output consistency and cover main key points that explain why **if** ABI parsing was performed consistently (so values are passed to the arithmetic routines are the same in any implementation of such routines) **then** outcome of those routines will always be consistent if we follow very narrow set of conventions (like "double-and-add" multiplications for elliptic curves instead of using NAF/wNAF forms).

- As covered in a document with explicit implementation details and formulas all the arithmetic in the precompile *must* be performed in Montgomery form as it's the only efficient way of implementation. Montgomery form requires to use a constant usually labeled as `R = 2^k` where `k` is usually alligned to be a multiple of machine word (64 bits in out case) and `k >= ceil(log2(modulus))`. It's also required that `gcd(R, modulus) = 1`. This requirement is covered due to ABI check that modulus is odd, while `R` obviously has only `1` and `2` as divisors.
- After field is constructed and Montgomery form coefficients are precomputed all arithmetic operations are **always** consistent up to existence of the inverse elements that is covered separately (if implemented correctly of course. This is well demonstrated by the Rust and C++ implementations that use two completely different approaches to implementation of the field operaitons).
- Extension fields are constructed in some precompile calls. Even while the same document with expicit formulas specifies how one *can* implement arithmetic in such extensions, one does not have to follow it because every *valid* implementation will give consistent results: extension field arithmetic is in it's nature largely a polynomial multiplication where coefficients of the polynomials are field elements and field arithmetic *is* consistent from the previous point.
- For elliptic curve operations we present explicit formulas and conventions how one *has to* implement them for ease of consistency. Even while every set of user input parameters defines *some* elliptic curve and supplied points are always checked to be on the curve (so addition and multiplication operations *are* well defined), we've decided to provide explicit formulas due to large number of them in existence. We also agree on convention that multiplication is performed as "double-and-add".
- Same argument holds for pairing operations where we do provide explicit formulas and conventions.

As a result **if** all the input values are the same and **if** implemented field arithmetic algorithms are valid **then** output will always be consistent between different implementations.

## Documents structure

There are three main documents that describe this EIP and allows to one to implement it from scratch:
- [ABI](https://github.com/matter-labs/eip1962/blob/master/documentation/ABI.md) that describes parsing and validation specification. Even though such validation performs some *field* arithmetic operations, consistency with this spec is a largest part of work.
- [Gas schedule](https://github.com/matter-labs/eip1962/blob/master/documentation/Gas_schedule.md) that describes gas schedule approach, specification with explicit formulas, description of the gas model files and examples.
- [Explicit arithmetic formulas](https://github.com/matter-labs/eip1962/blob/master/documentation/Algorithms_for_EIP1962.pdf) that contains in a single place formulas from various sources that were used to implement this EIP. It also contains remarks about conventions chosen (e.g. we use "double-and-add" multiplication and "square-and-multiply" powering). Please note that at the time of writing this document is updated constantly

## Downgrading

It's possible to downgrade this EIP for a gradual rollout in two stages:
- First use a whitelist of the most necessary curves (such whitelisting is prefix matching in our ABI)
- Then lift the restriction

Whitelisting at the first stage also would allow one to provide specialized implementations for each curve for simplicity, although the list of curves is large.

NOTE: it's also possible to even assign a separate precompile address to each specific whitelisted curve, but it's a large set of addresses! 6 addresses per pairing-friendly curve and 3 for normal ones.

List of the curves most wanted at the current state of affairs:
- BN254 (same as the current Ethereum precompile): get multiexponentiation, G2 operations for BLS signature aggregation, separate subgroup check from pairing cause it's not required in most of the cases
- BLS12-381 - de-factor standard for 120+ bits of security. Used by zCash, Eth2.0, other
- BLS12-377 - unique curve from ZEXE paper that has large number of roots of unity for both group order and base field (allows efficient wrapping)
- Some Weierstrass curve with k=6 for wrapping of BLS12-377: either an original or updated from Zexe, or recent works by EY (more efficient)
- MNT4/6 cycle with 750 bit fields: questionable. This is an old curve used by Coda, most likely is too heavy (used for resursive SNARKs)
- BN/Weierstrass half-pairing-friendly curve cycle: for Halo applications
- Secp256k1: for multisignatures, Schnorr signatures, Bulletproofs (no pairings for this one)
- Rational equivalent to Twisted Edwards curve used for Pedersen hashing by BN254 based SNARKS (no pairings for this one)
- Rational equivalent to Twisted Edwards curve used for Pedersen hashing by BLS12-381 based SNARKS (no pairings for this one)
- Rational equivalent to Twisted Edwards curve used for Pedersen hashing by BLS12-377 based SNARKS (no pairings for this one)

There are other flavors of curves that readers can check out and may be suggest to include (using e.g. a PR):
- Twisted Edwards curve over `2^255 - 19` field (`Ed25519`, used by `dalek-crypto`) - also de-factor standard
- Curves used by goverment service (e.g. Estonia's e-resident card) - to verify authority-issued signatures
