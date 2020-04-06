pub(crate) mod pairings;
pub(crate) mod g2_ops;
pub(crate) mod g1_ops;
pub(crate) mod parsers;
pub(crate) mod public_api;
pub(crate) mod spec_generator;
pub(crate) mod arithmetic_tests;

mod fields;
// mod fuzzing;
mod gas_meter;

use num_bigint::BigUint;
use num_traits::Zero;
use num_traits::cast::ToPrimitive;

use crate::errors::ApiError;

pub(crate) fn num_limbs_for_modulus(modulus: &BigUint) -> Result<usize, ApiError> {
    use crate::field::calculate_num_limbs;

    let modulus_limbs = calculate_num_limbs(modulus.bits()).map_err(|_| ApiError::InputError(format!("Modulus is too large, file {}, line {}", file!(), line!())) )?;

    Ok(modulus_limbs)
}

pub(crate) fn num_units_for_group_order(order: &BigUint) -> Result<usize, ApiError> {
    let limbs = (order.bits() + 63) / 64;
    if limbs > 16 {
        return Err(ApiError::InputError(format!("Group order is too large, file {}, line {}", file!(), line!())));
    }

    Ok(limbs)
}

pub(crate) fn calculate_num_limbs(modulus: &BigUint) -> Result<usize, ()> {
    let bitlength = modulus.bits();

    let mut num_limbs = (bitlength / 64) + 1;
    if num_limbs < 4 {
        num_limbs = 4;
    }

    if num_limbs > 16 {
        return Err(());
    }

    Ok(num_limbs)
}
    
pub(crate) fn biguint_to_u64_vec(mut v: BigUint) -> Vec<u64> {

    let m = BigUint::from(1u64) << 64;
    let mut ret = Vec::with_capacity((v.bits() / 64) + 1);

    while v > BigUint::zero() {
        ret.push((&v % &m).to_u64().expect("is guaranteed to fit"));
        v >>= 64;
    }

    ret
}

#[cfg(test)]
mod test {
    #[test]
    #[ignore]
    fn benchmark_ecrecover() {
        use std::thread;

        use std::sync::mpsc::{channel, TryRecvError};

        use rayon::prelude::*;

        use rand::{Rng, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        const RUNS_PER_WORK_UNIT: usize = 100000;
        const PARALLEL_WORKS: usize = 50;

        use parity_crypto::publickey::{Signature, recover as ec_recover};
        use ethereum_types::{H256};
        use keccak_hash::keccak;


        // let mut multiexp_len = vec![0, 2, 4, 8, 16, 32, 64, 128];
        let mut multiexp_len = vec![2, 4, 8, 16, 32, 64, 128];
        multiexp_len.reverse();

        use indicatif::{ProgressBar, ProgressStyle};

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        let mut parameters_space = vec![];
        let (tx, rx) = channel();

        for i in 0..PARALLEL_WORKS {
            parameters_space.push((i, rng.clone(), pb.clone(), tx.clone()));
        }

        drop(tx);

        pb.set_length((parameters_space.len() * RUNS_PER_WORK_UNIT) as u64);

        let handler = thread::spawn(move || {
            parameters_space.into_par_iter().for_each(|(i, mut rng, pb, tx)| {
                for _ in 0..RUNS_PER_WORK_UNIT {
                    let mut r: u64 = rng.gen();
                    r = r.wrapping_add(i as u64);
                    let i = if r & 1 == 1 {
                        hex::decode("38d18acb67d25c8bb9942764b62f18e17054f66a817bd4295423adf9ed98873e000000000000000000000000000000000000000000000000000000000000001b38d18acb67d25c8bb9942764b62f18e17054f66a817bd4295423adf9ed98873e789d1dd423d25f0772d2748d60f7e4b81bb14d086eba8e8e8efb6dcff8a4ae02").unwrap()
                    } else {
                        hex::decode("47173285a8d7341e5e972fc677286384f802f8ef42a5ec5f03bbfa254cb01fad000000000000000000000000000000000000000000000000000000000000001b650acf9d3f5f0a2c799776a1254355d5f4061762a237396a99a0e0e3fc2bcd6729514a0dacb2e623ac4abd157cb18163ff942280db4d5caad66ddf941ba12e03").unwrap()
                    };
                    // this is a ECRECOVER computation without writing a result
                    let len = std::cmp::min(i.len(), 128);

                    let mut input = [0; 128];
                    input[..len].copy_from_slice(&i[..len]);

                    let hash = H256::from_slice(&input[0..32]);
                    let v = H256::from_slice(&input[32..64]);
                    let r = H256::from_slice(&input[64..96]);
                    let s = H256::from_slice(&input[96..128]);

                    let bit = match v[31] {
                        27 | 28 if v.0[..31] == [0; 31] => v[31] - 27,
                        _ => { panic!("invalid input"); },
                    };

                    let start = std::time::Instant::now();
                    let s = Signature::from_rsv(&r, &s, bit);
                    let mut o = [0u8; 32];

                    use std::io::Write;

                    assert!(s.is_valid());
                    if let Ok(p) = ec_recover(&s, &hash) {
                        let r = keccak(p);
                        (&mut o[12..]).write(&r.as_bytes()[12..]).unwrap();
                    }

                    let elapsed_nanos = start.elapsed().as_nanos();
                    tx.send((elapsed_nanos, o)).unwrap();

                    pb.inc(1);
                }
            });
        });

        let mut sum = 0u128;

        loop {
            let subres = rx.try_recv();
            match subres {
                Ok(r) => {
                    sum += r.0;
                },
                Err(TryRecvError::Empty) => {
                    std::thread::sleep(std::time::Duration::from_millis(1000u64));
                },
                Err(TryRecvError::Disconnected) => {
                    handler.join().unwrap();
                    break;
                }
            }
        }

        pb.finish();

        let ns_per_run = sum / ((PARALLEL_WORKS * RUNS_PER_WORK_UNIT) as u128);

        const ECRECOVER_GAS: u128 = 3000;

        let gas_per_second = 1_000_000 * 1_000 * ECRECOVER_GAS / ns_per_run;

        println!("MGAS per second on the current PC = {}", (gas_per_second as f64) / 1_000_000f64);
    }

    #[test]
    #[ignore]
    fn benchmark_existing_pairing_precompile() {
        use std::thread;

        use std::sync::mpsc::{channel, TryRecvError};

        use rayon::prelude::*;

        use rand::{Rng, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        const RUNS_PER_WORK_UNIT: usize = 1000;
        const PARALLEL_WORKS: usize = 50;

        use ethereum_types::{U256};

        use indicatif::{ProgressBar, ProgressStyle};

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        let mut parameters_space = vec![];
        let (tx, rx) = channel();

        for i in 0..PARALLEL_WORKS {
            parameters_space.push((i, rng.clone(), pb.clone(), tx.clone()));
        }

        drop(tx);

        pb.set_length((parameters_space.len() * RUNS_PER_WORK_UNIT) as u64);

        let handler = thread::spawn(move || {
            parameters_space.into_par_iter().for_each(|(i, mut rng, pb, tx)| {
                for _ in 0..RUNS_PER_WORK_UNIT {
                    let mut r: u64 = rng.gen();
                    r = r.wrapping_add(i as u64);
                    let input = if r & 1 == 1 {
                        hex::decode("00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed275dc4a288d1afb3cbb1ac09187524c7db36395df7be3b99e673b13a075a65ec1d9befcd05a5323e6da4d435f3b617cdb3af83285c2df711ef39c01571827f9d00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed275dc4a288d1afb3cbb1ac09187524c7db36395df7be3b99e673b13a075a65ec1d9befcd05a5323e6da4d435f3b617cdb3af83285c2df711ef39c01571827f9d00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed275dc4a288d1afb3cbb1ac09187524c7db36395df7be3b99e673b13a075a65ec1d9befcd05a5323e6da4d435f3b617cdb3af83285c2df711ef39c01571827f9d00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed275dc4a288d1afb3cbb1ac09187524c7db36395df7be3b99e673b13a075a65ec1d9befcd05a5323e6da4d435f3b617cdb3af83285c2df711ef39c01571827f9d00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed275dc4a288d1afb3cbb1ac09187524c7db36395df7be3b99e673b13a075a65ec1d9befcd05a5323e6da4d435f3b617cdb3af83285c2df711ef39c01571827f9d").unwrap()
                    } else {
                        hex::decode("00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002203e205db4f19b37b60121b83a7333706db86431c6d835849957ed8c3928ad7927dc7234fd11d3e8c36c59277c3e6f149d5cd3cfa9a62aee49f8130962b4b3b9195e8aa5b7827463722b8c153931579d3505566b4edf48d498e185f0509de15204bb53b8977e5f92a0bc372742c4830944a59b4fe6b1c0466e2a6dad122b5d2e030644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd31a76dae6d3272396d0cbe61fced2bc532edac647851e3ac53ce1cc9c7e645a83198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002203e205db4f19b37b60121b83a7333706db86431c6d835849957ed8c3928ad7927dc7234fd11d3e8c36c59277c3e6f149d5cd3cfa9a62aee49f8130962b4b3b9195e8aa5b7827463722b8c153931579d3505566b4edf48d498e185f0509de15204bb53b8977e5f92a0bc372742c4830944a59b4fe6b1c0466e2a6dad122b5d2e030644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd31a76dae6d3272396d0cbe61fced2bc532edac647851e3ac53ce1cc9c7e645a83198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002203e205db4f19b37b60121b83a7333706db86431c6d835849957ed8c3928ad7927dc7234fd11d3e8c36c59277c3e6f149d5cd3cfa9a62aee49f8130962b4b3b9195e8aa5b7827463722b8c153931579d3505566b4edf48d498e185f0509de15204bb53b8977e5f92a0bc372742c4830944a59b4fe6b1c0466e2a6dad122b5d2e030644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd31a76dae6d3272396d0cbe61fced2bc532edac647851e3ac53ce1cc9c7e645a83198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002203e205db4f19b37b60121b83a7333706db86431c6d835849957ed8c3928ad7927dc7234fd11d3e8c36c59277c3e6f149d5cd3cfa9a62aee49f8130962b4b3b9195e8aa5b7827463722b8c153931579d3505566b4edf48d498e185f0509de15204bb53b8977e5f92a0bc372742c4830944a59b4fe6b1c0466e2a6dad122b5d2e030644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd31a76dae6d3272396d0cbe61fced2bc532edac647851e3ac53ce1cc9c7e645a83198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa00000000000000000000000000000000000000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000002203e205db4f19b37b60121b83a7333706db86431c6d835849957ed8c3928ad7927dc7234fd11d3e8c36c59277c3e6f149d5cd3cfa9a62aee49f8130962b4b3b9195e8aa5b7827463722b8c153931579d3505566b4edf48d498e185f0509de15204bb53b8977e5f92a0bc372742c4830944a59b4fe6b1c0466e2a6dad122b5d2e030644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd31a76dae6d3272396d0cbe61fced2bc532edac647851e3ac53ce1cc9c7e645a83198e9393920d483a7260bfb731fb5d25f1aa493335a9e71297e485b7aef312c21800deef121f1e76426a00665e5c4479674322d4f75edadd46debd5cd992f6ed090689d0585ff075ec9e99ad690c3395bc4b313370b38ef355acdadcd122975b12c85ea5db8c6deb4aab71808dcb408fe3d1e7690c43d37b4ce6cc0166fa7daa").unwrap()
                    };

                    use bn::{AffineG1, AffineG2, Fq, Fq2, pairing_batch, G1, G2, Gt, Group};

                    let start = std::time::Instant::now();

                    let ret_val = if input.is_empty() {
                        U256::one()
                    } else {
                        // (a, b_a, b_b - each 64-byte affine coordinates)
                        let elements = input.len() / 192;
                        let mut vals = Vec::new();
                        for idx in 0..elements {
                            let a_x = Fq::from_slice(&input[idx*192..idx*192+32]).unwrap();

                            let a_y = Fq::from_slice(&input[idx*192+32..idx*192+64]).unwrap();

                            let b_a_y = Fq::from_slice(&input[idx*192+64..idx*192+96]).unwrap();

                            let b_a_x = Fq::from_slice(&input[idx*192+96..idx*192+128]).unwrap();

                            let b_b_y = Fq::from_slice(&input[idx*192+128..idx*192+160]).unwrap();

                            let b_b_x = Fq::from_slice(&input[idx*192+160..idx*192+192]).unwrap();

                            let b_a = Fq2::new(b_a_x, b_a_y);
                            let b_b = Fq2::new(b_b_x, b_b_y);
                            let b = if b_a.is_zero() && b_b.is_zero() {
                                G2::zero()
                            } else {
                                G2::from(AffineG2::new(b_a, b_b).unwrap())
                            };
                            let a = if a_x.is_zero() && a_y.is_zero() {
                                G1::zero()
                            } else {
                                G1::from(AffineG1::new(a_x, a_y).unwrap())
                            };
                            vals.push((a, b));
                        };

                        let mul = pairing_batch(&vals);

                        if mul == Gt::one() {
                            U256::one()
                        } else {
                            U256::zero()
                        }
                    };

                    let mut buf = [0u8; 32];
                    ret_val.to_big_endian(&mut buf);
                    
                    let elapsed_nanos = start.elapsed().as_nanos();
                    tx.send((elapsed_nanos, buf)).unwrap();

                    pb.inc(1);
                }
            });
        });

        let mut sum = 0u128;

        loop {
            let subres = rx.try_recv();
            match subres {
                Ok(r) => {
                    sum += r.0;
                },
                Err(TryRecvError::Empty) => {
                    std::thread::sleep(std::time::Duration::from_millis(1000u64));
                },
                Err(TryRecvError::Disconnected) => {
                    handler.join().unwrap();
                    break;
                }
            }
        }

        pb.finish();

        let ns_per_run = sum / ((PARALLEL_WORKS * RUNS_PER_WORK_UNIT) as u128);

        const BN_PAIRING_10_POINTS_GAS: u128 = 45000 + 34000 * 10;

        let gas_per_second = 1_000_000 * 1_000 * BN_PAIRING_10_POINTS_GAS / ns_per_run;

        println!("MGAS per second on the current PC = {}", (gas_per_second as f64) / 1_000_000f64);
    }

    #[test]
    #[ignore]
    fn benchmark_sha256_precompile() {
        use std::thread;

        use std::sync::mpsc::{channel, TryRecvError};

        use rayon::prelude::*;

        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        const RUNS_PER_WORK_UNIT: usize = 10000;

        use parity_crypto::digest;

        use indicatif::{ProgressBar, ProgressStyle};

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        let mut parameters_space = vec![];
        let (tx, rx) = channel();

        let input_len = vec![0, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192];

        for i in input_len.into_iter() {
            parameters_space.push((i, rng.clone(), pb.clone(), tx.clone()));
        }

        drop(tx);

        pb.set_length((parameters_space.len() * RUNS_PER_WORK_UNIT) as u64);

        let handler = thread::spawn(move || {
            parameters_space.into_par_iter().for_each(|(i, mut rng, pb, tx)| {
                let mut input = vec![0u8; i];
                let mut output = [0u8; 32];
                for _ in 0..RUNS_PER_WORK_UNIT {

                    rng.fill_bytes(&mut input);

                    let start = std::time::Instant::now();
                    use std::io::Write;

                    let d = digest::sha256(&input);
                    (&mut output[..]).write(&d).unwrap();

                    let elapsed_nanos = start.elapsed().as_nanos();
                    tx.send((elapsed_nanos, i)).unwrap();

                    pb.inc(1);
                }
            });
        });

        let mut sums = std::collections::HashMap::new();

        loop {
            let subres = rx.try_recv();
            match subres {
                Ok((r, i)) => {
                    let entry = sums.entry(i).or_insert(0u128);
                    *entry += r as u128;
                },
                Err(TryRecvError::Empty) => {
                    std::thread::sleep(std::time::Duration::from_millis(1000u64));
                },
                Err(TryRecvError::Disconnected) => {
                    handler.join().unwrap();
                    break;
                }
            }
        }

        pb.finish();

        const ECRECOVER_GAS: u128 = 3000;

        const GAS_PER_SECOND: u128 = 35_000_000;

        println!("Using {} gas/second", GAS_PER_SECOND);

        for (k, v) in sums.into_iter() {
            let ns_average = (v as u128) / (RUNS_PER_WORK_UNIT as u128);
            let gas_average = GAS_PER_SECOND * ns_average / 1_000_000_000u128;
            println!("Hashed {} bytes for {} gas", k, gas_average);
        }
    }

    #[test]
    #[ignore]
    fn benchmark_ripemd160_precompile() {
        use std::thread;

        use std::sync::mpsc::{channel, TryRecvError};

        use rayon::prelude::*;

        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        const RUNS_PER_WORK_UNIT: usize = 10000;

        use parity_crypto::digest;

        use indicatif::{ProgressBar, ProgressStyle};

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        let mut parameters_space = vec![];
        let (tx, rx) = channel();

        let input_len = vec![0, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192];

        for i in input_len.into_iter() {
            parameters_space.push((i, rng.clone(), pb.clone(), tx.clone()));
        }

        drop(tx);

        pb.set_length((parameters_space.len() * RUNS_PER_WORK_UNIT) as u64);

        let handler = thread::spawn(move || {
            parameters_space.into_par_iter().for_each(|(i, mut rng, pb, tx)| {
                let mut input = vec![0u8; i];
                let mut output = [0u8; 20];
                for _ in 0..RUNS_PER_WORK_UNIT {

                    rng.fill_bytes(&mut input);

                    let start = std::time::Instant::now();
                    use std::io::Write;

                    let d = digest::ripemd160(&input);
                    (&mut output[..]).write(&d).unwrap();

                    let elapsed_nanos = start.elapsed().as_nanos();
                    tx.send((elapsed_nanos, i)).unwrap();

                    pb.inc(1);
                }
            });
        });

        let mut sums = std::collections::HashMap::new();

        loop {
            let subres = rx.try_recv();
            match subres {
                Ok((r, i)) => {
                    let entry = sums.entry(i).or_insert(0u128);
                    *entry += r as u128;
                },
                Err(TryRecvError::Empty) => {
                    std::thread::sleep(std::time::Duration::from_millis(1000u64));
                },
                Err(TryRecvError::Disconnected) => {
                    handler.join().unwrap();
                    break;
                }
            }
        }

        pb.finish();

        const ECRECOVER_GAS: u128 = 3000;

        const GAS_PER_SECOND: u128 = 35_000_000;

        println!("Using {} gas/second", GAS_PER_SECOND);

        for (k, v) in sums.into_iter() {
            let ns_average = (v as u128) / (RUNS_PER_WORK_UNIT as u128);
            let gas_average = GAS_PER_SECOND * ns_average / 1_000_000_000u128;
            println!("Hashed {} bytes for {} gas", k, gas_average);
        }
    }

    #[test]
    #[ignore]
    fn benchmark_blake2f_precompile() {
        use std::thread;

        use std::sync::mpsc::{channel, TryRecvError};

        use rayon::prelude::*;

        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        const RUNS_PER_WORK_UNIT: usize = 10000;

        use indicatif::{ProgressBar, ProgressStyle};

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        let mut parameters_space = vec![];
        let (tx, rx) = channel();

        let num_rounds: Vec<u32> = vec![0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192];

        for i in num_rounds.into_iter() {
            parameters_space.push((i, rng.clone(), pb.clone(), tx.clone()));
        }

        drop(tx);

        pb.set_length((parameters_space.len() * RUNS_PER_WORK_UNIT) as u64);

        let handler = thread::spawn(move || {
            parameters_space.into_par_iter().for_each(|(i, mut rng, pb, tx)| {
                use std::io::{Cursor, Write};
                use byteorder::{BigEndian, LittleEndian};
                use byteorder::{ReadBytesExt, WriteBytesExt};
                use eip_152::compress;

                const BLAKE2_F_ARG_LEN: usize = 213;
                const PROOF: &str = "Checked the length of the input above; qed";
        
                let mut input = vec![0u8; BLAKE2_F_ARG_LEN];
                (&mut input[0..4]).write_u32::<BigEndian>(i).expect("must write number of rounds");
                 
                let mut output = [0u8; 64];
                for _ in 0..RUNS_PER_WORK_UNIT {
                    rng.fill_bytes(&mut input[4..]);
                    let last_byte = input[BLAKE2_F_ARG_LEN-1];
                    input[BLAKE2_F_ARG_LEN-1] = last_byte & 1u8; // last byte is a binary flag

                    let start = std::time::Instant::now();

                    if input.len() != BLAKE2_F_ARG_LEN {
                        panic!("input length for Blake2 F precompile should be exactly 213 bytes");
                    }
            
                    let mut cursor = Cursor::new(&input);
                    let rounds = cursor.read_u32::<BigEndian>().expect(PROOF);
            
                    // state vector, h
                    let mut h = [0u64; 8];
                    for state_word in &mut h {
                        *state_word = cursor.read_u64::<LittleEndian>().expect(PROOF);
                    }
            
                    // message block vector, m
                    let mut m = [0u64; 16];
                    for msg_word in &mut m {
                        *msg_word = cursor.read_u64::<LittleEndian>().expect(PROOF);
                    }
            
                    // 2w-bit offset counter, t
                    let t = [
                        cursor.read_u64::<LittleEndian>().expect(PROOF),
                        cursor.read_u64::<LittleEndian>().expect(PROOF),
                    ];
            
                    // final block indicator flag, "f"
                    let f = match input.last() {
                            Some(1) => true,
                            Some(0) => false,
                            _ => {
                                panic!("incorrect final block indicator flag, was: {:?}", input.last());
                            }
                        };
            
                    compress(&mut h, m, t, f, rounds as usize);
            
                    let mut output_buf = [0u8; 8 * std::mem::size_of::<u64>()];
                    for (i, state_word) in h.iter().enumerate() {
                        output_buf[i*8..(i+1)*8].copy_from_slice(&state_word.to_le_bytes());
                    }

                    (&mut output[..]).write(&output_buf).unwrap();

                    let elapsed_nanos = start.elapsed().as_nanos();
                    tx.send((elapsed_nanos, i)).unwrap();

                    pb.inc(1);
                }
            });
        });

        let mut sums = std::collections::HashMap::new();

        loop {
            let subres = rx.try_recv();
            match subres {
                Ok((r, i)) => {
                    let entry = sums.entry(i).or_insert(0u128);
                    *entry += r as u128;
                },
                Err(TryRecvError::Empty) => {
                    std::thread::sleep(std::time::Duration::from_millis(1000u64));
                },
                Err(TryRecvError::Disconnected) => {
                    handler.join().unwrap();
                    break;
                }
            }
        }

        pb.finish();

        const ECRECOVER_GAS: u128 = 3000;

        const GAS_PER_SECOND: u128 = 35_000_000;

        println!("Using {} gas/second", GAS_PER_SECOND);

        for (k, v) in sums.into_iter() {
            let ns_average = (v as u128) / (RUNS_PER_WORK_UNIT as u128);
            let gas_average = GAS_PER_SECOND * ns_average / 1_000_000_000u128;
            println!("Hashed {} rounds for {} gas", k, gas_average);
        }
    }

    fn read_fr(reader: &[u8]) -> Result<bn::Fr, &'static str> {
        let mut buf = [0u8; 32];
        buf.copy_from_slice(&reader);

        bn::Fr::from_slice(&buf[0..32]).map_err(|_| "Invalid field element")
    }
    
    fn read_point(reader: &[u8]) -> Result<bn::G1, &'static str> {
        use bn::{Fq, AffineG1, G1, Group};
    
        let mut buf = [0u8; 32];

        buf.copy_from_slice(&reader[0..32]);
    
        let px = Fq::from_slice(&buf[0..32]).map_err(|_| "Invalid point x coordinate")?;

        buf.copy_from_slice(&reader[32..64]);

        let py = Fq::from_slice(&buf[0..32]).map_err(|_| "Invalid point y coordinate")?;
        Ok(
            if px == Fq::zero() && py == Fq::zero() {
                G1::zero()
            } else {
                AffineG1::new(px, py).map_err(|_| "Invalid curve point")?.into()
            }
        )
    }

    #[test]
    #[ignore]
    fn benchmark_bn_add_precompile() {
        use std::thread;

        use std::sync::mpsc::{channel, TryRecvError};

        use rayon::prelude::*;

        use rand::{Rng, SeedableRng, RngCore};
        use rand_xorshift::XorShiftRng;

        let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        const RUNS_PER_WORK_UNIT: usize = 1000;
        const PARALLEL_WORKS: usize = 50;

        use indicatif::{ProgressBar, ProgressStyle};

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        let mut parameters_space = vec![];
        let (tx, rx) = channel();

        for i in 0..PARALLEL_WORKS {
            let _: u8 = rng.gen();
            parameters_space.push((i, rng.clone(), pb.clone(), tx.clone()));
        }

        drop(tx);

        pb.set_length((parameters_space.len() * RUNS_PER_WORK_UNIT) as u64);

        let handler = thread::spawn(move || {
            parameters_space.into_par_iter().for_each(|(_i, mut rng, pb, tx)| {
                let mut scalar_buffer = vec![0u8; 32];
                let mut base_point_x = vec![0u8; 32];
                base_point_x[31] = 1;
                let mut base_point_y = vec![0u8; 32];
                base_point_y[31] = 2;
                let mut base_point_encoding = vec![];
                base_point_encoding.extend(base_point_x);
                base_point_encoding.extend(base_point_y);

                let mut input = vec![0u8; 128];

                let base_point = read_point(&mut base_point_encoding).unwrap();

                let mut output = vec![0u8; 32];

                for _ in 0..RUNS_PER_WORK_UNIT {
                    rng.fill_bytes(&mut scalar_buffer);
                    scalar_buffer[0] = 0;

                    let fr = read_fr(&mut scalar_buffer).unwrap();

                    let p1 = AffineG1::from_jacobian(base_point * fr).unwrap();
                    p1.x().to_big_endian(&mut input[0..32]).expect("Cannot fail since 0..32 is 32-byte length");
                    p1.y().to_big_endian(&mut input[32..64]).expect("Cannot fail since 32..64 is 32-byte length");

                    rng.fill_bytes(&mut scalar_buffer);
                    scalar_buffer[0] = 0;

                    let fr = read_fr(&mut scalar_buffer).unwrap();

                    let p2 = AffineG1::from_jacobian(base_point * fr).unwrap();
                    p2.x().to_big_endian(&mut input[64..96]).expect("Cannot fail since 0..32 is 32-byte length");
                    p2.y().to_big_endian(&mut input[96..128]).expect("Cannot fail since 32..64 is 32-byte length");


                    use bn::AffineG1;

                    let start = std::time::Instant::now();

                    let p1 = read_point(&mut input[0..64]).unwrap();
                    let p2 = read_point(&mut input[64..128]).unwrap();
            
                    let mut write_buf = [0u8; 64];
                    if let Some(sum) = AffineG1::from_jacobian(p1 + p2) {
                        // point not at infinity
                        sum.x().to_big_endian(&mut write_buf[0..32]).expect("Cannot fail since 0..32 is 32-byte length");
                        sum.y().to_big_endian(&mut write_buf[32..64]).expect("Cannot fail since 32..64 is 32-byte length");
                    }
                    use std::io::Write;

                    (&mut output).write(&write_buf).unwrap();
                    
                    let elapsed_nanos = start.elapsed().as_nanos();
                    tx.send(elapsed_nanos).unwrap();

                    pb.inc(1);
                }
            });
        });

        let mut sum = 0u128;

        loop {
            let subres = rx.try_recv();
            match subres {
                Ok(r) => {
                    sum += r;
                },
                Err(TryRecvError::Empty) => {
                    std::thread::sleep(std::time::Duration::from_millis(1000u64));
                },
                Err(TryRecvError::Disconnected) => {
                    handler.join().unwrap();
                    break;
                }
            }
        }

        pb.finish();

        const GAS_PER_SECOND: u128 = 35_000_000;

        println!("Using {} gas/second", GAS_PER_SECOND);

        let ns_average = (sum as u128) / (RUNS_PER_WORK_UNIT as u128) / (PARALLEL_WORKS as u128);
        let gas_average = GAS_PER_SECOND * ns_average / 1_000_000_000u128;
        println!("BN_ADD used {} gas on average", gas_average);
    }

    #[test]
    #[ignore]
    fn benchmark_bn_mul_precompile() {
        use std::thread;

        use std::sync::mpsc::{channel, TryRecvError};

        use rayon::prelude::*;

        use rand::{Rng, SeedableRng, RngCore};
        use rand_xorshift::XorShiftRng;

        let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        const RUNS_PER_WORK_UNIT: usize = 1000;
        const PARALLEL_WORKS: usize = 50;

        use indicatif::{ProgressBar, ProgressStyle};

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        let mut parameters_space = vec![];
        let (tx, rx) = channel();

        for i in 0..PARALLEL_WORKS {
            let _: u8 = rng.gen();
            parameters_space.push((i, rng.clone(), pb.clone(), tx.clone()));
        }

        drop(tx);

        pb.set_length((parameters_space.len() * RUNS_PER_WORK_UNIT) as u64);

        let handler = thread::spawn(move || {
            parameters_space.into_par_iter().for_each(|(_i, mut rng, pb, tx)| {
                let mut scalar_buffer = vec![0u8; 32];
                let mut base_point_x = vec![0u8; 32];
                base_point_x[31] = 1;
                let mut base_point_y = vec![0u8; 32];
                base_point_y[31] = 2;
                let mut base_point_encoding = vec![];
                base_point_encoding.extend(base_point_x);
                base_point_encoding.extend(base_point_y);

                let mut input = vec![0u8; 64];

                let worst_case_scalar = vec![255u8; 32];

                input.extend(worst_case_scalar);

                let base_point = read_point(&mut base_point_encoding).unwrap();

                let mut output = vec![0u8; 32];

                for _ in 0..RUNS_PER_WORK_UNIT {
                    rng.fill_bytes(&mut scalar_buffer);
                    scalar_buffer[0] = 0;

                    let fr = read_fr(&mut scalar_buffer).unwrap();

                    let p1 = AffineG1::from_jacobian(base_point * fr).unwrap();
                    p1.x().to_big_endian(&mut input[0..32]).expect("Cannot fail since 0..32 is 32-byte length");
                    p1.y().to_big_endian(&mut input[32..64]).expect("Cannot fail since 32..64 is 32-byte length");

                    use bn::AffineG1;

                    let start = std::time::Instant::now();

                    let p = read_point(&mut input[0..64]).unwrap();
                    let fr = read_fr(&mut input[64..96]).unwrap();
            
                    let mut write_buf = [0u8; 64];
                    if let Some(sum) = AffineG1::from_jacobian(p * fr) {
                        // point not at infinity
                        sum.x().to_big_endian(&mut write_buf[0..32]).expect("Cannot fail since 0..32 is 32-byte length");
                        sum.y().to_big_endian(&mut write_buf[32..64]).expect("Cannot fail since 32..64 is 32-byte length");
                    }
                    use std::io::Write;

                    (&mut output).write(&write_buf).unwrap();
                    
                    let elapsed_nanos = start.elapsed().as_nanos();
                    tx.send(elapsed_nanos).unwrap();

                    pb.inc(1);
                }
            });
        });

        let mut sum = 0u128;

        loop {
            let subres = rx.try_recv();
            match subres {
                Ok(r) => {
                    sum += r;
                },
                Err(TryRecvError::Empty) => {
                    std::thread::sleep(std::time::Duration::from_millis(1000u64));
                },
                Err(TryRecvError::Disconnected) => {
                    handler.join().unwrap();
                    break;
                }
            }
        }

        pb.finish();

        const GAS_PER_SECOND: u128 = 35_000_000;

        println!("Using {} gas/second", GAS_PER_SECOND);

        let ns_average = (sum as u128) / (RUNS_PER_WORK_UNIT as u128) / (PARALLEL_WORKS as u128);
        let gas_average = GAS_PER_SECOND * ns_average / 1_000_000_000u128;
        println!("BN_ADD used {} gas on average", gas_average);
    }

    #[test]
    #[ignore]
    fn benchmark_bn_pairing_precompile() {
        use std::thread;

        use std::sync::mpsc::{channel, TryRecvError};

        use rayon::prelude::*;

        use rand::{Rng, RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let mut rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        const RUNS_PER_WORK_UNIT: usize = 10000;

        use ethereum_types::U256;

        use indicatif::{ProgressBar, ProgressStyle};

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        let mut parameters_space = vec![];
        let (tx, rx) = channel();

        let num_pairs = vec![1, 2, 3, 4, 6, 8, 16, 32];

        for i in num_pairs.into_iter() {
            let _ : u8 = rng.gen();
            parameters_space.push((i, rng.clone(), pb.clone(), tx.clone()));
        }

        drop(tx);

        pb.set_length((parameters_space.len() * RUNS_PER_WORK_UNIT) as u64);

        let handler = thread::spawn(move || {
            parameters_space.into_par_iter().for_each(|(i, mut rng, pb, tx)| {
                use bn::{AffineG1, AffineG2, G1, G2, Group, Fq, Fq2, pairing_batch, Gt};

                let mut input = vec![0u8; i * (64 + 128)];
                let mut output = [0u8; 32];

                let mut scalar_buffer = vec![0u8; 32];

                use std::ops::Mul;

                for _ in 0..RUNS_PER_WORK_UNIT {
                    let mut offset = 0;
                    for _ in 0..i {
                        rng.fill_bytes(&mut scalar_buffer);
                        scalar_buffer[0] = 0;
    
                        let fr = read_fr(&mut scalar_buffer).unwrap();
                        
                        let p1 = bn::G1::one().mul(fr);
                        let p1 = AffineG1::from_jacobian(p1).unwrap();
                        p1.x().to_big_endian(&mut input[offset..(offset+32)]).expect("Cannot fail since 0..32 is 32-byte length");
                        offset += 32;
                        p1.y().to_big_endian(&mut input[offset..(offset+32)]).expect("Cannot fail since 32..64 is 32-byte length");
                        offset += 32;

                        rng.fill_bytes(&mut scalar_buffer);
                        scalar_buffer[0] = 0;

                        let fr = read_fr(&mut scalar_buffer).unwrap();

                        let p2 = bn::G2::one().mul(fr);
                        let p2 = AffineG2::from_jacobian(p2).unwrap();
                        p2.x().imaginary().to_big_endian(&mut input[offset..(offset+32)]).expect("Cannot fail since 0..32 is 32-byte length");
                        offset += 32;
                        p2.x().real().to_big_endian(&mut input[offset..(offset+32)]).expect("Cannot fail since 0..32 is 32-byte length");
                        offset += 32;

                        p2.y().imaginary().to_big_endian(&mut input[offset..(offset+32)]).expect("Cannot fail since 0..32 is 32-byte length");
                        offset += 32;
                        p2.y().real().to_big_endian(&mut input[offset..(offset+32)]).expect("Cannot fail since 0..32 is 32-byte length");
                        offset += 32;
                    }

                    let start = std::time::Instant::now();
                    use std::io::Write;

                    let ret_val = if input.is_empty() {
                        U256::one()
                    } else {
                        // (a, b_a, b_b - each 64-byte affine coordinates)
                        let elements = input.len() / 192;
                        let mut vals = Vec::new();
                        for idx in 0..elements {
                            let a_x = Fq::from_slice(&input[idx*192..idx*192+32])
                                .map_err(|_| "Invalid a argument x coordinate").unwrap();
            
                            let a_y = Fq::from_slice(&input[idx*192+32..idx*192+64])
                                .map_err(|_| "Invalid a argument y coordinate").unwrap();
            
                            let b_a_y = Fq::from_slice(&input[idx*192+64..idx*192+96])
                                .map_err(|_| "Invalid b argument imaginary coeff x coordinate").unwrap();
            
                            let b_a_x = Fq::from_slice(&input[idx*192+96..idx*192+128])
                                .map_err(|_| "Invalid b argument imaginary coeff y coordinate").unwrap();
            
                            let b_b_y = Fq::from_slice(&input[idx*192+128..idx*192+160])
                                .map_err(|_| "Invalid b argument real coeff x coordinate").unwrap();
            
                            let b_b_x = Fq::from_slice(&input[idx*192+160..idx*192+192])
                                .map_err(|_| "Invalid b argument real coeff y coordinate").unwrap();
            
                            let b_a = Fq2::new(b_a_x, b_a_y);
                            let b_b = Fq2::new(b_b_x, b_b_y);
                            let b = if b_a.is_zero() && b_b.is_zero() {
                                G2::zero()
                            } else {
                                G2::from(AffineG2::new(b_a, b_b).map_err(|_| "Invalid b argument - not on curve").unwrap())
                            };
                            let a = if a_x.is_zero() && a_y.is_zero() {
                                G1::zero()
                            } else {
                                G1::from(AffineG1::new(a_x, a_y).map_err(|_| "Invalid a argument - not on curve").unwrap())
                            };
                            vals.push((a, b));
                        };
            
                        let mul = pairing_batch(&vals);
            
                        if mul == Gt::one() {
                            U256::one()
                        } else {
                            U256::zero()
                        }
                    };
            
                    let mut buf = [0u8; 32];
                    ret_val.to_big_endian(&mut buf);
                    (&mut output[..]).write(&buf).unwrap();

                    let elapsed_nanos = start.elapsed().as_nanos();
                    tx.send((elapsed_nanos, i)).unwrap();

                    pb.inc(1);
                }
            });
        });

        let mut sums = std::collections::HashMap::new();

        loop {
            let subres = rx.try_recv();
            match subres {
                Ok((r, i)) => {
                    let entry = sums.entry(i).or_insert(0u128);
                    *entry += r as u128;
                },
                Err(TryRecvError::Empty) => {
                    std::thread::sleep(std::time::Duration::from_millis(1000u64));
                },
                Err(TryRecvError::Disconnected) => {
                    handler.join().unwrap();
                    break;
                }
            }
        }

        pb.finish();

        const ECRECOVER_GAS: u128 = 3000;

        const GAS_PER_SECOND: u128 = 35_000_000;

        println!("Using {} gas/second", GAS_PER_SECOND);

        for (k, v) in sums.into_iter() {
            let ns_average = (v as u128) / (RUNS_PER_WORK_UNIT as u128);
            let gas_average = GAS_PER_SECOND * ns_average / 1_000_000_000u128;
            println!("Paired {} pairs for {} gas", k, gas_average);
        }
    }

    #[test]
    #[ignore]
    fn benchmark_keccak_function() {
        use std::thread;

        use std::sync::mpsc::{channel, TryRecvError};

        use rayon::prelude::*;

        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        const RUNS_PER_WORK_UNIT: usize = 1000000;

        use keccak_hash;

        use indicatif::{ProgressBar, ProgressStyle};

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        let mut parameters_space = vec![];
        let (tx, rx) = channel();

        let input_len = vec![0, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192];

        for i in input_len.into_iter() {
            parameters_space.push((i, rng.clone(), pb.clone(), tx.clone()));
        }

        drop(tx);

        pb.set_length((parameters_space.len() * RUNS_PER_WORK_UNIT) as u64);

        let handler = thread::spawn(move || {
            parameters_space.into_par_iter().for_each(|(i, mut rng, pb, tx)| {
                let mut input = vec![0u8; i];
                let mut output = [0u8; 32];
                for _ in 0..RUNS_PER_WORK_UNIT {

                    rng.fill_bytes(&mut input);

                    let start = std::time::Instant::now();
                    use std::io::Write;

                    let d = keccak_hash::keccak(&input);
                    (&mut output[..]).write(&d.as_bytes()).unwrap();

                    let elapsed_nanos = start.elapsed().as_nanos();
                    tx.send((elapsed_nanos, i)).unwrap();

                    pb.inc(1);
                }
            });
        });

        let mut sums = std::collections::HashMap::new();

        loop {
            let subres = rx.try_recv();
            match subres {
                Ok((r, i)) => {
                    let entry = sums.entry(i).or_insert(0u128);
                    *entry += r as u128;
                },
                Err(TryRecvError::Empty) => {
                    std::thread::sleep(std::time::Duration::from_millis(1000u64));
                },
                Err(TryRecvError::Disconnected) => {
                    handler.join().unwrap();
                    break;
                }
            }
        }

        pb.finish();

        const ECRECOVER_GAS: u128 = 3000;

        const GAS_PER_SECOND: u128 = 35_000_000;

        println!("Using {} gas/second", GAS_PER_SECOND);

        for (k, v) in sums.into_iter() {
            let ns_average = (v as u128) / (RUNS_PER_WORK_UNIT as u128);
            let gas_average = GAS_PER_SECOND * ns_average / 1_000_000_000u128;
            println!("Hashed {} bytes for {} gas", k, gas_average);
        }
    }

    #[test]
    #[ignore]
    fn benchmark_keccak_sponge_price() {
        use std::thread;

        use std::sync::mpsc::{channel, TryRecvError};

        use rayon::prelude::*;

        use rand::{RngCore, SeedableRng};
        use rand_xorshift::XorShiftRng;

        let rng = XorShiftRng::from_seed([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]);

        const RUNS_PER_WORK_UNIT: usize = 10000;

        use keccak_hash;

        use indicatif::{ProgressBar, ProgressStyle};

        let pb = ProgressBar::new(1u64);

        pb.set_style(ProgressStyle::default_bar()
            .template("[{elapsed_precise}|{eta_precise}] {bar:50} {pos:>7}/{len:7} {msg}")
            .progress_chars("##-"));

        let mut parameters_space = vec![];
        let (tx, rx) = channel();

        let multiple = 8;

        for i in 0..=(768/multiple) {
            parameters_space.push((i*multiple, rng.clone(), pb.clone(), tx.clone()));
        }

        drop(tx);

        pb.set_length((parameters_space.len() * RUNS_PER_WORK_UNIT) as u64);

        let handler = thread::spawn(move || {
            parameters_space.into_par_iter().for_each(|(i, mut rng, pb, tx)| {
                let mut input = vec![0u8; i];
                let mut output = [0u8; 32];
                for _ in 0..RUNS_PER_WORK_UNIT {

                    rng.fill_bytes(&mut input);

                    let start = std::time::Instant::now();
                    use std::io::Write;

                    let d = keccak_hash::keccak(&input);
                    (&mut output[..]).write(&d.as_bytes()).unwrap();

                    let elapsed_nanos = start.elapsed().as_nanos();
                    tx.send((elapsed_nanos, i)).unwrap();

                    pb.inc(1);
                }
            });
        });

        let mut sums = std::collections::HashMap::new();

        loop {
            let subres = rx.try_recv();
            match subres {
                Ok((r, i)) => {
                    let entry = sums.entry(i).or_insert(0u128);
                    *entry += r as u128;
                },
                Err(TryRecvError::Empty) => {
                    std::thread::sleep(std::time::Duration::from_millis(1000u64));
                },
                Err(TryRecvError::Disconnected) => {
                    handler.join().unwrap();
                    break;
                }
            }
        }

        pb.finish();

        const ECRECOVER_GAS: u128 = 3000;

        const GAS_PER_SECOND: u128 = 35_000_000;

        println!("Using {} gas/second", GAS_PER_SECOND);

        let mut as_vec: Vec<_> = sums.into_iter().collect();
        as_vec.sort_by(|a, b| a.0.cmp(&b.0));

        for (k, v) in as_vec.into_iter() {
            let ns_average = (v as u128) / (RUNS_PER_WORK_UNIT as u128);
            let gas_average = GAS_PER_SECOND * ns_average / 1_000_000_000u128;
            println!("Hashed {} bytes for {} gas", k, gas_average);
        }
    }
}
