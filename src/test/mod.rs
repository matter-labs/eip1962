pub(crate) mod pairings;
pub(crate) mod g2_ops;
pub(crate) mod g1_ops;
pub(crate) mod parsers;

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

}
