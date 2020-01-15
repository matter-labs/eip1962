extern crate test as rust_test;
use self::rust_test::Bencher;

use crate::public_interface::PairingApi;

#[bench]
fn bench_bls12_381_pairing_4_through_the_api(b: &mut Bencher) {
    use crate::test::pairings::bls12::assemble_single;

    let calldata = assemble_single(4);

    b.iter(|| {
        crate::public_interface::PublicPairingApi::pair(&calldata).unwrap();
    });
}

#[bench]
fn bench_bls12_377_pairing_4_through_the_api(b: &mut Bencher) {
    use crate::test::pairings::bls12::assemble_single_bls12_377;

    let calldata = assemble_single_bls12_377(4);

    b.iter(|| {
        crate::public_interface::PublicPairingApi::pair(&calldata).unwrap();
    });
}

#[bench]
fn bench_bls12_381_pairing_1_through_the_api(b: &mut Bencher) {
    use crate::test::pairings::bls12::assemble_single;

    let calldata = assemble_single(1);

    b.iter(|| {
        crate::public_interface::PublicPairingApi::pair(&calldata).unwrap();
    });
}

#[bench]
fn bench_bls12_381_pairing_through_the_api_no_actual_pairing(b: &mut Bencher) {
    use crate::test::pairings::bls12::assemble_single;

    let calldata = assemble_single(0);

    b.iter(|| {
        crate::public_interface::PublicPairingApi::pair(&calldata).expect_err("must fail for 0 pairs");
    });
}