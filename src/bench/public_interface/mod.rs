extern crate test as rust_test;
use self::rust_test::Bencher;

use crate::public_interface::PairingApi;

#[bench]
fn bench_bls12_381_pairing_4_through_the_api(b: &mut Bencher) {
    use crate::test::pairings::bls12::assemble_single;

    let calldata = assemble_single(4);

    b.iter(|| {
        assert_eq!(crate::public_interface::PublicPairingApi::pair(&calldata).unwrap()[0], 1u8);
    });
}

#[bench]
fn bench_bls12_377_pairing_4_through_the_api(b: &mut Bencher) {
    use crate::test::pairings::bls12::assemble_single_bls12_377;

    let calldata = assemble_single_bls12_377(4);

    b.iter(|| {
        assert_eq!(crate::public_interface::PublicPairingApi::pair(&calldata).unwrap()[0], 1u8);
    });
}

#[bench]
fn bench_bn254_pairing_4_through_the_api(b: &mut Bencher) {
    use crate::test::pairings::bn::assemble_bn254;

    let calldata = assemble_bn254(4);

    b.iter(|| {
        assert_eq!(crate::public_interface::PublicPairingApi::pair(&calldata).unwrap()[0], 1u8);
    });
}

#[bench]
fn bench_mnt4_753_pairing_4_through_the_api(b: &mut Bencher) {
    use crate::test::pairings::mnt4::assemble_mnt4_753;

    let calldata = assemble_mnt4_753(4);

    b.iter(|| {
        // crate::public_interface::PublicPairingApi::pair(&calldata).unwrap();
        assert_eq!(crate::public_interface::PublicPairingApi::pair(&calldata).unwrap()[0], 1u8);
    });
}

#[bench]
fn bench_mnt4_753_pairing_1_through_the_api(b: &mut Bencher) {
    use crate::test::pairings::mnt4::assemble_mnt4_753;

    let calldata = assemble_mnt4_753(1);

    b.iter(|| {
        // crate::public_interface::PublicPairingApi::pair(&calldata).unwrap();
        assert_eq!(crate::public_interface::PublicPairingApi::pair(&calldata).unwrap()[0], 0u8);
    });
}

#[bench]
fn bench_bls12_381_pairing_1_through_the_api(b: &mut Bencher) {
    use crate::test::pairings::bls12::assemble_single;

    let calldata = assemble_single(1);

    b.iter(|| {
        assert_eq!(crate::public_interface::PublicPairingApi::pair(&calldata).unwrap()[0], 0u8);
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