#[test]
fn test_fp2_inversion() {
    use num_bigint::BigUint;
    use crate::field::{U320Repr, new_field};
    use crate::fp::Fp;
    use crate::traits::{FieldElement, ZeroAndOne};
    use crate::extension_towers::fp2::{Fp2, Extension2};
    // use crate::extension_towers::fp4_as_2_over_2::{Extension2Over2};
    use num_traits::Num;

    let modulus = BigUint::from_str_radix("475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081", 10).unwrap();
    let base_field = new_field::<U320Repr>("475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081", 10).unwrap();
    let nonres_repr = U320Repr::from(17);
    let fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

    let mut extension_2 = Extension2::new(fp_non_residue.clone());
    extension_2.calculate_frobenius_coeffs(modulus.clone()).expect("must work");

    let mut fp2 = Fp2::one(&extension_2);
    fp2.c1 = fp_non_residue;

    let inverse = fp2.inverse().unwrap();
    let mut maybe_one = fp2.clone();
    maybe_one.mul_assign(&inverse);

    assert!(maybe_one == Fp2::one(&extension_2));
}

#[test]
fn test_fp4_inversion() {
    use num_bigint::BigUint;
    use crate::field::{U320Repr, new_field};
    use crate::fp::Fp;
    use crate::traits::{FieldElement, ZeroAndOne};
    use crate::extension_towers::fp2::{Fp2, Extension2};
    use crate::extension_towers::fp4_as_2_over_2::{Fp4, Extension2Over2};
    use num_traits::Num;

    let modulus = BigUint::from_str_radix("475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081", 10).unwrap();
    let base_field = new_field::<U320Repr>("475922286169261325753349249653048451545124879242694725395555128576210262817955800483758081", 10).unwrap();
    let nonres_repr = U320Repr::from(17);
    let fp_non_residue = Fp::from_repr(&base_field, nonres_repr).unwrap();

    let mut extension_2 = Extension2::new(fp_non_residue.clone());
    extension_2.calculate_frobenius_coeffs(modulus.clone()).expect("must work");

    let mut fp2_non_residue = Fp2::zero(&extension_2); // non-residue is 13 + 0*u + 0*u^2
    fp2_non_residue.c0 = fp_non_residue.clone();

    let mut extension_4 = Extension2Over2::new(fp2_non_residue);
    extension_4.calculate_frobenius_coeffs(modulus.clone()).expect("must work");

    let mut fp2 = Fp2::one(&extension_2);
    fp2.c1 = fp_non_residue;

    let mut fp4 = Fp4::one(&extension_4);
    fp4.c1 = fp2;

    let inverse = fp4.inverse().unwrap();
    let mut maybe_one = fp4.clone();
    maybe_one.mul_assign(&inverse);

    assert!(maybe_one == Fp4::one(&extension_4));
}