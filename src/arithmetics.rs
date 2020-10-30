/// Calculate a - b - borrow, returning the result and modifying
/// the borrow value.
#[inline(always)]
pub fn sbb(a: u64, b: u64, borrow: &mut u64) -> u64 {
    let tmp = (1u128 << 64) + u128::from(a) - u128::from(b) - u128::from(*borrow);

    *borrow = if tmp >> 64 == 0 { 1 } else { 0 };

    tmp as u64
}

/// Calculate a + b + carry, returning the sum and modifying the
/// carry value.
#[inline(always)]
pub fn adc(a: u64, b: u64, carry: &mut u64) -> u64 {
    let tmp = u128::from(a) + u128::from(b) + u128::from(*carry);

    *carry = (tmp >> 64) as u64;

    tmp as u64
}

/// Calculate a + (b * c) + carry, returning the least significant digit
/// and setting carry to the most significant digit.
#[inline(always)]
pub fn mac_with_carry(a: u64, b: u64, c: u64, carry: &mut u64) -> u64 {
    let tmp = (u128::from(a)) + u128::from(b) * u128::from(c) + u128::from(*carry);

    *carry = (tmp >> 64) as u64;

    tmp as u64
}

#[cfg(all(target_feature = "bmi2", feature = "asm"))]
#[inline(always)]
pub fn mulx(a: u64, b: u64) -> (u64, u64) {
    let mut lo: u64 = a;
    let mut hi: u64 = b;

    unsafe {
                asm!(
                    "mov rdx, {lo}", // place a in rdx
                    "mulx {hi}, {lo}, {hi}", // a*b -> (hi, lo)
                    lo = inlateout(reg) lo,
                    hi = inlateout(reg) hi,
                    out("rdx") _,
                    options(pure, readonly, nostack)
                );
            }

    (hi, lo)
}

#[cfg(any(not(target_feature = "bmi2"), not(feature = "asm")))]
#[inline(always)]
pub fn mulx(a: u64, b: u64) -> (u64, u64) {
    let tmp = (a as u128) * (b as u128);
    let hi = (tmp >> 64) as u64;
    let lo = tmp as u64;

    (hi, lo)
}