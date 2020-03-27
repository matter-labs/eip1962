use crate::fp::*;
use crate::representation::*;
use crate::field::*;
use crate::traits::*;
use crate::extension_towers::fp2::Fp2;
use super::*;

use crate::weierstrass::*;
use crate::weierstrass::curve::*;
use crate::square_root::*;

#[derive(Clone)]
pub struct SwuParameters<C: CurveParameters> {
    pub z: C::BaseFieldElement,
    pub minus_b_by_a: C::BaseFieldElement,
    pub minus_z_inv: C::BaseFieldElement
}

pub(crate) fn simplified_swu_fp<
    'a,
    E: ElementRepr, 
    F: SizedPrimeField<Repr = E>, 
    C: CurveParameters<BaseFieldElement = Fp<'a, E, F>>
> (
    u: &C::BaseFieldElement,
    params: &SwuParameters<C>,
    curve: &WeierstrassCurve<'a, C>
) -> (C::BaseFieldElement, C::BaseFieldElement) {
    let one = Fp::one(u.field);

    // we do NOT use constant time operations here

    // 1.  tv1 = Z * u^2
    let mut tv1 = u.clone();
    tv1.square();
    tv1.mul_assign(&params.z);

    // 2.  tv2 = tv1^2
    let mut tv2 = tv1.clone();
    tv2.square();

    // 3.   x1 = tv1 + tv2
    let mut x1 = tv1.clone();
    x1.add_assign(&tv2);

    // 4.   x1 = inv0(x1)
    let mut x1 = x1.inverse().unwrap_or(Fp::zero(&u.field));

    // 5.   e1 = x1 == 0
    let e1 = x1.is_zero();

    // 6.   x1 = x1 + 1
    x1.add_assign(&one);

    // 7.   x1 = CMOV(x1, c2, e1)    # If (tv1 + tv2) == 0, set x1 = -1 / Z
    if e1 {
        x1 = params.minus_z_inv.clone();
    }

    // 8.   x1 = x1 * c1      # x1 = (-B / A) * (1 + (1 / (Z^2 * u^4 + Z * u^2)))
    x1.mul_assign(&params.minus_b_by_a);

    // 9.  gx1 = x1^2
    let mut gx1 = x1.clone();
    gx1.square();

    // 10. gx1 = gx1 + A
    gx1.add_assign(&curve.a);

    // 11. gx1 = gx1 * x1
    gx1.mul_assign(&x1);

    // 12. gx1 = gx1 + B             # gx1 = g(x1) = x1^3 + A * x1 + B
    gx1.add_assign(&curve.b);

    // 13.  x2 = tv1 * x1            # x2 = Z * u^2 * x1
    let mut x2 = tv1.clone();
    x2.mul_assign(&x1);

    // 14. tv2 = tv1 * tv2
    tv2.mul_assign(&tv1);

    // 15. gx2 = gx1 * tv2           # gx2 = (Z * u^2)^3 * gx1
    let mut gx2 = gx1.clone();
    gx2.mul_assign(&tv2);

    // 16.  e2 = is_square(gx1)
    let e2 = legendre_symbol_fp(&gx1);
    let e2 = e2 == LegendreSymbol::Zero || e2 == LegendreSymbol::QuadraticResidue;

    // 17.   x = CMOV(x2, x1, e2)    # If is_square(gx1), x = x1, else x = x2
    let x = if e2 {
        x1
    } else {
        x2
    };

    // 18.  y2 = CMOV(gx2, gx1, e2)  # If is_square(gx1), y2 = gx1, else y2 = gx2
    let y2 = if e2 {
        gx1
    } else {
        gx2
    };
    // 19.   y = sqrt(y2)
    let mut y = sqrt(&y2).expect("y2 is a square");

    // 20.  e3 = sgn0(u) == sgn0(y)  # Fix sign of y
    let u_sign = sign_of_fp_be(&u);
    let y_sign = sign_of_fp_be(&y);

    // 21.   y = CMOV(-y, y, e3)

    // our SignZero is equal to SignPlus
    match (u_sign, y_sign) {
        (Sign::Zero, Sign::Zero) |
        (Sign::Zero, Sign::SignPlus) | 
        (Sign::SignPlus, Sign::Zero) => {},
        (Sign::SignPlus, Sign::SignPlus) |
        (Sign::SignMinus, Sign::SignMinus) => {},
        (Sign::SignPlus, Sign::SignMinus) |
        (Sign::SignMinus, Sign::SignPlus) => {
            y.negate();
        },
        (Sign::Zero, Sign::SignMinus) | 
        (Sign::SignMinus, Sign::Zero) => {
            y.negate();
        },
    }

    // 22. return (x, y)
    (x, y)
}

pub(crate) fn simplified_swu_fp2<
    'a, 
    E: ElementRepr, 
    F: SizedPrimeField<Repr = E>, 
    C: CurveParameters<BaseFieldElement = Fp2<'a, E, F>>
> (
    u: &Fp2<'a, E, F>,
    params: &SwuParameters<C>,
    curve: &WeierstrassCurve<'a, C>
) -> (Fp2<'a, E, F>, Fp2<'a, E, F>)  {
    let one = Fp2::one(u.extension_field);

    // we do NOT use constant time operations here

    // 1.  tv1 = Z * u^2
    let mut tv1 = u.clone();
    tv1.square();
    tv1.mul_assign(&params.z);

    // 2.  tv2 = tv1^2
    let mut tv2 = tv1.clone();
    tv2.square();

    // 3.   x1 = tv1 + tv2
    let mut x1 = tv1.clone();
    x1.add_assign(&tv2);

    // 4.   x1 = inv0(x1)
    let mut x1 = x1.inverse().unwrap_or(Fp2::zero(&u.extension_field));

    // 5.   e1 = x1 == 0
    let e1 = x1.is_zero();

    // 6.   x1 = x1 + 1
    x1.add_assign(&one);

    // 7.   x1 = CMOV(x1, c2, e1)    # If (tv1 + tv2) == 0, set x1 = -1 / Z
    if e1 {
        x1 = params.minus_z_inv.clone();
    }

    // 8.   x1 = x1 * c1      # x1 = (-B / A) * (1 + (1 / (Z^2 * u^4 + Z * u^2)))
    x1.mul_assign(&params.minus_b_by_a);

    // 9.  gx1 = x1^2
    let mut gx1 = x1.clone();
    gx1.square();

    // 10. gx1 = gx1 + A
    gx1.add_assign(&curve.a);

    // 11. gx1 = gx1 * x1
    gx1.mul_assign(&x1);

    // 12. gx1 = gx1 + B             # gx1 = g(x1) = x1^3 + A * x1 + B
    gx1.add_assign(&curve.b);

    // 13.  x2 = tv1 * x1            # x2 = Z * u^2 * x1
    let mut x2 = tv1.clone();
    x2.mul_assign(&x1);

    // 14. tv2 = tv1 * tv2
    tv2.mul_assign(&tv1);

    // 15. gx2 = gx1 * tv2           # gx2 = (Z * u^2)^3 * gx1
    let mut gx2 = gx1.clone();
    gx2.mul_assign(&tv2);

    // 16.  e2 = is_square(gx1)
    let e2_legendre = legendre_symbol_fp2(&gx1);
    let e2 = e2_legendre == LegendreSymbol::Zero || e2_legendre == LegendreSymbol::QuadraticResidue;

    // 17.   x = CMOV(x2, x1, e2)    # If is_square(gx1), x = x1, else x = x2
    let x = if e2 {
        x1
    } else {
        x2
    };

    // 18.  y2 = CMOV(gx2, gx1, e2)  # If is_square(gx1), y2 = gx1, else y2 = gx2
    let y2 = if e2 {
        gx1
    } else {
        gx2
    };
    // 19.   y = sqrt(y2)
    let mut y = sqrt_ext2(&y2).expect("y2 is a square");

    // 20.  e3 = sgn0(u) == sgn0(y)  # Fix sign of y
    let u_sign = sign_of_fp2_be(&u);
    let y_sign = sign_of_fp2_be(&y);

    // 21.   y = CMOV(-y, y, e3)

    // our SignZero is equal to SignPlus
    match (u_sign, y_sign) {
        (Sign::Zero, Sign::Zero) |
        (Sign::Zero, Sign::SignPlus) | 
        (Sign::SignPlus, Sign::Zero) => {

        },
        (Sign::SignPlus, Sign::SignPlus) |
        (Sign::SignMinus, Sign::SignMinus) => {

        },
        (Sign::SignPlus, Sign::SignMinus) |
        (Sign::SignMinus, Sign::SignPlus) => {
            y.negate();
        },
        (Sign::Zero, Sign::SignMinus) | 
        (Sign::SignMinus, Sign::Zero) => {
            y.negate();
        },
    }

    // 22. return (x, y)
    (x, y)
}
