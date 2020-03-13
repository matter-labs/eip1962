use crate::traits::*;

use crate::weierstrass::*;

#[derive(Clone)]
pub struct IsogenyParameters<C: CurveParameters> {
    pub map_degree: usize,
    pub k1: Vec<C::BaseFieldElement>,
    pub k2: Vec<C::BaseFieldElement>,
    pub k3: Vec<C::BaseFieldElement>,
    pub k4: Vec<C::BaseFieldElement>,
}

pub(crate) fn apply_isogeny_map<
    'a, 
    C: CurveParameters + 'a
> (
    x: &C::BaseFieldElement,
    y: &C::BaseFieldElement,
    params: &IsogenyParameters<C>,
    curve_params: &C
) -> (C::BaseFieldElement, C::BaseFieldElement) {
    debug_assert_eq!(params.map_degree + 1, params.k1.len());
    debug_assert_eq!(params.map_degree + 1, params.k2.len());
    debug_assert_eq!(params.map_degree + 1, params.k3.len());
    debug_assert_eq!(params.map_degree + 1, params.k4.len());

    let mut x_num = params.k1[params.map_degree].clone();
    let mut x_den = params.k2[params.map_degree].clone();
    let mut y_num = params.k3[params.map_degree].clone();
    let mut y_den = params.k4[params.map_degree].clone();

    // horner rule
    for i in (0..params.map_degree).rev() {
        x_num.mul_assign(&x);
        x_den.mul_assign(&x);
        y_num.mul_assign(&x);
        y_den.mul_assign(&x);

        x_num.add_assign(&params.k1[i]);
        x_den.add_assign(&params.k2[i]);
        y_num.add_assign(&params.k3[i]);
        y_den.add_assign(&params.k4[i]);
    }

    let x_den = x_den.inverse().unwrap_or(C::BaseFieldElement::zero(curve_params.params()));
    let y_den = y_den.inverse().unwrap_or(C::BaseFieldElement::zero(curve_params.params()));

    x_num.mul_assign(&x_den);
    y_num.mul_assign(&y_den);

    y_num.mul_assign(&y);

    (x_num, y_num)
}

