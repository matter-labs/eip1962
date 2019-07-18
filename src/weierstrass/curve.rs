use crate::traits::{FieldElement, BitIterator};
use super::{CurveType, Group};
use super::CurveParameters;
use crate::traits::ZeroAndOne;

pub struct WeierstrassCurve<'a, C: CurveParameters> {
    // pub(crate) base_field: CurveParameters::BaseFieldParameters,
    pub(crate) a: C::BaseFieldElement,
    pub(crate) b: C::BaseFieldElement,
    pub(crate) curve_type: CurveType,
    pub(crate) subgroup_order_repr: Vec<u64>,
    pub(crate) params: &'a C
}

// pub struct WeierstrassCurve<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, C: CurveParameters> {
//     pub(crate) base_field: &'a F,
//     pub(crate) a: Fp<'a, FE, F>,
//     pub(crate) b: Fp<'a, FE, F>,
//     pub(crate) curve_type: CurveType,
//     pub(crate) subgroup_order_repr: Vec<u64>,
//     pub(crate) params: 
// }


impl<'a, C: CurveParameters> WeierstrassCurve<'a, C> {
    pub(crate) fn new(
        subgroup_order: Vec<u64>,
        a: C::BaseFieldElement, 
        b: C::BaseFieldElement,
        params: &'a C
    ) -> Result<Self, ()> {
        let mut curve_type = CurveType::Generic;
        if a.is_zero() {
            if b.is_zero() {
                return Err(());
            }
            curve_type = CurveType::AIsZero;
        }
        if b.is_zero() {
            return Err(());
        }

        Ok(Self {
            a: a,
            b: b,
            curve_type: curve_type,
            subgroup_order_repr: subgroup_order,
            params: params
        })
    }
}

pub struct CurvePoint<'a, C: CurveParameters> {
    pub(crate) curve: &'a WeierstrassCurve<'a, C>,
    pub(crate) x: C::BaseFieldElement,
    pub(crate) y: C::BaseFieldElement,
    pub(crate) z: C::BaseFieldElement,
}

impl<'a, C: CurveParameters> Clone for CurvePoint<'a, C> {
    #[inline(always)]
    fn clone(&self) -> Self {
        Self {
            curve: &self.curve,
            x: self.x.clone(),
            y: self.y.clone(),
            z: self.z.clone()
        }
    }
}

impl<'a, C: CurveParameters> CurvePoint<'a, C> {    
    pub fn zero(curve: &'a WeierstrassCurve<C>) -> Self {
        Self {
            curve: curve,
            x: C::BaseFieldElement::zero(curve.params.params()),
            y: C::BaseFieldElement::one(curve.params.params()),
            z: C::BaseFieldElement::zero(curve.params.params()),
        }
    }

    pub fn is_on_curve(&self) -> bool {
        if self.is_zero() {
            return true;
        }
        
        debug_assert!(self.is_normalized());

        let mut rhs = self.y.clone();
        rhs.square();

        let mut lhs = self.curve.b.clone();
        let mut ax = self.x.clone();
        ax.mul_assign(&self.curve.a);
        lhs.add_assign(&ax);

        let mut x_3 = self.x.clone();
        x_3.square();
        x_3.mul_assign(&self.x);
        lhs.add_assign(&x_3);

        rhs == lhs
    }

    pub fn point_from_xy(
        curve: &'a WeierstrassCurve<'a, C>,
        x: C::BaseFieldElement, 
        y: C::BaseFieldElement
    ) -> CurvePoint<'a, C> {
        if x.is_zero() && y.is_zero() {
            return Self::zero(curve);
        }

        CurvePoint {
            curve: curve,
            x: x,
            y: y,
            z: C::BaseFieldElement::one(curve.params.params())
        }
    }

    pub fn is_normalized(&self) -> bool {
        if self.is_zero() {
            return true;
        }

        let one = C::BaseFieldElement::one(self.curve.params.params());
        
        self.z == one
    }

    pub fn normalize(&mut self) {
        if self.is_zero() {
            return;
        }
        let one = C::BaseFieldElement::one(self.curve.params.params());
        if self.z == one {
            return;
        }

        match self.z.inverse() {
            Some(z_inv) => {
                let mut zinv_powered = z_inv.clone();
                zinv_powered.square();

                // X/Z^2
                self.x.mul_assign(&zinv_powered);

                // Y/Z^3
                zinv_powered.mul_assign(&z_inv);
                self.y.mul_assign(&zinv_powered);

                self.z = one;
            },
            None => {
                self.z = C::BaseFieldElement::zero(self.curve.params.params());
            }
        }
    }

    pub fn into_xy(&self) -> (C::BaseFieldElement, C::BaseFieldElement) {
        if self.is_zero() {
            return (C::BaseFieldElement::zero(self.curve.params.params()),
                    C::BaseFieldElement::zero(self.curve.params.params()));
        }

        let mut point = self.clone();

        if point.is_normalized() {
            (point.x, point.y)
        } else {
            point.normalize();
            point.into_xy()
        }
    }

    pub fn into_xy_from_homogenious(&self) -> (C::BaseFieldElement, C::BaseFieldElement) {
        if self.is_zero() {
            return (C::BaseFieldElement::zero(self.curve.params.params()),
                    C::BaseFieldElement::zero(self.curve.params.params()) );
        }

        let z_inv = self.z.inverse().unwrap_or_else(|| C::BaseFieldElement::zero(self.curve.params.params()));

        let mut x = self.x.clone();
        x.mul_assign(&z_inv);

        let mut y = self.y.clone();
        y.mul_assign(&z_inv);

        (x, y)
    }
    
    fn add_assign_generic_impl(&mut self, other: &Self) {
        if self.is_zero() {
            self.x = other.x.clone();
            self.y = other.y.clone();
            self.z = other.z.clone();
            return;
        }

        if other.is_zero() {
            return;
        }

        let one = C::BaseFieldElement::one(self.curve.params.params());
        if other.z == one {
            self.add_assign_mixed_generic_impl(&other);
            return;
        }

        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl

        // Z1Z1 = Z1^2
        let mut z1z1 = self.z.clone();
        z1z1.square();

        // Z2Z2 = Z2^2
        let mut z2z2 = other.z.clone();
        z2z2.square();

        // U1 = X1*Z2Z2
        let mut u1 = self.x.clone();
        u1.mul_assign(&z2z2);

        // U2 = X2*Z1Z1
        let mut u2 = other.x.clone();
        u2.mul_assign(&z1z1);

        // S1 = Y1*Z2*Z2Z2
        let mut s1 = self.y.clone();
        s1.mul_assign(&other.z);
        s1.mul_assign(&z2z2);

        // S2 = Y2*Z1*Z1Z1
        let mut s2 = other.y.clone();
        s2.mul_assign(&self.z);
        s2.mul_assign(&z1z1);

        if u1 == u2 && s1 == s2 {
            // The two points are equal, so we double.
            self.double();
        } else {
            // If we're adding -a and a together, self.z becomes zero as H becomes zero.

            if u1 == u2 {
                *self = Self::zero(self.curve);
                // self.x = C::BaseFieldElement::zero(self.curve.params.params());
                // self.y = C::BaseFieldElement::one(self.curve.params.params());
                // self.z = C::BaseFieldElement::zero(self.curve.params.params());
                return;
            }

            // H = U2-U1
            let mut h = u2.clone();
            h.sub_assign(&u1);

            // I = (2*H)^2
            let mut i = h.clone();
            i.double();
            i.square();

            // J = H*I
            let mut j = h.clone();
            j.mul_assign(&i);

            // r = 2*(S2-S1)
            let mut r = s2.clone();
            r.sub_assign(&s1);
            r.double();

            // V = U1*I
            let mut v = u1.clone();
            v.mul_assign(&i);

            // X3 = r^2 - J - 2*V
            self.x = r.clone();
            self.x.square();
            self.x.sub_assign(&j);
            self.x.sub_assign(&v);
            self.x.sub_assign(&v);

            // Y3 = r*(V - X3) - 2*S1*J
            self.y = v.clone();
            self.y.sub_assign(&self.x);
            self.y.mul_assign(&r);
            s1.mul_assign(&j); // S1 = S1 * J * 2
            s1.double();
            self.y.sub_assign(&s1);

            // Z3 = ((Z1+Z2)^2 - Z1Z1 - Z2Z2)*H
            self.z.add_assign(&other.z);
            self.z.square();
            self.z.sub_assign(&z1z1);
            self.z.sub_assign(&z2z2);
            self.z.mul_assign(&h);
        }
    }

    fn add_assign_mixed_generic_impl(&mut self, other: &Self) {
        if other.is_zero() {
            return;
        }

        if self.is_zero() {
            *self = other.clone();
            // self.x = other.x.clone();
            // self.y = other.y.clone();
            // self.z = other.z.clone();
            return;
        }

        let one = C::BaseFieldElement::one(self.curve.params.params());
        if other.z != one {
            self.add_assign_generic_impl(&other);
            return;
        }

        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl

        // Z1Z1 = Z1^2
        let mut z1z1 = self.z.clone();
        z1z1.square();

        // U2 = X2*Z1Z1
        let mut u2 = other.x.clone();
        u2.mul_assign(&z1z1);

        // S2 = Y2*Z1*Z1Z1
        let mut s2 = other.y.clone();
        s2.mul_assign(&self.z);
        s2.mul_assign(&z1z1);

        if self.x == u2 && self.y == s2 {
            // The two points are equal, so we double.
            self.double();
        } else {
            // If we're adding -a and a together, self.z becomes zero as H becomes zero.

            // H = U2-X1
            let mut h = u2.clone();
            h.sub_assign(&self.x);

            // HH = H^2
            let mut hh = h.clone();
            hh.square();

            // I = 4*HH
            let mut i = hh.clone();
            i.double();
            i.double();

            // J = H*I
            let mut j = h.clone();
            j.mul_assign(&i);

            // r = 2*(S2-Y1)
            let mut r = s2.clone();
            r.sub_assign(&self.y);
            r.double();

            // V = X1*I
            let mut v = self.x.clone();
            v.mul_assign(&i);

            // X3 = r^2 - J - 2*V
            self.x = r.clone();
            self.x.square();
            self.x.sub_assign(&j);
            self.x.sub_assign(&v);
            self.x.sub_assign(&v);

            // Y3 = r*(V-X3)-2*Y1*J
            j.mul_assign(&self.y); // J = 2*Y1*J
            j.double();
            self.y = v.clone();
            self.y.sub_assign(&self.x);
            self.y.mul_assign(&r);
            self.y.sub_assign(&j);

            // Z3 = (Z1+H)^2-Z1Z1-HH
            self.z.add_assign(&h);
            self.z.square();
            self.z.sub_assign(&z1z1);
            self.z.sub_assign(&hh);

            if self.z.is_zero() {
                *self = Self::zero(self.curve);
            }
        }
    }

    fn negate_impl(&mut self) {
        if !self.is_zero() {
            self.y.negate()
        }
    }

    fn mul_impl<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let one = C::BaseFieldElement::one(self.curve.params.params());
        if self.z == one {
            return self.mul_impl_mixed_addition(exp);
        }

        let mut res = Self::zero(&self.curve);

        let mut found_one = false;

        for i in BitIterator::new(exp)
        {
            if found_one {
                res.double();
            } else {
                found_one = i;
            }

            if i {
                res.add_assign(self);
            }
        }

        res
    }

    pub fn wnaf_mul_impl<S: crate::representation::IntoWnaf>(&self, exp: S) -> Self {
        // let one = Fp::<'a, FE, F>::one(&self.curve.field);
        // if self.z == one {
        //     return self.mul_impl_mixed_addition(exp);
        // }

        const WINDOW_SIZE: u32 = 3;

        let mut precomp_table = vec![Self::zero(&self.curve); (1 << (WINDOW_SIZE-1)) as usize];

        let index_for_positive = (1 << (WINDOW_SIZE-2)) as usize;

        let mut two_self = self.clone();
        two_self.double();

        let mut precomp = self.clone();
        precomp_table[index_for_positive] = precomp.clone();
        let mut neg_precomp = precomp.clone();
        neg_precomp.negate();
        precomp_table[index_for_positive-1] = neg_precomp;

        for i in 1..index_for_positive {
            precomp.add_assign(&two_self);
            precomp_table[index_for_positive+i] = precomp.clone();
            let mut neg_precomp = precomp.clone();
            neg_precomp.negate();
            precomp_table[index_for_positive-1-i] = neg_precomp;
        }

        let wnaf = exp.wnaf(WINDOW_SIZE);

        let mut res = Self::zero(&self.curve);
        let mut found_nonzero = false;

        for w in wnaf.into_iter().rev() {
            if found_nonzero {
                res.double();
            }
            if w != 0 {
                found_nonzero = true;
                if w > 0 {
                    let idx = (w >> 1) as usize;
                    res.add_assign(&precomp_table[index_for_positive + idx]);
                } else {
                    let idx = ((-w) >> 1) as usize;
                    res.add_assign(&precomp_table[index_for_positive - 1 - idx]);
                }
            }
        }
        
        res
    }

    fn mul_impl_mixed_addition<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Self::zero(&self.curve);

        let mut found_one = false;

        for i in BitIterator::new(exp)
        {
            if found_one {
                res.double();
            } else {
                found_one = i;
            }

            if i {
                res.add_assign_mixed(self);
            }
        }

        res
    }

    fn is_zero_generic_impl(&self) -> bool {
        return self.z.is_zero();
    }

    fn double_generic_impl(&mut self) {
        if self.is_zero() {
            return;
        }

        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian.html#doubling-dbl-2007-bl

        // A = X1^2
        let mut a = self.x.clone();
        a.square();

        // B = Y1^2
        let mut b = self.y.clone();
        b.square();

        // C = B^2 = Y1^4
        let mut c = b.clone();
        c.square();

        let mut z_2 = self.z.clone();
        z_2.square();

        // D = 2*((X1+B)2-A-C)
        let mut d = self.x.clone();
        d.add_assign(&b);
        d.square();
        d.sub_assign(&a);
        d.sub_assign(&c);
        d.double();

        // E = 3*A + curve_a*z^4
        let mut e = a.clone();
        e.double();
        e.add_assign(&a);

        // curve_a*z^4
        let mut a_z_4 = z_2.clone();
        a_z_4.square();
        a_z_4.mul_assign(&self.curve.a);

        e.add_assign(&a_z_4);

        // T = D^2
        let mut t = d.clone();
        t.double();

        // F = E^2 - 2*D
        let mut f = e.clone();
        f.square();
        f.sub_assign(&t);

        self.x = f;

        // Z3 = (Y1+Z1)^2-B-Z^2
        self.z.add_assign(&self.y);
        self.z.square();
        self.z.sub_assign(&b);
        self.z.sub_assign(&z_2);

        // Y3 = E*(D-X3)-8*C 
        self.y = d;
        self.y.sub_assign(&self.x);
        self.y.mul_assign(&e);
        c.double();
        c.double();
        c.double();
        self.y.sub_assign(&c);
    }

    fn double_a_is_zero_impl(&mut self) {
        if self.is_zero() {
            return;
        }

        // Other than the point at infinity, no points on E or E'
        // can double to equal the point at infinity, as y=0 is
        // never true for points on the curve.

        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l

        // A = X1^2
        let mut a = self.x.clone();
        a.square();

        // B = Y1^2
        let mut b = self.y.clone();
        b.square();

        // C = B^2
        let mut c = b.clone();
        c.square();

        // D = 2*((X1+B)2-A-C)
        let mut d = self.x.clone();
        d.add_assign(&b);
        d.square();
        d.sub_assign(&a);
        d.sub_assign(&c);
        d.double();

        // E = 3*A
        let mut e = a.clone();
        e.double();
        e.add_assign(&a);

        // F = E^2
        let mut f = e.clone();
        f.square();

        // Z3 = 2*Y1*Z1
        self.z.mul_assign(&self.y);
        self.z.double();

        // X3 = F-2*D
        self.x = f;
        self.x.sub_assign(&d);
        self.x.sub_assign(&d);

        // Y3 = E*(D-X3)-8*C
        self.y = d;
        self.y.sub_assign(&self.x);
        self.y.mul_assign(&e);
        c.double();
        c.double();
        c.double();
        self.y.sub_assign(&c);
    }

    fn check_correct_subgroup_impl(&self) -> bool {
        let p = self.mul_impl(&self.curve.subgroup_order_repr);

        p.is_zero_generic_impl()
    }
}

impl<'a, C: CurveParameters> Group for CurvePoint<'a, C> {
    fn add_assign(&mut self, other: &Self) {
        match self.curve.curve_type {
            _ => {
                self.add_assign_generic_impl(&other);
            },
        }
    }

    fn add_assign_mixed(&mut self, other: &Self) {
        match self.curve.curve_type {
            _ => {
                self.add_assign_mixed_generic_impl(&other);
            },
        }
    }

    fn sub_assign(&mut self, other: &Self) {
        let mut other_neg = other.clone();
        other_neg.negate();
        self.add_assign(&other_neg);
    }

    fn negate(&mut self) {
        match self.curve.curve_type {
            _ => {
                self.negate_impl();
            },
        }
    }

    fn mul<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        match self.curve.curve_type {
            _ => {
                return self.mul_impl(exp);
            },
        }
    }

    fn is_zero(&self) -> bool {
        match self.curve.curve_type {
            _ => {
                return self.is_zero_generic_impl();
            },
        }
    }

    fn double(&mut self) {
        match self.curve.curve_type {
            CurveType::Generic => {
                self.double_generic_impl();
            },
            CurveType::AIsZero => {
                self.double_a_is_zero_impl();
            }
            _ => {unimplemented!()}
        }
    }

    fn wnaf_mul<S: crate::representation::IntoWnaf>(&self, exp: S) -> Self {
        match self.curve.curve_type {
            _ => {
                return self.wnaf_mul_impl(exp);
            },
        }
    }

    fn check_correct_subgroup(&self) -> bool {
        match self.curve.curve_type {
            _ => {
                return self.check_correct_subgroup_impl();
            },
        }
    }
}
