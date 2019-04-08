use super::field::{PrimeFieldElement, SizedPrimeField};
use super::representation::ElementRepr;
use crate::traits::{FieldElement, BitIterator};

#[derive(Eq, PartialEq, Clone, Copy)]
pub enum CurveType {
    Generic,
    AIsMinus3,
    AIsZero,
    BIsZero,
}

pub struct WeierstrassCurve<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> {
    pub(crate) field: &'a F,
    pub(crate) group: &'a G,
    a: PrimeFieldElement<'a, FE, F>,
    b: PrimeFieldElement<'a, FE, F>,
    curve_type: CurveType
}


impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> WeierstrassCurve<'a, FE, F, GE, G> {
    pub fn new(
        group: &'a G,
        a: PrimeFieldElement<'a, FE, F>, 
        b: PrimeFieldElement<'a, FE, F>,
        curve_type: CurveType
    ) -> Self {
        Self {
            field: &a.field,
            group: group,
            a: a,
            b: b,
            curve_type: curve_type
        }
    }
}

pub struct CurvePoint<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> {
    curve: &'a WeierstrassCurve<'a, FE, F, GE, G>,
    x: PrimeFieldElement<'a, FE, F>,
    y: PrimeFieldElement<'a, FE, F>,
    z: PrimeFieldElement<'a, FE, F>,
}

pub trait Group: Sized {
    fn add_assign(&mut self, other: &Self);
    fn add_assign_mixed(&mut self, other: &Self);
    fn sub_assign(&mut self, other: &Self);
    fn negate(&mut self);
    fn double(&mut self);
    fn mul<S: AsRef<[u64]>>(&self, exp: S) -> Self;
    fn is_zero(&self) -> bool;
}

impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> Clone for CurvePoint<'a, FE, F, GE, G> {
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

impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> CurvePoint<'a, FE, F, GE, G> {    
    pub fn zero(curve: &'a WeierstrassCurve<'a, FE, F, GE, G>) -> Self {
        Self {
            curve: curve,
            x: PrimeFieldElement::<'a, FE, F>::zero(&curve.field),
            y: PrimeFieldElement::<'a, FE, F>::one(&curve.field),
            z: PrimeFieldElement::<'a, FE, F>::zero(&curve.field),
        }
    }

    pub fn point_from_xy(
        curve: &'a WeierstrassCurve<'a, FE, F, GE, G>,
        x: PrimeFieldElement<'a, FE, F>, 
        y: PrimeFieldElement<'a, FE, F>
    ) -> CurvePoint<'a, FE, F, GE, G> {
        CurvePoint {
            curve: curve,
            x: x,
            y: y,
            z: PrimeFieldElement::<'a, FE, F>::one(&curve.field)
        }
    }

    pub fn normalize(&mut self) {
        if self.is_zero() {
            return;
        }
        let one = PrimeFieldElement::one(self.curve.field);
        if self.z == one {
            return;
        }

        let z_inv = self.z.mont_inverse().unwrap();
        let mut zinv_powered = z_inv.clone();
        zinv_powered.square();

        // X/Z^2
        self.x.mul_assign(&zinv_powered);

        // Y/Z^3
        zinv_powered.mul_assign(&z_inv);
        self.y.mul_assign(&zinv_powered);
    }

    pub fn into_xy(&self) -> (PrimeFieldElement<'a, FE, F>, PrimeFieldElement<'a, FE, F>) {
        if self.is_zero() {
            return (PrimeFieldElement::zero(self.curve.field), PrimeFieldElement::zero(self.curve.field));
        }

        let mut point = self.clone();
        point.normalize();

        (point.x, point.y)
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

        let one = PrimeFieldElement::<'a, FE, F>::one(&self.curve.field);
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
                // The two points are equal, so we double.
                self.x =  PrimeFieldElement::<'a, FE, F>::zero(&self.curve.field);
                self.y = PrimeFieldElement::<'a, FE, F>::one(&self.curve.field);
                self.z = PrimeFieldElement::<'a, FE, F>::zero(&self.curve.field);
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
            self.x = other.x.clone();
            self.y = other.y.clone();
            self.z = other.z.clone();
            return;
        }

        let one = PrimeFieldElement::one(self.curve.field);
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
        }
    }

    fn negate_impl(&mut self) {
        if !self.is_zero() {
            self.y.negate()
        }
    }

    fn mul_impl<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let one = PrimeFieldElement::<'a, FE, F>::one(&self.curve.field);
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
}

impl<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>, GE: ElementRepr, G: SizedPrimeField<Repr = GE>> Group for CurvePoint<'a, FE, F, GE, G> {
    fn add_assign(&mut self, other: &Self) {
        match self.curve.curve_type {
            CurveType::Generic => {
                self.add_assign_generic_impl(&other);
            },
            _ => {unimplemented!()}
        }
    }

    fn add_assign_mixed(&mut self, other: &Self) {
        match self.curve.curve_type {
            CurveType::Generic => {
                self.add_assign_mixed_generic_impl(&other);
            },
            _ => {unimplemented!()}
        }
    }

    fn sub_assign(&mut self, other: &Self) {
        let mut other_neg = other.clone();
        other_neg.negate();
        self.add_assign(&other_neg);
    }

    fn negate(&mut self) {
        match self.curve.curve_type {
            CurveType::Generic => {
                self.negate_impl();
            },
            _ => {unimplemented!()}
        }
    }

    fn mul<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        match self.curve.curve_type {
            CurveType::Generic => {
                return self.mul_impl(exp);
            },
            _ => {unimplemented!()}
        }
    }

    fn is_zero(&self) -> bool {
        match self.curve.curve_type {
            CurveType::Generic => {
                return self.is_zero_generic_impl();
            },
            _ => {unimplemented!()}
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
}
