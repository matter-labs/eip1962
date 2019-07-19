#ifndef H_B_ENGINE
#define H_B_ENGINE

#include "../common.h"
#include "../curve.h"
#include "../fp.h"
#include "../extension_towers/fp2.h"
#include "../extension_towers/fp6_3.h"
#include "../extension_towers/fp12.h"

template <usize N>
using ThreePoint = std::tuple<Fp2<N>, Fp2<N>, Fp2<N>>;

template <usize N>
class Bengine
{
protected:
    std::vector<u64> u;
    bool u_is_negative;
    TwistType twist_type;
    WeierstrassCurve<Fp2<N>> const &curve_twist;

public:
    Bengine(std::vector<u64> u,
            bool u_is_negative,
            TwistType twist_type,
            WeierstrassCurve<Fp2<N>> const &curve_twist) : u(u), u_is_negative(u_is_negative), twist_type(twist_type), curve_twist(curve_twist) {}

    std::optional<Fp12<N>>
    pair(std::vector<std::tuple<CurvePoint<Fp<N>>, CurvePoint<Fp2<N>>>> const &points, FieldExtension2over3over2<N> const &context) const
    {
        if (points.size() == 0)
        {
            return {};
        }
        auto res = miller_loop(points, context);
        return final_exponentiation(res);
    }

protected:
    virtual Fp12<N> miller_loop(std::vector<std::tuple<CurvePoint<Fp<N>>, CurvePoint<Fp2<N>>>> const &points, FieldExtension2over3over2<N> const &context) const = 0;

    virtual std::vector<ThreePoint<N>> prepare(CurvePoint<Fp2<N>> const &twist_point, FieldExtension2over3over2<N> const &context) const = 0;

    virtual std::optional<Fp12<N>> final_exponentiation(Fp12<N> const &f) const = 0;

    ThreePoint<N> doubling_step(
        CurvePoint<Fp2<N>> &r,
        Fp<N> const &two_inv) const
    {
        // Use adapted formulas from ZEXE instead

        // X*Y/2
        auto a = r.x;
        a.mul(r.y);
        a.mul_by_fp(two_inv);

        // Y^2
        auto b = r.y;
        b.square();

        // Z^2
        auto c = r.z;
        c.square();

        auto e = curve_twist.get_b();

        // 3*Z^2
        auto t0 = c;
        t0.mul2();
        t0.add(c);

        // 3*b*Z^2
        e.mul(t0);

        // 9*b*Z^2
        auto f = e;
        f.mul2();
        f.add(e);

        // (Y^2 + 9*b*Z^2)/2
        auto g = b;
        g.add(f);
        g.mul_by_fp(two_inv);

        // (Y + Z)^2
        auto h = r.y;
        h.add(r.z);
        h.square();

        // (Y^2 + Z^2)
        auto t1 = b;
        t1.add(c);

        // 2*Y*Z
        h.sub(t1);

        // 3*b*Z^2 - Y^2
        auto i = e;
        i.sub(b);

        // X^2
        auto j = r.x;
        j.square();

        // (3*b*Z^2)^2
        auto e_square = e;
        e_square.square();

        // X = (Y^2 - 9*b*Z^2)*X*Y/2
        r.x = b;
        r.x.sub(f);
        r.x.mul(a);

        // 27*b^2*Z^4
        auto e_square_by_3 = e_square;
        e_square_by_3.mul2();
        e_square_by_3.add(e_square);

        // Y = ((Y^2 + 9*b*Z^2)/2)^2 - 27*b^2*Z^4
        r.y = g;
        r.y.square();
        r.y.sub(e_square_by_3);

        // Z = 2*Y^3*Z
        r.z = b;
        r.z.mul(h);

        // 3*X^2
        auto j_by_three = j;
        j_by_three.mul2();
        j_by_three.add(j);

        // - 2*Y*Z
        h.negate();

        switch (twist_type)
        {
        case M:
            return std::tuple(i, j_by_three, h);
        case D: // (0, 3, 4) = (-2*Y*Z, 3*X^2, 3*b*Z^2 - Y^2)
            return std::tuple(h, j_by_three, i);
        }
        unreachable("");
    }

    ThreePoint<N> addition_step(
        CurvePoint<Fp2<N>> &r,
        CurvePoint<Fp2<N>> const &q) const
    {
        assert(q.is_normalized());
        // use adapted zexe formulas too instead of ones from pairing crate
        // Capitals are coors of R (homogenious), normals are coordinates of Q (affine)
        // Y - y*Z
        auto theta = q.y;
        theta.mul(r.z);
        theta.negate();
        theta.add(r.y);

        // X - x*Z
        auto lambda = q.x;
        lambda.mul(r.z);
        lambda.negate();
        lambda.add(r.x);

        // Theta^2
        auto c = theta;
        c.square();

        // Lambda^2
        auto d = lambda;
        d.square();

        // Lambda^3
        auto e = lambda;
        e.mul(d);

        // Theta^2 * Z
        auto f = r.z;
        f.mul(c);

        // Lambda^2 * X
        auto g = r.x;
        g.mul(d);

        // Lambda^3 + Theta^2 * Z - 2*Lambda^2 * X
        auto h = g;
        h.mul2();
        h.negate();
        h.add(e);
        h.add(f);

        r.x = lambda;
        r.x.mul(h);

        // (Lambda^2 * X - H)*Theta
        auto t0 = g;
        t0.sub(h);
        t0.mul(theta);

        // Y = (Lambda^2 * X - H)*Theta - Lambda^3 * Y
        r.y.mul(e);
        r.y.negate();
        r.y.add(t0);

        // Z = Lambda^3 * Z
        r.z.mul(e);

        // Lambda*y
        auto t1 = lambda;
        t1.mul(q.y);

        // Theta*x - Lambda*y
        auto j = theta;
        j.mul(q.x);
        j.sub(t1);

        theta.negate();

        // lambda.negate();
        switch (twist_type)
        {
        case M:
            return std::tuple(j, theta, lambda);
        case D: // (0, 3, 4) = (lambda, -theta, Theta*x - Lambda*y)
            return std::tuple(lambda, theta, j);
        }
        unreachable("");
    }

    void for_ell(Fp12<N> &f, usize n, std::vector<CurvePoint<Fp<N>>> const &g1_references, std::vector<std::vector<ThreePoint<N>>> const &prepared_coeffs, std::vector<usize> &pc_indexes) const
    {
        for (usize j = 0; j < n; j++)
        {
            auto const p = g1_references[j];
            auto const coeffs = prepared_coeffs[j][pc_indexes[j]];
            pc_indexes[j]++;
            ell(f, coeffs, p);
        }
    }

    void ell(
        Fp12<N> &f,
        ThreePoint<N> const &coeffs,
        CurvePoint<Fp<N>> const &p) const
    {
        assert(p.is_normalized());
        auto c0 = std::get<0>(coeffs);
        auto c1 = std::get<1>(coeffs);
        auto c2 = std::get<2>(coeffs);

        switch (twist_type)
        {
        case M:
        {
            c2.mul_by_fp(p.y);
            c1.mul_by_fp(p.x);
            f.mul_by_014(c0, c1, c2);
            break;
        }
        case D:
        {
            c0.mul_by_fp(p.y);
            c1.mul_by_fp(p.x);
            f.mul_by_034(c0, c1, c2);
            break;
        }
        }
    }

    void exp_by_x(Fp12<N> &f) const
    {
        f = f.cyclotomic_exp(this->u);
        if (u_is_negative)
        {
            f.conjugate();
        }
    }
};

#endif
