#ifndef H_MNTengine
#define H_MNTengine

#include "../common.h"
#include "../curve.h"
#include "../fp.h"
#include "../extension_towers/fp2.h"
#include "../extension_towers/fp4.h"

template <class F, usize N>
struct AteDoubleCoefficients
{
    F c_h, c_4c, c_j, c_l;
};

template <class F, usize N>
struct AteAdditionCoefficients
{
    F c_l1, c_rz;
};

template <class F, usize N>
struct ExtendedCoordinates
{
    F x, y, z, t;
};

template <class F, usize N>
struct PrecomputedG1
{
    Fp<N> x, y;
    F x_by_twist, y_by_twist;
};

template <class F, usize N>
struct PrecomputedG2
{
    F x, y, x_over_twist, y_over_twist;
    std::vector<AteDoubleCoefficients<F, N>> double_coefficients;
    std::vector<AteAdditionCoefficients<F, N>> addition_coefficients;
};

template <class F1, class F2, usize N>
class MNTengine
{
    std::vector<u64> x;
    bool x_is_negative;
    std::vector<u64> exp_w0, exp_w1;
    bool exp_w0_is_negative;
    WeierstrassCurve<F1> const &curve_twist;
    F1 twist;

public:
    MNTengine(std::vector<u64> x, bool x_is_negative, std::vector<u64> exp_w0, std::vector<u64> exp_w1, bool exp_w0_is_negative,
              WeierstrassCurve<F1> const &curve_twist, F1 twist) : x(x), x_is_negative(x_is_negative), exp_w0(exp_w0), exp_w1(exp_w1), exp_w0_is_negative(exp_w0_is_negative), curve_twist(curve_twist), twist(twist)
    {
    }

    template <class C>
    std::optional<F2>
    pair(std::vector<std::tuple<CurvePoint<Fp<N>>, CurvePoint<F1>>> const &points, C const &context) const
    {
        if (points.size() == 0)
        {
            return {};
        }
        auto res = miller_loop(points, context);
        return final_exponentiation(res);
    }

private:
    template <class C>
    F2 miller_loop(std::vector<std::tuple<CurvePoint<Fp<N>>, CurvePoint<F1>>> const &points, C const &context) const
    {
        auto f = F2::one(context);
        for (auto it = points.cbegin(); it != points.cend(); it++)
        {
            f.mul(ate_pairing_loop(std::get<0>(*it), std::get<1>(*it), context));
        }
        return f;
    }

    template <class C>
    F2 ate_pairing_loop(
        CurvePoint<Fp<N>> const &point,
        CurvePoint<F1> const &twist_point, C const &context) const
    {
        assert(point.is_normalized());
        assert(twist_point.is_normalized());
        auto const twist_inv = this->twist.inverse().value();

        auto const p = precompute_g1(point);
        auto const q = precompute_g2(twist_point, twist_inv, context);
        auto l1_coeff = F1::zero(context);
        l1_coeff.c0 = p.x;
        l1_coeff.sub(q.x_over_twist);

        auto f = F2::one(context);

        usize dbl_idx = 0;
        usize add_idx = 0;
        // The for loop is executed for all bits (EXCEPT the MSB itself) of
        auto it = RevBitIterator(this->x);
        it.before(); // skip 1
        while (it.before())
        {
            auto const bit = *it;

            auto const dc = q.double_coefficients[dbl_idx];
            dbl_idx += 1;

            auto g_rr_at_p = F2::zero(context);

            auto t0 = dc.c_j;
            t0.mul(p.x_by_twist);
            t0.negate();
            t0.add(dc.c_l);
            t0.sub(dc.c_4c);

            auto t1 = dc.c_h;
            t1.mul(p.y_by_twist);

            g_rr_at_p.c0 = t0;
            g_rr_at_p.c1 = t1;

            f.square();
            f.mul(g_rr_at_p);

            if (bit)
            {
                auto const ac = q.addition_coefficients[add_idx];
                add_idx += 1;

                auto g_rq_at_p = F2::zero(context);

                auto t0 = ac.c_rz;
                t0.mul(p.y_by_twist);

                auto t = l1_coeff;
                t.mul(ac.c_l1);

                auto t1 = q.y_over_twist;
                t1.mul(ac.c_rz);
                t1.add(t);
                t1.negate();

                g_rq_at_p.c0 = t0;
                g_rq_at_p.c1 = t1;

                f.mul(g_rq_at_p);
            }
        }

        if (this->x_is_negative)
        {
            auto const ac = q.addition_coefficients[add_idx];

            auto g_rnegr_at_p = F2::zero(context);

            auto t0 = ac.c_rz;
            t0.mul(p.y_by_twist);

            auto t = l1_coeff;
            t.mul(ac.c_l1);

            auto t1 = q.y_over_twist;
            t1.mul(ac.c_rz);
            t1.add(t);
            t1.negate();

            g_rnegr_at_p.c0 = t0;
            g_rnegr_at_p.c1 = t1;

            f.mul(g_rnegr_at_p);
            f = f.inverse().value();
        }

        return f;
    }

    PrecomputedG1<F1, N> precompute_g1(CurvePoint<Fp<N>> const &g1_point) const
    {
        // not asserting normalization, it will be asserted in the loop
        auto x_twist = this->twist;
        x_twist.mul_by_fp(g1_point.x);

        auto y_twist = this->twist;
        y_twist.mul_by_fp(g1_point.y);

        return PrecomputedG1<F1, N>{
            g1_point.x,
            g1_point.y,
            x_twist,
            y_twist,
        };
    }

    template <class C>
    PrecomputedG2<F1, N> precompute_g2(CurvePoint<F1> const &g2_point, F1 const &twist_inv, C const &context) const
    {
        // not asserting normalization, it will be asserted in the loop
        // precompute addition and doubling coefficients
        auto x_over_twist = g2_point.x;
        x_over_twist.mul(twist_inv);

        auto y_over_twist = g2_point.y;
        y_over_twist.mul(twist_inv);

        auto g2_p = PrecomputedG2<F1, N>{
            g2_point.x,
            g2_point.y,
            x_over_twist,
            y_over_twist,
            std::vector<AteDoubleCoefficients<F1, N>>(),
            std::vector<AteAdditionCoefficients<F1, N>>(),
        };

        auto r = ExtendedCoordinates<F1, N>{
            g2_point.x,
            g2_point.y,
            F1::one(context),
            F1::one(context),
        };

        auto it = RevBitIterator(this->x);
        it.before(); // skip 1
        while (it.before())
        {
            auto const coeff = this->doubling_step(r);
            g2_p.double_coefficients.push_back(coeff);

            if (*it)
            {
                auto const coeff = this->addition_step(g2_point.x, g2_point.y, r);
                g2_p.addition_coefficients.push_back(coeff);
            }
        }

        if (this->x_is_negative)
        {
            auto const rz_inv = r.z.inverse().value();
            auto rz2_inv = rz_inv;
            rz2_inv.square();
            auto rz3_inv = rz_inv;
            rz3_inv.mul(rz2_inv);

            auto minus_r_affine_x = rz2_inv;
            minus_r_affine_x.mul(r.x);
            auto minus_r_affine_y = rz3_inv;
            minus_r_affine_y.mul(r.y);
            minus_r_affine_y.negate();

            auto const coeff = this->addition_step(
                minus_r_affine_x,
                minus_r_affine_y,
                r);

            g2_p.addition_coefficients.push_back(coeff);
        }

        return g2_p;
    }

    AteDoubleCoefficients<F1, N> doubling_step(ExtendedCoordinates<F1, N> &r) const
    {
        auto a = r.t;
        a.square();
        auto b = r.x;
        b.square();
        auto c = r.y;
        c.square();
        auto d = c;
        d.square();

        auto e = r.x;
        e.add(c);
        e.square();
        e.sub(b);
        e.sub(d);

        auto f = this->curve_twist.get_a();
        f.mul(a);
        f.add(b);
        f.add(b);
        f.add(b);

        auto g = f;
        g.square();

        auto d_eight = d;
        d_eight.mul2();
        d_eight.mul2();
        d_eight.mul2();

        auto t0 = e;
        t0.mul2();
        t0.mul2();

        auto x = g;
        x.sub(t0);

        auto y = e;
        y.mul2();
        y.sub(x);
        y.mul(f);
        y.sub(d_eight);

        auto h0 = r.z;
        h0.square();

        auto z = r.y;
        z.add(r.z);
        z.square();
        z.sub(c);
        z.sub(h0);

        auto t = z;
        t.square();

        auto c_h = z;
        c_h.add(r.t);
        c_h.square();
        c_h.sub(t);
        c_h.sub(a);

        auto c_4c = c;
        c_4c.mul2();
        c_4c.mul2();

        auto c_j = f;
        c_j.add(r.t);
        c_j.square();
        c_j.sub(g);
        c_j.sub(a);

        auto c_l = f;
        c_l.add(r.x);
        c_l.square();
        c_l.sub(g);
        c_l.sub(b);

        auto const coeff = AteDoubleCoefficients<F1, N>{c_h, c_4c, c_j, c_l};

        r.x = x;
        r.y = y;
        r.z = z;
        r.t = t;

        return coeff;
    }

    AteAdditionCoefficients<F1, N> addition_step(F1 const &x, F1 const &y, ExtendedCoordinates<F1, N> &r) const
    {
        auto a = y;
        a.square();
        auto b = r.t;
        b.mul(x);

        auto d = r.z;
        d.add(y);
        d.square();
        d.sub(a);
        d.sub(r.t);
        d.mul(r.t);

        auto h = b;
        h.sub(r.x);

        auto i = h;
        i.square();

        auto e = i;
        e.mul2();
        e.mul2();

        auto j = h;
        j.mul(e);

        auto v = r.x;
        v.mul(e);

        auto l1 = d;
        l1.sub(r.y);
        l1.sub(r.y);

        auto x0 = l1;
        x0.square();
        x0.sub(j);
        x0.sub(v);
        x0.sub(v);

        auto t0 = r.y;
        t0.mul2();
        t0.mul(j);

        auto y0 = v;
        y0.sub(x0);
        y0.mul(l1);
        y0.sub(t0);

        auto z = r.z;
        z.add(h);
        z.square();
        z.sub(r.t);
        z.sub(i);

        auto t = z;
        t.square();

        auto const coeff = AteAdditionCoefficients<F1, N>{l1, z};

        r.x = x0;
        r.y = y0;
        r.z = z;
        r.t = t;

        return coeff;
    }

    std::optional<F2> final_exponentiation(F2 const &f) const
    {
        auto const ovalue_inv = f.inverse();
        if (ovalue_inv)
        {
            return {};
        }
        auto const value_inv = ovalue_inv.value();
        auto const value_to_first_chunk = this->final_exponentiation_part_one(f, value_inv);
        auto const value_inv_to_first_chunk = this->final_exponentiation_part_one(value_inv, f);

        return this->final_exponentiation_part_two(value_to_first_chunk, value_inv_to_first_chunk);
    }

    virtual F2 final_exponentiation_part_one(F2 const &elt, F2 const &elt_inv) const = 0;

    F2 final_exponentiation_part_two(F2 const &elt, F2 const &elt_inv) const
    {
        auto elt_q = elt;
        elt_q.frobenius_map(1);

        auto w1_part = elt_q.cyclotomic_exp(this->exp_w1);
        if (this->exp_w0_is_negative)
        {
            auto w0_part = elt_inv.cyclotomic_exp(this->exp_w0);
            w1_part.mul(w0_part);
        }
        else
        {
            auto w0_part = elt.cyclotomic_exp(this->exp_w0);
            w1_part.mul(w0_part);
        }

        return w1_part;
    }
};

#endif