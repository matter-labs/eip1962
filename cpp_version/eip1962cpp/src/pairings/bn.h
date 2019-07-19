#ifndef H_BN
#define H_BN

#include "b_engine.h"

template <usize N>
class BNengine : public Bengine<N>
{
    std::vector<u64> six_u_plus_2;
    Fp2<N> non_residue_in_p_minus_one_over_2;

public:
    BNengine(std::vector<u64> u,
             bool u_is_negative,
             TwistType twist_type,
             WeierstrassCurve<Fp2<N>> const &curve_twist,
             Fp2<N> const &non_residue) : Bengine<N>(u, u_is_negative, twist_type, curve_twist),
                                          non_residue_in_p_minus_one_over_2(non_residue)
    {
        // Calculate six_u_plus_two
        six_u_plus_2 = this->u;
        mul_scalar(six_u_plus_2, 6);
        add_scalar(six_u_plus_2, 2);
        if (calculate_hamming_weight(six_u_plus_2) > MAX_BN_SIX_U_PLUS_TWO_HAMMING)
        {
            input_err("6*U + 2 has too large hamming weight");
        }

        // Calculate non_residue_in_p_minus_one_over_2
        constexpr Repr<N> one = {1};
        auto const p_minus_one_over_2 = cbn::shift_right(non_residue.field.mod() - one, 1);
        non_residue_in_p_minus_one_over_2 = non_residue.pow(p_minus_one_over_2);
    }

protected:
    Fp12<N> miller_loop(std::vector<std::tuple<CurvePoint<Fp<N>>, CurvePoint<Fp2<N>>>> const &points, FieldExtension2over3over2<N> const &context) const
    {

        std::vector<CurvePoint<Fp<N>>> g1_references;
        std::vector<std::vector<ThreePoint<N>>> prepared_coeffs;

        for (auto it = points.cbegin(); it != points.cend(); it++)
        {
            auto const p = std::get<0>(*it);
            auto const q = std::get<1>(*it);
            if (!p.is_zero() && !q.is_zero())
            {
                auto coeffs = prepare(q, context);
                prepared_coeffs.push_back(coeffs);
                g1_references.push_back(p);
            }
        }

        auto const n = prepared_coeffs.size();

        std::vector<usize> pc_indexes;
        pc_indexes.resize(n);

        auto f = Fp12<N>::one(context);
        auto it = RevBitIterator(six_u_plus_2);
        it.before(); //skip 1
        for (; it.before();)
        {
            auto const i = *it;
            f.square();

            this->for_ell(f, n, g1_references, prepared_coeffs, pc_indexes);

            if (i)
            {
                this->for_ell(f, n, g1_references, prepared_coeffs, pc_indexes);
            }
        }

        if (this->u_is_negative)
        {
            f.conjugate();
        }

        this->for_ell(f, n, g1_references, prepared_coeffs, pc_indexes);
        this->for_ell(f, n, g1_references, prepared_coeffs, pc_indexes);

        for (usize j = 0; j < n; j++)
        {
            assert(pc_indexes[j] == prepared_coeffs[j].size());
        }

        return f;
    }

    std::vector<ThreePoint<N>> prepare(CurvePoint<Fp2<N>> const &twist_point, FieldExtension2over3over2<N> const &context) const
    {
        assert(twist_point.is_normalized());

        auto two_inv_ = Fp<N>::one(context);
        two_inv_.mul2();
        auto const two_inv = two_inv_.inverse().value();

        std::vector<ThreePoint<N>> ell_coeffs;

        if (twist_point.is_zero())
        {
            return ell_coeffs;
        }

        auto r = CurvePoint<Fp2<N>>(twist_point.x, twist_point.y);

        auto it = RevBitIterator(six_u_plus_2);
        it.before(); //skip 1
        for (; it.before();)
        {
            ell_coeffs.push_back(this->doubling_step(r, two_inv));

            if (*it)
            {
                ell_coeffs.push_back(this->addition_step(r, twist_point));
            }
        }

        if (this->u_is_negative)
        {
            r.negate();
        }

        auto q = twist_point;

        q.x.c1.negate();
        FieldExtension3over2<N> const &field_3_2 = context;
        q.x.mul(field_3_2.frobenius_coeffs_c1[1]);

        q.y.c1.negate();
        q.y.mul(non_residue_in_p_minus_one_over_2);

        ell_coeffs.push_back(this->addition_step(r, q));

        auto minusq2 = twist_point;
        minusq2.x.mul(field_3_2.frobenius_coeffs_c1[2]);

        ell_coeffs.push_back(this->addition_step(r, minusq2));

        return ell_coeffs;
    }

    std::optional<Fp12<N>> final_exponentiation(Fp12<N> const &f) const
    {
        // use Zexe and pairing crate fused
        // https://eprint.iacr.org/2012/232.pdf

        // f1 = r.conjugate() = f^(p^6)
        auto f1 = f;
        f1.frobenius_map(6);

        if (auto of2 = f.inverse())
        {
            auto f2 = of2.value();
            auto r = f1;
            r.mul(f2);

            f2 = r;

            r.frobenius_map(2);
            r.mul(f2);

            auto fp = r;
            fp.frobenius_map(1);

            auto fp2 = r;
            fp2.frobenius_map(2);
            auto fp3 = fp2;
            fp3.frobenius_map(1);

            auto fu = r;
            this->exp_by_x(fu);
            // exp_by_x(fu, x);

            auto fu2 = fu;
            this->exp_by_x(fu2);
            // exp_by_x(fu2, x);

            auto fu3 = fu2;
            this->exp_by_x(fu3);
            // exp_by_x(fu3, x);

            auto y3 = fu;
            y3.frobenius_map(1);

            auto fu2p = fu2;
            fu2p.frobenius_map(1);

            auto fu3p = fu3;
            fu3p.frobenius_map(1);

            auto y2 = fu2;
            y2.frobenius_map(2);

            auto y0 = fp;
            y0.mul(fp2);
            y0.mul(fp3);

            auto y1 = r;
            y1.conjugate();

            auto y5 = fu2;
            y5.conjugate();

            y3.conjugate();

            auto y4 = fu;
            y4.mul(fu2p);
            y4.conjugate();

            auto y6 = fu3;
            y6.mul(fu3p);
            y6.conjugate();

            y6.square();
            y6.mul(y4);
            y6.mul(y5);

            auto t1 = y3;
            t1.mul(y5);
            t1.mul(y6);

            y6.mul(y2);

            t1.square();
            t1.mul(y6);
            t1.square();

            auto t0 = t1;
            t0.mul(y1);

            t1.mul(y0);

            t0.square();
            t0.mul(t1);

            return t0;
        }
        else
        {
            return {};
        }
    }
};

#endif
