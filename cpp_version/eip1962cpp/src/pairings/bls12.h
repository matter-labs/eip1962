#ifndef H_BLS12
#define H_BLS12

#include "b_engine.h"
#include "../constants.h"

template <usize N>
class BLS12engine : public Bengine<N>
{

public:
    BLS12engine(std::vector<u64> u,
                bool u_is_negative,
                TwistType twist_type,
                WeierstrassCurve<Fp2<N>> const &curve_twist, Fp2<N> const &non_residue) : Bengine<N>(u, u_is_negative, twist_type, curve_twist)
    {
        UNUSED(non_residue);
        if (calculate_hamming_weight(this->u) > MAX_BLS12_X_HAMMING)
        {
            input_err("X has too large hamming weight");
        }
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
        auto it = RevBitIterator(this->u);
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

        auto it = RevBitIterator(this->u);
        it.before(); //skip 1
        for (; it.before();)
        {
            ell_coeffs.push_back(this->doubling_step(r, two_inv));

            if (*it)
            {
                ell_coeffs.push_back(this->addition_step(r, twist_point));
            }
        }

        return ell_coeffs;
    }

    std::optional<Fp12<N>> final_exponentiation(Fp12<N> const &f) const
    {
        // Computing the final exponentation following
        // https://eprint.iacr.org/2016/130.pdf.
        // We don't use their "faster" formula because it is difficult to make
        // it work for curves with odd `P::X`.
        // Hence we implement the algorithm from Table 1 below.

        // f1 = r.conjugate() = f^(p^6)
        auto f1 = f;
        f1.frobenius_map(6);

        if (auto of2 = f.inverse())
        {
            auto f2 = of2.value();
            // f2 = f^(-1);
            // r = f^(p^6 - 1)
            auto r = f1;
            r.mul(f2);

            // f2 = f^(p^6 - 1)
            f2 = r;
            // r = f^((p^6 - 1)(p^2))
            r.frobenius_map(2);

            // r = f^((p^6 - 1)(p^2) + (p^6 - 1))
            // r = f^((p^6 - 1)(p^2 + 1))
            r.mul(f2);

            // Hard part of the final exponentation is below:
            // From https://eprint.iacr.org/2016/130.pdf, Table 1
            auto y0 = r;
            y0.cyclotomic_square();
            y0.conjugate();

            auto y5 = r;
            this->exp_by_x(y5);

            auto y1 = y5;
            y1.cyclotomic_square();

            auto y3 = y0;
            y3.mul(y5);

            auto e0 = y3;
            this->exp_by_x(e0);

            auto y2 = e0;
            this->exp_by_x(y2);

            auto y4 = y2;
            this->exp_by_x(y4);
            y4.mul(y1);

            auto e1 = y4;
            this->exp_by_x(e1);

            y3.conjugate();
            e1.mul(y3);
            e1.mul(r);

            auto e3 = r;
            e3.conjugate();
            e0.mul(r);
            e0.frobenius_map(3);

            y4.mul(e3);
            y4.frobenius_map(1);

            y5.mul(y2);
            y5.frobenius_map(2);

            y5.mul(e0);
            y5.mul(y4);
            y5.mul(e1);

            return y5;
        }
        else
        {
            return {};
        }
    }
};

#endif