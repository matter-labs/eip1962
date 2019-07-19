#ifndef H_FP12
#define H_FP12

#include "../common.h"
#include "../element.h"
#include "fp6_3.h"
#include "fpM2.h"
#include "../field.h"
#include "common.h"

template <usize N>
class FieldExtension2over3over2 : public FieldExtension3over2<N>
{
public:
    std::array<Fp2<N>, 12> frobenius_coeffs_c1;

    FieldExtension2over3over2(FieldExtension3over2<N> const &field, WindowExpBase<Fp2<N>> const &base) : FieldExtension3over2<N>(field), frobenius_coeffs_c1({Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field)})
    {

        // Fq2(u + 1)**(((q^0) - 1) / 6)
        auto const f_0 = Fp2<N>::one(field);

        // Fq2(u + 1)**(((q^1) - 1) / 6)

        // 1
        auto const q_power_1 = field.mod();
        auto const f_1 = base.exponentiate(calc_frobenius_factor_2(q_power_1, 6, "Fp12"));

        // 2
        auto const q_power_2 = q_power_1 * field.mod();
        auto const f_2 = base.exponentiate(calc_frobenius_factor_2(q_power_2, 6, "Fp12"));

        // 3
        auto const q_power_3 = q_power_2 * field.mod();
        auto const f_3 = base.exponentiate(calc_frobenius_factor_2(q_power_3, 6, "Fp12"));

        // 6
        auto const q_power_6 = q_power_3 * q_power_3;
        auto const f_6 = base.exponentiate(calc_frobenius_factor_2(q_power_6, 6, "Fp12"));

        auto const f_4 = Fp2<N>::zero(field);
        auto const f_5 = Fp2<N>::zero(field);

        auto const f_7 = Fp2<N>::zero(field);
        auto const f_8 = Fp2<N>::zero(field);
        auto const f_9 = Fp2<N>::zero(field);
        auto const f_10 = Fp2<N>::zero(field);
        auto const f_11 = Fp2<N>::zero(field);

        std::array<Fp2<N>, 12> calc_frobenius_coeffs_c1 = {f_0, f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9, f_10, f_11};
        frobenius_coeffs_c1 = calc_frobenius_coeffs_c1;
    }

    void mul_by_nonresidue(Fp6_3<N> &el) const
    {
        // IMPORTANT: This only works cause the structure of extension field for Fp12
        // is w^2 - v = 0!
        // take an element in Fp6 that is 3 over 2 and multiply by non-residue
        // (c0 + c1 * v + c2 * v^2)*v with v^3 - xi = 0 -> (c2*xi + c0 * v + c1 * v^2)
        auto new_c0 = el.c2;
        el.field.mul_by_nonresidue(new_c0);
        el.c2 = el.c1;
        el.c1 = el.c0;
        el.c0 = new_c0;
    }
};

template <usize N>
class Fp12 : public FpM2<Fp6_3<N>, FieldExtension2over3over2<N>, Fp12<N>, N>
{
public:
    Fp12(Fp6_3<N> c0, Fp6_3<N> c1, FieldExtension2over3over2<N> const &field) : FpM2<Fp6_3<N>, FieldExtension2over3over2<N>, Fp12<N>, N>(c0, c1, field) {}

    auto operator=(Fp12<N> const &other)
    {
        this->c0 = other.c0;
        this->c1 = other.c1;
    }

    void cyclotomic_square()
    {
        FieldExtension3over2<N> const &field_2 = this->field;

        auto const z0 = this->c0.c0;
        auto const z4 = this->c0.c1;
        auto const z3 = this->c0.c2;
        auto const z2 = this->c1.c0;
        auto const z1 = this->c1.c1;
        auto const z5 = this->c1.c2;

        // t0 + t1*y = (z0 + z1*y)^2 = a^2
        auto tmp = z0;
        tmp.mul(z1);

        auto a0_0 = z0;
        a0_0.add(z1);
        auto a1_0 = z1;
        field_2.mul_by_nonresidue(a1_0);
        a1_0.add(z0);

        auto a2_0 = tmp;
        field_2.mul_by_nonresidue(a2_0);

        auto t0 = a0_0;
        t0.mul(a1_0);
        t0.sub(tmp);
        t0.sub(a2_0);
        auto t1 = tmp;
        t1.mul2();

        // t2 + t3*y = (z2 + z3*y)^2 = b^2
        auto tmp2 = z2;
        tmp2.mul(z3);

        auto a0_1 = z2;
        a0_1.add(z3);
        auto a1_1 = z3;
        field_2.mul_by_nonresidue(a1_1);
        a1_1.add(z2);

        auto a2_1 = tmp2;
        field_2.mul_by_nonresidue(a2_1);

        auto t2 = a0_1;
        t2.mul(a1_1);
        t2.sub(tmp2);
        t2.sub(a2_1);

        auto t3 = tmp2;
        t3.mul2();

        // t4 + t5*y = (z4 + z5*y)^2 = c^2
        auto tmp3 = z4;
        tmp3.mul(z5);

        auto a0_2 = z4;
        a0_2.add(z5);
        auto a1_2 = z5;
        field_2.mul_by_nonresidue(a1_2);
        a1_2.add(z4);

        auto a2_2 = tmp3;
        field_2.mul_by_nonresidue(a2_2);

        auto t4 = a0_2;
        t4.mul(a1_2);
        t4.sub(tmp3);
        t4.sub(a2_2);

        auto t5 = tmp3;
        t5.mul2();

        // for A

        // g0 = 3 * t0 - 2 * z0
        auto g0 = t0;
        g0.sub(z0);
        g0.mul2();
        g0.add(t0);

        this->c0.c0 = g0;

        // g1 = 3 * t1 + 2 * z1
        auto g1 = t1;
        g1.add(z1);
        g1.mul2();
        g1.add(t1);
        this->c1.c1 = g1;

        // for B

        // g2 = 3 * (xi * t5) + 2 * z2
        auto tmp4 = t5;
        field_2.mul_by_nonresidue(tmp4);
        auto g2 = tmp4;
        g2.add(z2);
        g2.mul2();
        g2.add(tmp4);
        this->c1.c0 = g2;

        // g3 = 3 * t4 - 2 * z3
        auto g3 = t4;
        g3.sub(z3);
        g3.mul2();
        g3.add(t4);
        this->c0.c2 = g3;

        // for C

        // g4 = 3 * t2 - 2 * z4
        auto g4 = t2;
        g4.sub(z4);
        g4.mul2();
        g4.add(t2);
        this->c0.c1 = g4;

        // g5 = 3 * t3 + 2 * z5
        auto g5 = t3;
        g5.add(z5);
        g5.mul2();
        g5.add(t3);
        this->c1.c2 = g5;
    }

    Fp12<N> cyclotomic_exp(std::vector<u64> const &exp) const
    {
        auto res = one();

        auto found_one = false;

        for (auto it = RevBitIterator(exp); it.before();)
        {
            auto const i = *it;
            if (found_one)
            {
                res.cyclotomic_square();
            }
            else
            {
                found_one = i;
            }

            if (i)
            {
                res.mul(*this);
            }
        }

        return res;
    }

    void frobenius_map(usize power)
    {
        if (!(power == 1 || power == 2 || power == 3 || power == 6))
        {
            unreachable(stringf("can not reach power %u", power));
        }
        this->c0.frobenius_map(power);
        this->c1.frobenius_map(power);

        this->c1
            .c0
            .mul(this->field.frobenius_coeffs_c1[power % 12]);
        this->c1
            .c1
            .mul(this->field.frobenius_coeffs_c1[power % 12]);
        this->c1
            .c2
            .mul(this->field.frobenius_coeffs_c1[power % 12]);
    }

    void mul_by_034(Fp2<N> const &c0, Fp2<N> const &c3, Fp2<N> const &c4)
    {
        auto a = this->c0;
        a.c0.mul(c0);
        a.c1.mul(c0);
        a.c2.mul(c0);

        auto b = this->c1;
        b.mul_by_01(c3, c4);

        auto t0 = c0;
        t0.add(c3);

        auto e = this->c0;
        e.add(this->c1);
        e.mul_by_01(t0, c4);

        this->c1 = e;
        this->c1.sub(a);
        this->c1.sub(b);

        auto t1 = b;
        this->field.mul_by_nonresidue(t1);
        this->c0 = a;
        this->c0.add(t1);
    }

    void mul_by_014(Fp2<N> const &c0, Fp2<N> const &c1, Fp2<N> const &c4)
    {
        auto aa = this->c0;
        aa.mul_by_01(c0, c1);
        auto bb = this->c1;
        bb.mul_by_1(c4);
        auto o = c1;
        o.add(c4);
        this->c1.add(this->c0);
        this->c1.mul_by_01(c0, o);
        this->c1.sub(aa);
        this->c1.sub(bb);
        this->c0 = bb;
        this->field.mul_by_nonresidue(this->c0);
        this->c0.add(aa);
    }

    // ************************* ELEMENT impl ********************************* //

    template <class C>
    static Fp12<N> one(C const &context)
    {
        FieldExtension2over3over2<N> const &field = context;
        return Fp12<N>(Fp6_3<N>::one(context), Fp6_3<N>::zero(context), field);
    }

    template <class C>
    static Fp12<N> zero(C const &context)
    {
        FieldExtension2over3over2<N> const &field = context;
        return Fp12<N>(Fp6_3<N>::zero(context), Fp6_3<N>::zero(context), field);
    }

    Fp12<N> one() const
    {
        return Fp12::one(this->field);
    }

    Fp12<N> zero() const
    {
        return Fp12::zero(this->field);
    }

    Fp12<N> &self()
    {
        return *this;
    }

    Fp12<N> const &self() const
    {
        return *this;
    }

    bool operator==(Fp12<N> const &other) const
    {
        return this->c0 == other.c0 && this->c1 == other.c1;
    }

    bool operator!=(Fp12<N> const &other) const
    {
        return !(*this == other);
    }
};

#endif