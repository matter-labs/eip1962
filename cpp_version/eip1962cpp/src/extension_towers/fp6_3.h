#ifndef H_FP6_3
#define H_FP6_3

#include "../common.h"
#include "../element.h"
#include "fp2.h"
#include "fp3.h"
#include "../field.h"
#include "common.h"

template <usize N>
class FieldExtension3over2 : public FieldExtension2<N>
{
    Fp2<N> _non_residue;

public:
    std::array<Fp2<N>, 6> frobenius_coeffs_c1, frobenius_coeffs_c2;
    FieldExtension3over2(Fp2<N> _non_residue, FieldExtension2<N> const &field, WindowExpBase<Fp2<N>> const &base) : FieldExtension2<N>(field), _non_residue(_non_residue), frobenius_coeffs_c1({Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field)}), frobenius_coeffs_c2({Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field), Fp2<N>::zero(field)})
    {

        // NON_RESIDUE**(((q^0) - 1) / 3)
        auto const f_0 = Fp2<N>::one(field);

        auto const q_power_1 = field.mod();
        auto const f_1 = base.exponentiate(calc_frobenius_factor_2(q_power_1, 3, "Fp6"));

        auto const q_power_2 = q_power_1 * field.mod();
        auto const f_2 = base.exponentiate(calc_frobenius_factor_2(q_power_2, 3, "Fp6"));

        auto const q_power_3 = q_power_2 * field.mod();
        auto const f_3 = base.exponentiate(calc_frobenius_factor_2(q_power_3, 3, "Fp6"));

        auto const f_4 = Fp2<N>::zero(field);
        auto const f_5 = Fp2<N>::zero(field);

        auto const f_0_c2 = f_0;

        auto f_1_c2 = f_1;
        f_1_c2.square();
        auto f_2_c2 = f_2;
        f_2_c2.square();
        auto f_3_c2 = f_3;
        f_3_c2.square();

        auto const f_4_c2 = f_4;
        auto const f_5_c2 = f_5;

        std::array<Fp2<N>, 6> calc_frobenius_coeffs_c1 = {f_0, f_1, f_2, f_3, f_4, f_5};
        frobenius_coeffs_c1 = calc_frobenius_coeffs_c1;
        std::array<Fp2<N>, 6> calc_frobenius_coeffs_c2 = {f_0_c2, f_1_c2, f_2_c2, f_3_c2, f_4_c2, f_5_c2};
        frobenius_coeffs_c2 = calc_frobenius_coeffs_c2;
    }

    void mul_by_nonresidue(Fp2<N> &num) const
    {
        num.mul(_non_residue);
    }

    Fp2<N> const &non_residue() const
    {
        return _non_residue;
    }
};

template <usize N>
class Fp6_3 : public Element<Fp6_3<N>>
{

public:
    FieldExtension3over2<N> const &field;
    Fp2<N> c0, c1, c2;

    Fp6_3(Fp2<N> c0, Fp2<N> c1, Fp2<N> c2, FieldExtension3over2<N> const &field) : field(field), c0(c0), c1(c1), c2(c2) {}

    auto operator=(Fp6_3<N> const &other)
    {
        c0 = other.c0;
        c1 = other.c1;
        c2 = other.c2;
    }

    void frobenius_map(usize power)
    {
        if (!(power == 0 || power == 1 || power == 2 || power == 3 || power == 6))
        {
            unreachable(stringf("can not reach power %u", power));
        }

        c0.frobenius_map(power);
        c1.frobenius_map(power);
        c2.frobenius_map(power);

        c1.mul(field.frobenius_coeffs_c1[power % 6]);
        c2.mul(field.frobenius_coeffs_c2[power % 6]);
    }

    void mul_by_1(Fp2<N> const &c1)
    {
        auto b_b = this->c1;
        b_b.mul(c1);

        auto t1 = c1;
        {
            auto tmp = this->c1;
            tmp.add(this->c2);

            t1.mul(tmp);
            t1.sub(b_b);
            field.mul_by_nonresidue(t1);
        }

        auto t2 = c1;
        {
            auto tmp = this->c0;
            tmp.add(this->c1);

            t2.mul(tmp);
            t2.sub(b_b);
        }

        this->c0 = t1;
        this->c1 = t2;
        this->c2 = b_b;
    }

    void mul_by_01(Fp2<N> const &c0, Fp2<N> const &c1)
    {
        auto a_a = this->c0;
        auto b_b = this->c1;
        a_a.mul(c0);
        b_b.mul(c1);

        auto t1 = c1;
        {
            auto tmp = this->c1;
            tmp.add(this->c2);

            t1.mul(tmp);
            t1.sub(b_b);
            field.mul_by_nonresidue(t1);
            t1.add(a_a);
        }

        auto t3 = c0;
        {
            auto tmp = this->c0;
            tmp.add(this->c2);

            t3.mul(tmp);
            t3.sub(a_a);
            t3.add(b_b);
        }

        auto t2 = c0;
        t2.add(c1);
        {
            auto tmp = this->c0;
            tmp.add(this->c1);

            t2.mul(tmp);
            t2.sub(a_a);
            t2.sub(b_b);
        }

        this->c0 = t1;
        this->c1 = t2;
        this->c2 = t3;
    }

    // ************************* ELEMENT impl ********************************* //

    template <class C>
    static Fp6_3<N> one(C const &context)
    {
        FieldExtension3over2<N> const &field = context;
        return Fp6_3<N>(Fp2<N>::one(context), Fp2<N>::zero(context), Fp2<N>::zero(context), field);
    }

    template <class C>
    static Fp6_3<N> zero(C const &context)
    {
        FieldExtension3over2<N> const &field = context;
        return Fp6_3<N>(Fp2<N>::zero(context), Fp2<N>::zero(context), Fp2<N>::zero(context), field);
    }

    Fp6_3<N> one() const
    {
        return Fp6_3::one(field);
    }

    Fp6_3<N> zero() const
    {
        return Fp6_3::zero(field);
    }

    Fp6_3<N> &self()
    {
        return *this;
    }

    Fp6_3<N> const &self() const
    {
        return *this;
    }

    void serialize(u8 mod_byte_len, std::vector<u8> &data) const
    {
        c0.serialize(mod_byte_len, data);
        c1.serialize(mod_byte_len, data);
        c2.serialize(mod_byte_len, data);
    }

    // Computes the multiplicative inverse of this element, if nonzero.
    std::optional<Fp6_3<N>> inverse() const
    {
        auto e0 = this->c2;
        field.mul_by_nonresidue(e0);
        e0.mul(this->c1);
        e0.negate();
        {
            auto e0s = this->c0;
            e0s.square();
            e0.add(e0s);
        }
        auto e1 = this->c2;
        e1.square();
        field.mul_by_nonresidue(e1);
        {
            auto e01 = this->c0;
            e01.mul(this->c1);
            e1.sub(e01);
        }
        auto e2 = this->c1;
        e2.square();
        {
            auto e02 = this->c0;
            e02.mul(this->c2);
            e2.sub(e02);
        }

        auto tmp1 = this->c2;
        tmp1.mul(e1);
        auto tmp2 = this->c1;
        tmp2.mul(e2);
        tmp1.add(tmp2);
        field.mul_by_nonresidue(tmp1);
        tmp2 = this->c0;
        tmp2.mul(e0);
        tmp1.add(tmp2);

        if (auto const t = tmp1.inverse())
        {
            auto tmp = Fp6_3<N>(t.value(), t.value(), t.value(), field);
            tmp.c0.mul(e0);
            tmp.c1.mul(e1);
            tmp.c2.mul(e2);
            return tmp;
        }
        else
        {
            return {};
        }
    }

    void square()
    {
        auto s0 = this->c0;
        s0.square();
        auto ab = this->c0;
        ab.mul(this->c1);
        auto s1 = ab;
        s1.mul2();
        auto s2 = this->c0;
        s2.sub(this->c1);
        s2.add(this->c2);
        s2.square();
        auto bc = this->c1;
        bc.mul(this->c2);
        auto s3 = bc;
        s3.mul2();
        auto s4 = this->c2;
        s4.square();

        this->c0 = s3;
        field.mul_by_nonresidue(this->c0);
        this->c0.add(s0);

        this->c1 = s4;
        field.mul_by_nonresidue(this->c1);
        this->c1.add(s1);

        this->c2 = s1;
        this->c2.add(s2);
        this->c2.add(s3);
        this->c2.sub(s0);
        this->c2.sub(s4);
    }

    void mul2()
    {
        c0.mul2();
        c1.mul2();
        c2.mul2();
    }

    void mul(Fp6_3<N> const &other)
    {
        auto a_a = this->c0;
        auto b_b = this->c1;
        auto c_c = this->c2;
        a_a.mul(other.c0);
        b_b.mul(other.c1);
        c_c.mul(other.c2);

        auto t1 = other.c1;
        t1.add(other.c2);
        {
            auto tmp = this->c1;
            tmp.add(this->c2);

            t1.mul(tmp);
            t1.sub(b_b);
            t1.sub(c_c);
            field.mul_by_nonresidue(t1);
            t1.add(a_a);
        }

        auto t3 = other.c0;
        t3.add(other.c2);
        {
            auto tmp = this->c0;
            tmp.add(this->c2);

            t3.mul(tmp);
            t3.sub(a_a);
            t3.add(b_b);
            t3.sub(c_c);
        }

        auto t2 = other.c0;
        t2.add(other.c1);
        {
            auto tmp = this->c0;
            tmp.add(this->c1);

            t2.mul(tmp);
            t2.sub(a_a);
            t2.sub(b_b);
            field.mul_by_nonresidue(c_c);
            t2.add(c_c);
        }

        this->c0 = t1;
        this->c1 = t2;
        this->c2 = t3;
    }

    void sub(Fp6_3<N> const &e)
    {
        c0.sub(e.c0);
        c1.sub(e.c1);
        c2.sub(e.c2);
    }

    void add(Fp6_3<N> const &e)
    {
        c0.add(e.c0);
        c1.add(e.c1);
        c2.add(e.c2);
    }

    void negate()
    {
        c0.negate();
        c1.negate();
        c2.negate();
    }

    bool is_zero() const
    {
        return c0.is_zero() && c1.is_zero() && c2.is_zero();
    }

    bool operator==(Fp6_3<N> const &other) const
    {
        return c0 == other.c0 && c1 == other.c1 && c2 == other.c2;
    }

    bool operator!=(Fp6_3<N> const &other) const
    {
        return !(*this == other);
    }
};

#endif