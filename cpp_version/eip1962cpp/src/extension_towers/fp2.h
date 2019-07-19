#ifndef H_FP2
#define H_FP2

#include "../common.h"
#include "../element.h"
#include "../fp.h"
#include "../field.h"
#include "common.h"

using namespace cbn::literals;

template <usize N>
class FieldExtension2 : public PrimeField<N>
{
    Fp<N> _non_residue;

public:
    std::array<Fp<N>, 2> frobenius_coeffs_c1;

    FieldExtension2(Fp<N> non_residue, PrimeField<N> const &field) : PrimeField<N>(field), _non_residue(non_residue), frobenius_coeffs_c1({Fp<N>::zero(field), Fp<N>::zero(field)})
    {
        // calculate_frobenius_coeffs

        // NONRESIDUE**(((q^0) - 1) / 2)
        auto const f_0 = Fp<N>::one(field);

        // NONRESIDUE**(((q^1) - 1) / 2)
        auto const f_1 = calc_frobenius_factor(non_residue, field.mod(), 2, "Fp2");

        std::array<Fp<N>, 2> calc_frobenius_coeffs_c1 = {f_0, f_1};
        frobenius_coeffs_c1 = calc_frobenius_coeffs_c1;
    }

    void mul_by_nonresidue(Fp<N> &num) const
    {
        num.mul(_non_residue);
    }

    Fp<N> const &non_residue() const
    {
        return _non_residue;
    }
};

template <usize N>
class Fp2 : public Element<Fp2<N>>
{
public:
    FieldExtension2<N> const &field;
    Fp<N> c0, c1;

    Fp2(Fp<N> c0, Fp<N> c1, FieldExtension2<N> const &field) : field(field), c0(c0), c1(c1) {}

    auto operator=(Fp2<N> const &other)
    {
        c0 = other.c0;
        c1 = other.c1;
    }

    void mul_by_fp(Fp<N> const &element)
    {
        c0.mul(element);
        c1.mul(element);
    }

    // ************************* ELEMENT impl ********************************* //

    template <class C>
    static Fp2<N> one(C const &context)
    {
        FieldExtension2<N> const &field = context;
        return Fp2<N>(Fp<N>::one(context), Fp<N>::zero(context), field);
    }

    template <class C>
    static Fp2<N> zero(C const &context)
    {
        FieldExtension2<N> const &field = context;
        return Fp2<N>(Fp<N>::zero(context), Fp<N>::zero(context), field);
    }

    Fp2<N> one() const
    {
        return Fp2::one(field);
    }

    Fp2<N> zero() const
    {
        return Fp2::zero(field);
    }

    Fp2<N> &self()
    {
        return *this;
    }

    Fp2<N> const &self() const
    {
        return *this;
    }

    void serialize(u8 mod_byte_len, std::vector<u8> &data) const
    {
        c0.serialize(mod_byte_len, data);
        c1.serialize(mod_byte_len, data);
    }

    Option<Fp2<N>> inverse() const
    {
        if (is_zero())
        {
            return {};
        }
        else
        {
            // Guide to Pairing-based Cryptography, Algorithm 5.19.
            // v0 = c0.square()
            auto v0 = c0;
            v0.square();
            // v1 = c1.square()
            auto v1 = c1;
            v1.square();
            // v0 = v0 - beta * v1
            auto v1_by_nonresidue = v1;
            field.mul_by_nonresidue(v1_by_nonresidue);
            v0.sub(v1_by_nonresidue);
            auto o = v0.inverse();
            if (o)
            {
                auto v1 = o.value();
                auto e0 = c0;
                e0.mul(v1);
                auto e1 = c1;
                e1.mul(v1);
                e1.negate();

                return Fp2(e0, e1, field);
            }
            else
            {
                return {};
            }
        }
    }

    void square()
    {
        // v0 = c0 - c1
        auto v0 = c0;
        v0.sub(c1);
        // v3 = c0 - beta * c1
        auto v3 = c0;
        auto t0 = c1;
        field.mul_by_nonresidue(t0);
        v3.sub(t0);
        // v2 = c0 * c1
        auto v2 = c0;
        v2.mul(c1);

        // v0 = (v0 * v3) + v2
        v0.mul(v3);
        v0.add(v2);

        c1 = v2;
        c1.mul2();
        c0 = v0;
        field.mul_by_nonresidue(v2);
        c0.add(v2);
    }

    void mul2()
    {
        c0.mul2();
        c1.mul2();
    }

    void mul(Fp2<N> const &other)
    {
        auto v0 = c0;
        v0.mul(other.c0);
        auto v1 = c1;
        v1.mul(other.c1);

        c1.add(c0);
        auto t0 = other.c0;
        t0.add(other.c1);
        c1.mul(t0);
        c1.sub(v0);
        c1.sub(v1);
        c0 = v0;
        field.mul_by_nonresidue(v1);
        c0.add(v1);
    }

    void sub(Fp2<N> const &e)
    {
        c0.sub(e.c0);
        c1.sub(e.c1);
    }

    void add(Fp2<N> const &e)
    {
        c0.add(e.c0);
        c1.add(e.c1);
    }

    void negate()
    {
        c0.negate();
        c1.negate();
    }

    void frobenius_map(usize power)
    {
        c1.mul(field.frobenius_coeffs_c1[power % 2]);
    }

    bool is_zero() const
    {
        return c0.is_zero() && c1.is_zero();
    }

    bool operator==(Fp2<N> const &other) const
    {
        return c0 == other.c0 && c1 == other.c1;
    }

    bool operator!=(Fp2<N> const &other) const
    {
        return !(*this == other);
    }

    // *************** impl ************ //

    bool is_non_nth_root(u64 n) const
    {
        if (is_zero())
        {
            return false;
        }
        return this->is_non_nth_root_with(n, field.mod() * field.mod());
    }
};

#endif