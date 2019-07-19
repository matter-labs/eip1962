#ifndef H_FP3
#define H_FP3

#include "../common.h"
#include "../element.h"
#include "../fp.h"
#include "../field.h"
#include "common.h"

using namespace cbn::literals;

template <usize N>
class FieldExtension3 : public PrimeField<N>
{
    Fp<N> _non_residue;

public:
    std::array<Fp<N>, 3> frobenius_coeffs_c1, frobenius_coeffs_c2;

    FieldExtension3(Fp<N> non_residue, PrimeField<N> const &field) : PrimeField<N>(field), _non_residue(non_residue), frobenius_coeffs_c1({Fp<N>::zero(field), Fp<N>::zero(field), Fp<N>::zero(field)}), frobenius_coeffs_c2({Fp<N>::zero(field), Fp<N>::zero(field), Fp<N>::zero(field)})
    {
        // NON_RESIDUE**(((q^0) - 1) / 3)
        auto const f_0 = Fp<N>::one(field);

        // NON_RESIDUE**(((q^1) - 1) / 3)
        auto const f_1 = calc_frobenius_factor(non_residue, field.mod(), 3, "Fp3");

        // NON_RESIDUE**(((q^2) - 1) / 3)
        auto const f_2 = calc_frobenius_factor(non_residue, field.mod() * field.mod(), 3, "Fp3");

        auto const f_0_c2 = f_0;

        auto f_1_c2 = f_1;
        auto f_2_c2 = f_2;

        f_1_c2.square();
        f_2_c2.square();

        std::array<Fp<N>, 3> calc_frobenius_coeffs_c1 = {f_0, f_1, f_2};
        std::array<Fp<N>, 3> calc_frobenius_coeffs_c2 = {f_0_c2, f_1_c2, f_2_c2};

        frobenius_coeffs_c1 = calc_frobenius_coeffs_c1;
        frobenius_coeffs_c2 = calc_frobenius_coeffs_c2;
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
class Fp3 : public Element<Fp3<N>>
{
    FieldExtension3<N> const &field;

public:
    Fp<N> c0, c1, c2;

    Fp3(Fp<N> c0, Fp<N> c1, Fp<N> c2, FieldExtension3<N> const &field) : field(field), c0(c0), c1(c1), c2(c2) {}

    auto operator=(Fp3<N> const &other)
    {
        c0 = other.c0;
        c1 = other.c1;
        c2 = other.c2;
    }

    void mul_by_fp(Fp<N> const &element)
    {
        c0.mul(element);
        c1.mul(element);
        c2.mul(element);
    }

    // ************************* ELEMENT impl ********************************* //

    template <class C>
    static Fp3<N> one(C const &context)
    {
        FieldExtension3<N> const &field = context;
        return Fp3<N>(Fp<N>::one(context), Fp<N>::zero(context), Fp<N>::zero(context), field);
    }

    template <class C>
    static Fp3<N> zero(C const &context)
    {
        FieldExtension3<N> const &field = context;
        return Fp3<N>(Fp<N>::zero(context), Fp<N>::zero(context), Fp<N>::zero(context), field);
    }

    Fp3<N> one() const
    {
        return Fp3::one(field);
    }

    Fp3<N> zero() const
    {
        return Fp3::zero(field);
    }

    Fp3<N> &self()
    {
        return *this;
    }

    Fp3<N> const &self() const
    {
        return *this;
    }

    void serialize(u8 mod_byte_len, std::vector<u8> &data) const
    {
        c0.serialize(mod_byte_len, data);
        c1.serialize(mod_byte_len, data);
        c2.serialize(mod_byte_len, data);
    }

    Option<Fp3<N>> inverse() const
    {
        if (this->is_zero())
        {
            return {};
        }

        auto t0 = this->c0;
        t0.square();
        auto t1 = this->c1;
        t1.square();
        auto t2 = this->c2;
        t2.square();
        auto t3 = this->c0;
        t3.mul(this->c1);
        auto t4 = this->c0;
        t4.mul(this->c2);
        auto t5 = this->c1;
        t5.mul(this->c2);
        auto n5 = t5;
        field.mul_by_nonresidue(n5);

        auto s0 = t0;
        s0.sub(n5);
        auto s1 = t2;
        field.mul_by_nonresidue(s1);
        s1.sub(t3);
        auto s2 = t1;
        s2.sub(t4); // typo in paper referenced above. should be "-" as per Scott, but is "*"

        auto a1 = this->c2;
        a1.mul(s1);
        auto a2 = this->c1;
        a2.mul(s2);
        auto a3 = a1;
        a3.add(a2);
        field.mul_by_nonresidue(a3);
        auto t6 = this->c0;
        t6.mul(s0);
        t6.add(a3);
        auto const o_t6 = t6.inverse();
        if (!o_t6)
        {
            return {};
        }

        t6 = o_t6.value();

        auto x0 = t6;
        x0.mul(s0);
        auto x1 = t6;
        x1.mul(s1);
        auto x2 = t6;
        x2.mul(s2);

        return Fp3(x0, x1, x2, field);
    }

    void square()
    {
        auto const a = c0;
        auto const b = c1;
        auto const c = c2;

        auto s0 = a;
        s0.square();
        auto ab = a;
        ab.mul(b);
        auto s1 = ab;
        s1.mul2();
        auto s2 = a;
        s2.sub(b);
        s2.add(c);
        s2.square();
        auto bc = b;
        bc.mul(c);
        auto s3 = bc;
        s3.mul2();
        auto s4 = c;
        s4.square();

        c0 = s0;
        auto t0 = s3;
        field.mul_by_nonresidue(t0);
        c0.add(t0);

        c1 = s1;
        auto t1 = s4;
        field.mul_by_nonresidue(t1);
        c1.add(t1);

        c2 = s1;
        c2.add(s2);
        c2.add(s3);
        c2.sub(s0);
        c2.sub(s4);
    }

    void mul2()
    {
        c0.mul2();
        c1.mul2();
        c2.mul2();
    }

    void mul(Fp3<N> const &other)
    {
        auto const a = other.c0;
        auto const b = other.c1;
        auto const c = other.c2;

        auto const d = c0;
        auto const e = c1;
        auto const f = c2;

        auto ad = d;
        ad.mul(a);
        auto be = e;
        be.mul(b);
        auto cf = f;
        cf.mul(c);

        auto t0 = b;
        t0.add(c);

        auto x = e;
        x.add(f);
        x.mul(t0);
        x.sub(be);
        x.sub(cf);

        auto t1 = a;
        t1.add(b);

        auto y = d;
        y.add(e);
        y.mul(t1);
        y.sub(ad);
        y.sub(be);

        auto t2 = a;
        t2.add(c);

        auto z = d;
        z.add(f);
        z.mul(t2);
        z.sub(ad);
        z.add(be);
        z.sub(cf);

        auto t3 = x;
        field.mul_by_nonresidue(t3);

        c0 = t3;
        c0.add(ad);

        auto t4 = cf;
        field.mul_by_nonresidue(t4);

        c1 = t4;
        c1.add(y);

        c2 = z;
    }

    void sub(Fp3<N> const &e)
    {
        c0.sub(e.c0);
        c1.sub(e.c1);
        c2.sub(e.c2);
    }

    void add(Fp3<N> const &e)
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
    void frobenius_map(usize power)
    {
        c1.mul(field.frobenius_coeffs_c1[power % 3]);
        c2.mul(field.frobenius_coeffs_c2[power % 3]);
    }

    bool is_zero() const
    {
        return c0.is_zero() && c1.is_zero() && c2.is_zero();
    }

    bool operator==(Fp3<N> const &other) const
    {
        return c0 == other.c0 && c1 == other.c1 && c2 == other.c2;
    }

    bool operator!=(Fp3<N> const &other) const
    {
        return !(*this == other);
    }
};

#endif