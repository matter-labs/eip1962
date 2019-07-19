#ifndef H_FP
#define H_FP

#include "common.h"
#include "element.h"
#include "repr.h"
#include "field.h"

using namespace cbn::literals;

template <usize N>
class Fp : public Element<Fp<N>>
{
    PrimeField<N> const &field;
    Repr<N> repr;

public:
    Fp(Repr<N> repr, PrimeField<N> const &field) : field(field), repr(repr) {}

    static Fp<N> from_repr(Repr<N> repr, PrimeField<N> const &field)
    {
        auto fpo = Fp::from_repr_try(repr, field);
        if (fpo)
        {
            return fpo.value();
        }
        else
        {
            api_err("not an element of the field");
        }
    }

    Fp(Fp<N> const &other) : Fp(other.repr, other.field) {}

    auto operator=(Fp<N> const &other)
    {
        this->repr = other.repr;
    }

    Repr<N> const &representation() const
    {
        return repr;
    }

    ~Fp() {}

    // ************************* ELEMENT impl ********************************* //
    template <class C>
    static Fp<N> one(C const &context)
    {
        PrimeField<N> const &field = context;
        return Fp(field.mont_r(), field);
    }

    template <class C>
    static Fp<N> zero(C const &context)
    {
        constexpr Repr<N> zero = {0};
        PrimeField<N> const &field = context;
        return Fp(zero, field);
    }

    Fp<N> one() const
    {
        return Fp::one(field);
    }

    Fp<N> zero() const
    {
        return Fp::zero(field);
    }

    Fp<N> &self()
    {
        return *this;
    }

    Fp<N> const &self() const
    {
        return *this;
    }

    // Serializes bytes from number to BigEndian u8 format.
    void serialize(u8 mod_byte_len, std::vector<u8> &data) const
    {
        auto const normal_repr = into_repr();
        for (i32 i = i32(mod_byte_len) - 1; i >= 0; i--)
        {
            auto const j = i / sizeof(u64);
            if (j < N)
            {

                auto const off = (i - j * sizeof(u64)) * 8;
                data.push_back(normal_repr[j] >> off);
            }
            else
            {
                data.push_back(0);
            }
        }
    }

    Option<Fp<N>> inverse() const
    {
        return mont_inverse();
    }

    void square()
    {
        repr = cbn::montgomery_mul(repr, repr, field.mod(), field.mont_inv());
    }

    void mul2()
    {
        repr = cbn::shift_left(repr, 1) % field.mod();
    }

    void mul(Fp<N> const &e)
    {
        repr = cbn::montgomery_mul(repr, e.repr, field.mod(), field.mont_inv());
    }

    void sub(Fp<N> const &e)
    {
        repr = cbn::mod_sub(repr, e.repr, field.mod());
    }

    void add(Fp<N> const &e)
    {
        repr = cbn::mod_add(repr, e.repr, field.mod());
    }

    void negate()
    {
        if (!is_zero())
        {
            repr = cbn::subtract_ignore_carry(field.mod(), repr);
        }
    }

    bool is_zero() const
    {
        return cbn::is_zero(repr);
    }

    bool operator==(Fp<N> const &other) const
    {
        return repr == other.repr;
    }

    bool operator!=(Fp<N> const &other) const
    {
        return repr != other.repr;
    }

    // *************** impl ************ //

    bool is_non_nth_root(u64 n) const
    {
        return this->is_non_nth_root_with(n, field.mod());
    }

private:
    static Option<Fp<N>> from_repr_try(Repr<N> repr, PrimeField<N> const &field)
    {
        if (field.is_valid(repr))
        {
            Fp<N> r1 = Fp(repr, field);
            Fp<N> r2 = Fp(field.mont_r2(), field);

            r1.mul(r2);

            return r1;
        }
        else
        {
            return {};
        }
    }

    Repr<N> into_repr() const
    {
        return cbn::montgomery_reduction(cbn::detail::pad<N>(repr), field.mod(), field.mont_inv());
    }

    Option<Fp<N>> mont_inverse() const
    {
        if (is_zero())
        {
            return {};
        }

        // The Montgomery Modular Inverse - Revisited

        // Phase 1
        auto const modulus = field.mod();
        auto u = modulus;
        auto v = repr;
        Repr<N> r = {0};
        Repr<N> s = {1};
        u64 k = 0;

        auto found = false;
        for (usize i = 0; i < N * 128; i++)
        {
            if (cbn::is_zero(v))
            {
                found = true;
                break;
            }
            if (cbn::is_even(u))
            {
                u = cbn::div2(u);
                s = cbn::mul2(s);
            }
            else if (cbn::is_even(v))
            {
                v = cbn::div2(v);
                r = cbn::mul2(r);
            }
            else if (u > v)
            {
                u = cbn::subtract_ignore_carry(u, v);
                u = cbn::div2(u);
                r = cbn::add_ignore_carry(r, s);
                s = cbn::mul2(s);
            }
            else if (v >= u)
            {
                v = cbn::subtract_ignore_carry(v, u);
                v = cbn::div2(v);
                s = cbn::add_ignore_carry(s, r);
                r = cbn::mul2(r);
            }

            k += 1;
        }

        if (!found)
        {
            return {};
        }

        if (r >= modulus)
        {
            r = cbn::subtract_ignore_carry(r, modulus);
        }

        r = cbn::subtract_ignore_carry(modulus, r);

        // phase 2

        auto const mont_power_param = field.mont_power();
        if (k < mont_power_param)
        {
            return {};
        }

        for (usize i = 0; i < (k - mont_power_param); i++)
        {
            if (cbn::is_even(r))
            {
                r = cbn::div2(r);
            }
            else
            {
                r = cbn::add_ignore_carry(r, modulus);
                r = cbn::div2(r);
            }
        }

        auto const el = Fp::from_repr_try(r, field);
        if (el)
        {
            return el.value();
        }
        else
        {
            return {};
        }
    }
};

#endif