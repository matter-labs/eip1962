#ifndef H_ELEMENT
#define H_ELEMENT

#include "repr.h"

// Compute unit a layer above representation
template <class E>
class Element
{
public:
    template <class C>
    static E one(C const &context);

    template <class C>
    static E zero(C const &context);

    virtual E one() const = 0;

    virtual E zero() const = 0;

    virtual E &self() = 0;

    virtual E const &self() const = 0;

    virtual void serialize(u8 mod_byte_len, std::vector<u8> &data) const = 0;

    // Computes the multiplicative inverse of this element, if nonzero.
    virtual Option<E> inverse() const = 0;

    virtual void square() = 0;

    virtual void mul2() = 0;

    virtual void mul(E const &e) = 0;

    virtual void sub(E const &e) = 0;

    virtual void add(E const &e) = 0;

    virtual void negate() = 0;

    virtual bool is_zero() const = 0;

    virtual bool operator==(E const &other) const = 0;

    virtual bool operator!=(E const &other) const = 0;

    // ********************* DEFAULT IMPLEMENTED ***************** //
    template <usize N>
    E pow(Repr<N> const &e) const
    {
        auto res = this->one();
        auto found_one = false;

        for (auto it = RevBitIterator(e); it.before();)
        {
            auto i = *it;
            if (found_one)
            {
                res.square();
            }
            else
            {
                found_one = i;
            }

            if (i)
            {
                res.mul(self());
            }
        }

        return res;
    }

    template <usize N>
    bool is_non_nth_root_with(u64 n, Repr<N> power) const
    {
        if (is_zero())
        {
            return false;
        }

        constexpr Repr<N> one = {1};
        power = cbn::subtract_ignore_carry(power, one);
        Repr<N> divisor = {n};
        if (!cbn::is_zero(power % divisor))
        {
            return false;
        }
        power = power / divisor;

        auto l = this->pow(power);
        auto e_one = this->one();

        return l != e_one;
    }
};

#endif