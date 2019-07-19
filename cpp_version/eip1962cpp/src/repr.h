#ifndef H_REPR
#define H_REPR

#include "common.h"
#include "ctbignum/ctbignum.hpp"

const static usize LIMB_BITS = sizeof(u64) * 8;

template <usize N>
using Repr = cbn::big_int<N>;

// *********************** FUNCTIONS on Repr ******************* //
namespace cbn
{
template <usize N>
bool is_even(Repr<N> const &repr)
{
    return (repr[0] & 0x1) == 0;
}

template <usize N>
Repr<N> mul2(Repr<N> repr)
{
    return cbn::detail::first<N>(cbn::shift_left(repr, 1));
}

template <usize N>
Repr<N> div2(Repr<N> repr)
{
    return cbn::shift_right(repr, 1);
}

template <usize N>
bool is_zero(Repr<N> const &repr)
{
    constexpr Repr<N> zero = {0};
    return repr == zero;
}

} // namespace cbn

// *********************** FUNCTIONS on u64 ******************* //

// Calculate a - b - borrow, returning the result and modifying
// the borrow value.
u64 sbb(u64 a, u64 b, u64 &borrow);

// Calculate a + b + carry, returning the sum and modifying the
// carry value.
u64 adc(u64 a, u64 b, u64 &carry);

/* Function to get no of set bits in binary 
representation of positive integer n */
u32 count_ones(u64 n);

// *********************** FUNCTIONS on std::vector<u64> ******************* //

bool is_zero(std::vector<u64> const &repr);

bool is_odd(std::vector<u64> const &repr);

void div2(std::vector<u64> &repr);

void sub_noborrow(std::vector<u64> &repr, u64 value);

void add_nocarry(std::vector<u64> &repr, u64 value);

void add_scalar(std::vector<u64> &repr, u64 value);

// a >= b
// Where a and b are numbers
bool greater_or_equal(std::vector<u64> const &a, std::vector<u64> const &b);

// Minimal number of bits necessary to represent number in repr
u32 num_bits(std::vector<u64> const &repr);

void mul_scalar(std::vector<u64> &repr, u64 scalar);

void right_shift(std::vector<u64> &repr, u64 shift);

std::vector<i64> into_ternary_wnaf(std::vector<u64> const &repr);

u32 calculate_hamming_weight(std::vector<u64> const &repr);

// ********************** ITERATORS ******************* //
// A is required to have to methods: size() and operator[]->u64
template <class A>
class RevBitIterator
{
    A const &repr;
    usize at;

public:
    // Skips higher zeros
    RevBitIterator(A const &repr) : repr(repr), at(repr.size() * LIMB_BITS)
    {
        // Skip higher zeros
        while (this->before())
        {
            if (**this)
            {
                at++;
                break;
            }
        }
    }

    bool operator*()
    {
        auto i = at / LIMB_BITS;
        auto off = at - (i * LIMB_BITS);
        return (repr[i] >> off) & 0x1;
    }

    /// True if moved
    bool before()
    {
        if (at > 0)
        {
            at--;
            return true;
        }
        else
        {
            return false;
        }
    }
};

// A is required to have to methods: size() and operator[]->u64
template <class A>
class BitIterator
{
    A const &repr;
    usize at;

public:
    BitIterator(A const &repr) : repr(repr), at(0)
    {
    }

    bool operator*()
    {
        auto i = at / LIMB_BITS;
        auto off = at - (i * LIMB_BITS);
        return (repr[i] >> off) & 0x1;
    }

    bool ok() const
    {
        return at < repr.size() * LIMB_BITS;
    }

    void operator++()
    {
        at++;
    }
};

#endif