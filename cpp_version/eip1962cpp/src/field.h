#ifndef H_FIELD
#define H_FIELD

#include "common.h"
#include "repr.h"
#include "ctbignum/slicing.hpp"

template <usize N>
class PrimeField
{
    Repr<N> modulus;
    u64 mont_power_;
    Repr<N> mont_r_;
    Repr<N> mont_r2_;
    u64 mont_inv_;

public:
    PrimeField(Repr<N> modulus) : modulus(modulus), mont_power_(N * LIMB_BITS)
    {
        // Compute -m^-1 mod 2**64 by exponentiating by totient(2**64) - 1
        u64 inv = 1;
        for (auto i = 0; i < 63; i++)
        {
            inv = inv * inv;
            inv = inv * modulus[0];
        }
        inv = (std::numeric_limits<u64>::max() - inv) + 2 + std::numeric_limits<u64>::max();
        mont_inv_ = inv;

        Repr<N + 1> pow_N_LIMB_BITS = {0};
        pow_N_LIMB_BITS[N] = 1;
        mont_r_ = pow_N_LIMB_BITS % modulus;

        mont_r2_ = (mont_r_ * mont_r_) % modulus;
    }

    Repr<N> mod() const
    {
        return modulus;
    }

    Repr<N> mont_r() const
    {
        return mont_r_;
    }

    Repr<N> mont_r2() const
    {
        return mont_r2_;
    }

    u64 mont_power() const
    {
        return mont_power_;
    }

    // Montgomery parametare for multiplication
    u64 mont_inv() const
    {
        return mont_inv_;
    }

    bool is_valid(Repr<N> const &repr) const
    {
        return repr < modulus;
    }
};

#endif