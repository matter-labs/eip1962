#ifndef H_FROBENIUS
#define H_FROBENIUS

#include "../common.h"
#include "../element.h"
#include "../fp.h"
#include "../field.h"

template <usize M>
Repr<M + 1> calc_frobenius_factor_2(Repr<M> const &base, usize div, std::string const &err)
{
    // NON_RESIDUE**(((base) - 1) / div)
    auto const q_power = base;
    constexpr Repr<M> one = {1};
    auto power = q_power - one;
    Repr<M> const rdiv = {div};
    if (!cbn::is_zero(power % rdiv))
    {
        unknown_parameter_err("Failed to calculate Frobenius coeffs for " + err);
    }
    power = power / rdiv;
    return power;
}

template <usize N, usize M>
Fp<N> calc_frobenius_factor(Fp<N> const &non_residue, Repr<M> const &base, usize div, std::string const &err)
{
    auto const f = non_residue.pow(calc_frobenius_factor_2(base, div, err));
    return f;
}

template <class F>
std::vector<F> calculate_window_table(F base, usize window)
{
    std::vector<F> table;
    table.reserve(1 << (window - 1));

    auto acc = base;
    table.push_back(acc);
    auto square = acc;
    square.square();

    // pushed 1*G, 3*G, 5*G, etc (notation abuse, it's actually exp)
    for (auto i = 1; i < (1 << (window - 1)); i++)
    {
        acc.mul(square);
        table.push_back(acc);
    }

    return table;
}

template <class E>
class WindowExpBase
{
    u32 window_size;
    E one;
    std::vector<E> bases;

public:
    WindowExpBase(E base, E one, usize window) : window_size(window), one(one)
    {
        bases = calculate_window_table(base, window);
    }

    template <usize N>
    E exponentiate(Repr<N> scalar) const
    {
        auto const wnaf = windows(scalar);

        auto res = one;
        auto found_nonzero = false;

        for (auto it = wnaf.crbegin(); it != wnaf.crend(); it++)
        {
            auto const w = *it;
            if (w == 0 && found_nonzero)
            {
                res.square();
            }
            else if (w != 0)
            {
                found_nonzero = true;
                for (usize i = 0; i < window_size; i++)
                {
                    res.square();
                }
                usize const idx = w >> 1;
                res.mul(bases[idx]);
            }
        }

        return res;
    }

    template <usize N>
    std::vector<u64> windows(Repr<N> scalar) const
    {
        std::vector<u64> result;
        auto found_begining = false;
        u64 w = 0;
        u64 bit_count = 0;
        for (auto it = BitIterator(scalar); it.ok(); ++it)
        {
            auto const b = *it;
            if (b)
            {
                if (found_begining)
                {
                    w |= u64(1) << bit_count;
                    bit_count += 1;
                }
                else
                {
                    found_begining = true;
                    w |= u64(1) << bit_count;
                    bit_count += 1;
                }
            }
            else
            {
                if (found_begining)
                {
                    bit_count += 1;
                }
                else
                {
                    result.push_back(0);
                    continue;
                }
            }
            if (found_begining && bit_count == window_size)
            {
                result.push_back(w);
                w = 0;
                found_begining = false;
                bit_count = 0;
            }
        }

        if (w != 0)
        {
            // this is a last chunk if bit length is not divisible by window size
            result.push_back(w);
        }

        auto const n = result.size();
        for (usize i = 0; i < n; i++)
        {
            if (result.size() > 0)
            {
                if (result.back() == 0)
                {
                    result.pop_back();
                }
                else
                {
                    break;
                }
            }
        }

        return result;
    }
};

#endif