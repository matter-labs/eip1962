#ifndef H_MULTIEXP
#define H_MULTIEXP

#include "curve.h"
#include "common.h"

template <class E, class C>
CurvePoint<E> peepinger(std::vector<std::tuple<CurvePoint<E>, std::vector<u64>>> pairs, WeierstrassCurve<E> const &wc, C const &context)
{
    u32 c;
    if (pairs.size() < 32)
    {
        c = 3;
    }
    else
    {
        c = ceil(log((double)pairs.size()));
    };

    std::vector<CurvePoint<E>> windows;
    std::vector<CurvePoint<E>> buckets;

    u64 mask = (u64(1) << c) - u64(1);
    u32 cur = 0;
    auto const n_bits = num_bits(wc.subgroup_order());
    auto const zero_point = CurvePoint<E>::zero(context);

    while (cur <= n_bits)
    {
        auto acc = zero_point;

        buckets.resize(0, zero_point);
        buckets.resize((1 << c) - 1, zero_point);

        for (auto it = pairs.begin(); it != pairs.end(); it++)
        {
            CurvePoint<E> const &g = std::get<0>(*it);
            std::vector<u64> &s = std::get<1>(*it);
            usize const index = s[0] & mask;

            if (index != 0)
            {
                buckets[index - 1].add_mixed(g, wc, context);
            }

            right_shift(s, c);
        }

        auto running_sum = zero_point;
        for (auto it = buckets.crbegin(); it != buckets.crend(); it++)
        {
            running_sum.add(*it, wc, context);
            acc.add(running_sum, wc, context);
        }

        windows.push_back(acc);

        cur += c;
    }

    auto acc = zero_point;

    for (auto it = windows.crbegin(); it != windows.crend(); it++)
    {
        for (u32 i = 0; i < c; i++)
        {
            acc.mul2(wc);
        }

        acc.add(*it, wc, context);
    }

    return acc;
}

#endif
