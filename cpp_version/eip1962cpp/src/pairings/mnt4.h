#ifndef H_MNT4engine
#define H_MNT4engine

#include "../common.h"
#include "../curve.h"
#include "../fp.h"
#include "../extension_towers/fp2.h"
#include "../extension_towers/fp4.h"
#include "mnt.h"

template <usize N>
class MNT4engine : public MNTengine<Fp2<N>, Fp4<N>, N>
{
public:
    MNT4engine(std::vector<u64> x, bool x_is_negative, std::vector<u64> exp_w0, std::vector<u64> exp_w1, bool exp_w0_is_negative,
               WeierstrassCurve<Fp2<N>> const &curve_twist, Fp2<N> twist) : MNTengine<Fp2<N>, Fp4<N>, N>(x, x_is_negative, exp_w0, exp_w1, exp_w0_is_negative, curve_twist, twist)
    {
    }

    Fp4<N> final_exponentiation_part_one(Fp4<N> const &elt, Fp4<N> const &elt_inv) const
    {
        /* (q^2-1) */

        /* elt_q2 = elt^(q^2) */
        auto elt_q2 = elt;
        elt_q2.frobenius_map(2);
        /* elt_q2_over_elt = elt^(q^2-1) */
        auto elt_q2_over_elt = elt_q2;
        elt_q2_over_elt.mul(elt_inv);

        return elt_q2_over_elt;
    }
};

#endif