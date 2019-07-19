#ifndef H_MNT6engine
#define H_MNT6engine

#include "../common.h"
#include "../curve.h"
#include "../fp.h"
#include "../extension_towers/fp3.h"
#include "../extension_towers/fp6_2.h"
#include "mnt.h"

template <usize N>
class MNT6engine : public MNTengine<Fp3<N>, Fp6_2<N>, N>
{
public:
    MNT6engine(std::vector<u64> x, bool x_is_negative, std::vector<u64> exp_w0, std::vector<u64> exp_w1, bool exp_w0_is_negative,
               WeierstrassCurve<Fp3<N>> const &curve_twist, Fp3<N> twist) : MNTengine<Fp3<N>, Fp6_2<N>, N>(x, x_is_negative, exp_w0, exp_w1, exp_w0_is_negative, curve_twist, twist)
    {
    }

    Fp6_2<N> final_exponentiation_part_one(Fp6_2<N> const &elt, Fp6_2<N> const &elt_inv) const
    {
        // (q^3-1)*(q+1)

        // elt_q3 = elt^(q^3)
        auto elt_q3 = elt;
        elt_q3.frobenius_map(3);
        // elt_q3_over_elt = elt^(q^3-1)
        auto elt_q3_over_elt = elt_q3;
        elt_q3_over_elt.mul(elt_inv);
        // alpha = elt^((q^3-1) * q)
        auto alpha = elt_q3_over_elt;
        alpha.frobenius_map(1);
        // beta = elt^((q^3-1)*(q+1)
        alpha.mul(elt_q3_over_elt);

        return alpha;
    }
};

#endif