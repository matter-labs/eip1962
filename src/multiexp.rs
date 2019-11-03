use crate::weierstrass::Group;
use crate::weierstrass::curve::CurvePoint;
use crate::weierstrass::CurveParameters;
use crate::constants::MaxGroupSizeUint;

pub(crate) fn peppinger<'a, C: CurveParameters>
    (bases: &[CurvePoint<'a, C>], mut scalars: Vec<MaxGroupSizeUint>) -> CurvePoint<'a, C>
{
    use crate::representation::*;
    debug_assert!(bases.len() == scalars.len());

    let c = if bases.len() < 32 {
        3u32
    } else {
        (f64::from(bases.len() as u32)).ln().ceil() as u32
    };

    let mask = (1u64 << c) - 1u64;
    let mut cur = 0;
    let num_bits = num_bits(&bases[0].curve.subgroup_order_repr);
    let zero_point = CurvePoint::zero(bases[0].curve);

    let mut windows = Vec::with_capacity((num_bits / c + 1) as usize);
    let mut buckets = Vec::with_capacity((1 << c) - 1);

    while cur <= num_bits {
        let mut acc = zero_point.clone();

        buckets.truncate(0);
        buckets.resize((1 << c) - 1, zero_point.clone());

        for (s, g) in scalars.iter_mut().zip(bases.iter()) {
            let index = (s.as_ref()[0] & mask) as usize;

            if index != 0 {
                buckets[index - 1].add_assign_mixed(&g);
            }

            *s >>= c;

            // right_shift_representation(s, c as u64);
        }

        let mut running_sum = zero_point.clone();
        for exp in buckets.iter().rev() {
            running_sum.add_assign(exp);
            acc.add_assign(&running_sum);
        }

        windows.push(acc);

        cur += c;
    }

    let mut acc = zero_point.clone();

    for window in windows.into_iter().rev() {
        for _ in 0..c {
            acc.double();
        }

        acc.add_assign(&window);
    }

    acc
}