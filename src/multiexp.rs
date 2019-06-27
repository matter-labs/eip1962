use crate::weierstrass::Group;
use crate::representation::ElementRepr;
use crate::representation::IntoWnaf;

// // start with naive implementation. It'll not work faster than
// // a trivial one before wNAF is implemented
// pub(crate) fn ben_coster<G: Group, E: ElementRepr>(pairs: Vec<(G, E)>) -> G {
//     // sort the pairs together
//     if pairs.len() == 1 {
//         return pairs[0].0.mul(pairs[0].1);
//     }

//     let mut pairs = pairs;

//     pairs.sort_by(|a, b| a.1.cmp(&b.1));
//     let mut reduced = vec![];
//     let num_iter = pairs.len() - 1;
//     for _ in 0..num_iter {
//         // this element has a largest scalar
//         let mut last = pairs.pop().unwrap();
//         // subtract a scalar
//         last.1.sub_noborrow(&pairs.last().unwrap().1);
//         // add a point to the point of the "next scalar"
//         pairs.last_mut().unwrap().0.add_assign_mixed(&last.0);
//         // push the result
//         reduced.push(last);
//     }

//     let acc = pairs.pop().unwrap();
//     let mut acc = acc.0.mul(&acc.1);

//     let mut results: Vec<G> = reduced.into_iter().map(|pair| pair.0.mul(pair.1)).collect();
//     // let mut results: Vec<G> = reduced.into_iter().map(|pair| pair.0.wnaf_mul_impl(pair.1)).collect();

//     while let Some(p) = results.pop() {
//         acc.add_assign(&p);
//     }

//     acc
// }

// // start with naive implementation. It'll not work faster than
// // a trivial one before wNAF is implemented
// pub(crate) fn ben_coster_wnaf<G: Group, E: ElementRepr + IntoWnaf>(pairs: Vec<(G, E)>) -> G {
//     // sort the pairs together
//     if pairs.len() == 1 {
//         return pairs[0].0.mul(pairs[0].1);
//     }

//     let mut pairs = pairs;

//     pairs.sort_by(|a, b| a.1.cmp(&b.1));
//     let mut reduced = vec![];
//     let num_iter = pairs.len() - 1;
//     for _ in 0..num_iter {
//         // this element has a largest scalar
//         let mut last = pairs.pop().unwrap();
//         // subtract a scalar
//         last.1.sub_noborrow(&pairs.last().unwrap().1);
//         // add a point to the point of the "next scalar"
//         pairs.last_mut().unwrap().0.add_assign_mixed(&last.0);
//         // push the result
//         reduced.push(last);
//     }

//     let acc = pairs.pop().unwrap();
//     let mut acc = acc.0.wnaf_mul(acc.1);

//     let mut results: Vec<G> = reduced.into_iter().map(|pair| pair.0.wnaf_mul(pair.1)).collect();

//     while let Some(p) = results.pop() {
//         acc.add_assign(&p);
//     }

//     acc
// }


use crate::weierstrass::curve::CurvePoint;
use crate::field::SizedPrimeField;
pub(crate) fn peppinger<'a, FE: ElementRepr, F: SizedPrimeField<Repr = FE>>
    ( pairs: Vec<(CurvePoint<'a, FE, F>, Vec<u64>)> ) -> CurvePoint<'a, FE, F>
{
    use crate::representation::*;
    let mut g = vec![];
    let mut s: Vec<Vec<u64>> = vec![];

    for (point, scalar) in pairs.into_iter() {
        g.push(point);
        s.push(scalar);
    }

    let c = if s.len() < 32 {
        3u32
    } else {
        (f64::from(s.len() as u32)).ln().ceil() as u32
    };

    let mut windows = vec![];
    let mut buckets = vec![];

    let mask = (1u64 << c) - 1u64;
    let mut cur = 0;
    let num_bits = num_bits(&g[0].curve.subgroup_order_repr);
    let zero_point = CurvePoint::zero(g[0].curve);
    while cur <= num_bits {
        let mut acc = zero_point.clone();

        buckets.truncate(0);
        buckets.resize((1 << c) - 1, zero_point.clone());

        let g = g.clone();

        for (s, g) in s.iter_mut().zip(g) {
            let index = (s[0] & mask) as usize;

            if index != 0 {
                buckets[index - 1].add_assign_mixed(&g);
            }

            left_shift_representation(s, c as u64);

            // s.shr(c as u32);
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