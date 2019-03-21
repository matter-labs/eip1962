use crate::weierstrass::Group;
use crate::representation::ElementRepr;

// find what will work fast enough

// start with naive implementation
pub fn ben_coster<G: Group, E: ElementRepr>(pairs: Vec<(G, E)>) -> G {
    // sort the pairs together
    if pairs.len() == 1 {
        return pairs[0].0.mul(pairs[0].1);
    }

    let mut pairs = pairs;

    // pairs.sort_by(|a, b| b.1.cmp(&a.1));
    pairs.sort_by(|a, b| a.1.cmp(&b.1));
    let mut reduced = vec![];
    let num_iter = pairs.len() - 1;
    for _ in 0..num_iter {
        // this element has a largest scalar
        let mut last = pairs.pop().unwrap();
        // subtract a scalar
        last.1.sub_noborrow(&pairs.last().unwrap().1);
        // add a point to the point of the "next scalar"
        pairs.last_mut().unwrap().0.add_assign_mixed(&last.0);
        // push the result
        reduced.push(last);
    }

    let acc = pairs.pop().unwrap();
    let mut acc = acc.0.mul(&acc.1);

    let mut results: Vec<G> = reduced.into_iter().map(|pair| pair.0.mul(pair.1)).collect();

    while let Some(p) = results.pop() {
        acc.add_assign(&p);
    }

    acc
}