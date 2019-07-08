use crate::errors::ApiError;

pub(crate) fn calculate_hamming_weight(representation: &[u64]) -> u32 {
    let mut weight = 0;
    for el in representation.iter() {
        weight += el.count_ones();
    }

    weight
}

pub(crate) fn generate_powers(parameter: u64, max_power: usize) -> Result<Vec<u64>, ApiError> {
    let mut powers = vec![];
    powers.push(1u64);
    let mut p = 1u64;
    for _ in 1..=max_power {
        p = p.checked_mul(parameter).ok_or(ApiError::Overflow)?;
        powers.push(p);
    }

    Ok(powers)
}

pub(crate) fn evaluate_poly(coeffs: &[u64], parameter: u64) -> Result<u64, ApiError> {
    if coeffs.len() == 0 {
        return Err(ApiError::InputError(format!("can not evaluate empty combination, file {}, line {}", file!(), line!())));
    }
    let mut result = 0u64;
    result += coeffs[0];
    let mut p = 1u64;
    for i in 1..coeffs.len() {
        p = p.checked_mul(parameter).ok_or(ApiError::Overflow)?;
        let term = p.checked_mul(coeffs[i]).ok_or(ApiError::Overflow)?;
        result = result.checked_add(term).ok_or(ApiError::Overflow)?;
    }

    Ok(result)
}