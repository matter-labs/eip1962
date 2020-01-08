pub(crate) fn in_gas_metering() -> bool {
    #[cfg(feature = "gas_metering_mode")]
    return true;

    return crate::features::in_gas_metering();
}

pub(crate) fn in_fuzzing() -> bool {
    #[cfg(feature = "fuzzing_mode")]
    return true;

    return false;
}