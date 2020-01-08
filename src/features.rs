pub(crate) fn in_gas_metering() -> bool {
    #[cfg(feature = "gas_metering_mode")]
    return true;

    return std::option_env!("GAS_METERING").is_some();
}

pub(crate) fn in_fuzzing() -> bool {
    #[cfg(feature = "fuzzing_mode")]
    return true;

    return false;
}

pub(crate) fn in_fuzzing_or_gas_metering() -> bool {
    in_fuzzing() || in_gas_metering()
}