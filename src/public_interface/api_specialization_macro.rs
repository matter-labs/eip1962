#[macro_export]
macro_rules! expand_for_modulus_limbs {
    ($modulus_limbs: expr, $implementation: tt, $argument: expr, $func: tt) => {
        match $modulus_limbs {
            4 => {
                $implementation::<U256Repr>::$func(&$argument)
            },
            5 => {
                $implementation::<U320Repr>::$func(&$argument)
            },
            6 => {
                $implementation::<U384Repr>::$func(&$argument)
            },
            7 => {
                $implementation::<U448Repr>::$func(&$argument)
            },
            8 => {
                $implementation::<U512Repr>::$func(&$argument)
            },
            9 => {
                $implementation::<U576Repr>::$func(&$argument)
            },
            10 => {
                $implementation::<U640Repr>::$func(&$argument)
            },
            11 => {
                $implementation::<U704Repr>::$func(&$argument)
            },
            12 => {
                $implementation::<U768Repr>::$func($argument)
            },
            13 => {
                $implementation::<U832Repr>::$func(&$argument)
            },
            14 => {
                $implementation::<U896Repr>::$func(&$argument)
            },
            15 => {
                $implementation::<U960Repr>::$func(&$argument)
            },
            16 => {
                $implementation::<U1024Repr>::$func($argument)
            },

            field_limbs => {
                unimplemented!("unimplemented for {} modulus limbs", field_limbs);
            }
        }
    }
}