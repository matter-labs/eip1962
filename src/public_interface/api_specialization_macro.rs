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

// macro_rules! expand_for_modulus_and_group_limbs {
//     ($modulus_limbs: expr, $group_limbs: expr, $implementation: tt, $argument: expr, $func: tt) => {
//         match $modulus_limbs {
//             4 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U256Repr, $argument, $func)
//             },
//             5 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U320Repr, $argument, $func)
//             },
//             6 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U384Repr, $argument, $func)
//             },
//             7 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U448Repr, $argument, $func)
//             },
//             8 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U512Repr, $argument, $func)
//             },
//             9 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U576Repr, $argument, $func)
//             },
//             10 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U640Repr, $argument, $func)
//             },
//             11 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U704Repr, $argument, $func)
//             },
//             12 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U768Repr, $argument, $func)
//             },
//             13 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U832Repr, $argument, $func)
//             },
//             14 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U896Repr, $argument, $func)
//             },
//             15 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U960Repr, $argument, $func)
//             },
//             16 => {
//                 expand_for_group_limbs!($group_limbs, $implementation, U1024Repr, $argument, $func)
//             },

//             field_limbs => {
//                 unimplemented!("unimplemented for {} modulus limbs", field_limbs);
//             }
//         }
//     }
// }

// macro_rules! expand_for_group_limbs {
//     ($limbs: expr, $implementation: tt, $modulus_limbs: tt,  $argument: expr, $func: tt) => {
//         match $limbs {
//             4 => {
//                 $implementation::<$modulus_limbs, U256Repr>::$func(&$argument)
//             },
//             5 => {
//                 $implementation::<$modulus_limbs, U320Repr>::$func(&$argument)
//             },
//             6 => {
//                 $implementation::<$modulus_limbs, U384Repr>::$func(&$argument)
//             },
//             7 => {
//                 $implementation::<$modulus_limbs, U448Repr>::$func(&$argument)
//             },
//             8 => {
//                 $implementation::<$modulus_limbs, U512Repr>::$func(&$argument)
//             },
//             9 => {
//                 $implementation::<$modulus_limbs, U576Repr>::$func(&$argument)
//             },
//             10 => {
//                 $implementation::<$modulus_limbs, U640Repr>::$func(&$argument)
//             },
//             11 => {
//                 $implementation::<$modulus_limbs, U704Repr>::$func(&$argument)
//             },
//             12 => {
//                 $implementation::<$modulus_limbs, U768Repr>::$func(&$argument)
//             },
//             13 => {
//                 $implementation::<$modulus_limbs, U832Repr>::$func(&$argument)
//             },
//             14 => {
//                 $implementation::<$modulus_limbs, U896Repr>::$func(&$argument)
//             },
//             15 => {
//                 $implementation::<$modulus_limbs, U960Repr>::$func(&$argument)
//             },
//             16 => {
//                 $implementation::<$modulus_limbs, U1024Repr>::$func(&$argument)
//             },

//             group_limbs => {
//                 unimplemented!("unimplemented for {} group limbs", group_limbs);
//             }
//         }
//     }
// }