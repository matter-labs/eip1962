macro_rules! decl_fp {
    ($repr_type:ty) => ( Fp::<'static, $repr_type, PrimeField<$repr_type>> );
}

macro_rules! repr_into_fp {
    ($repr:expr, $repr_type:ty, $field:ident) => ({
        Fp::<'static, $repr_type, PrimeField<$repr_type>> {
            field: &$field,
            repr: $repr
        }
    });
}

macro_rules! decl_fp2 {
    ($repr_type:ty) => ( Fp2::<'static, $repr_type, PrimeField<$repr_type>> );
}

macro_rules! decl_fp6 {
    ($repr_type:ty) => ( Fp6::<'static, $repr_type, PrimeField<$repr_type>> );
}

macro_rules! decl_fp12 {
    ($repr_type:ty) => ( Fp12::<'static, $repr_type, PrimeField<$repr_type>> );
}

macro_rules! repr_into_fp2 {
    ($repr_c0:expr, $repr_c1:expr, $repr_type:ty, $field:ident) => ({
        Fp2::<'static, $repr_type, PrimeField<$repr_type>> {
            c0: $repr_c0,
            c1: $repr_c1,
            extension_field: &$field
        }
    });
}

macro_rules! repr_into_fp6 {
    ($repr_c0:expr, $repr_c1:expr, $repr_c2:expr, $repr_type:ty, $field:ident) => ({
        Fp6::<'static, $repr_type, PrimeField<$repr_type>> {
            c0: $repr_c0,
            c1: $repr_c1,
            c2: $repr_c2,
            extension_field: &$field
        }
    });
}

macro_rules! repr_into_fp12 {
    ($repr_c0:expr, $repr_c1:expr, $repr_type:ty, $field:ident) => ({
        Fp12::<'static, $repr_type, PrimeField<$repr_type>> {
            c0: $repr_c0,
            c1: $repr_c1,
            extension_field: &$field
        }
    });
}