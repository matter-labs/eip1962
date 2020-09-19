#![cfg_attr(not(feature = "std"), no_std)]
#![recursion_limit = "1024"]

extern crate byteorder;

extern crate proc_macro;
extern crate proc_macro2;
extern crate syn;

#[cfg(not(feature = "std"))]
#[macro_use]
extern crate alloc;

#[cfg(not(feature = "std"))]
mod std {
    pub use sp_std::*;
}

use crate::std::alloc::string::String;
use std::str::FromStr;

#[macro_use]
extern crate quote;

use quote::TokenStreamExt;

#[proc_macro_derive(ElementRepresentation, attributes(NumberOfLimbs))]
pub fn element_repr(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    // Parse the type definition
    let ast: syn::DeriveInput = syn::parse(input).unwrap();

    // The struct we're deriving for is a wrapper around a "Repr" type we must construct.
    let repr_ident = fetch_wrapped_ident(&ast.data)
        .expect("PrimeField derive only operates over tuple structs of a single item");

    let limbs_str = fetch_attr("NumberOfLimbs", &ast.attrs)
        .expect("Please supply a representation length in terms of 64 bit limbs");

    let limbs = usize::from_str(&limbs_str).expect("Number of limbs must be a number");

    let mut gen = proc_macro2::TokenStream::new();

    gen.extend(prime_field_repr_impl(&repr_ident, limbs));

    // Return the generated impl
    gen.into()
}

/// Fetches the ident being wrapped by the type we're deriving.
fn fetch_wrapped_ident(body: &syn::Data) -> Option<syn::Ident> {
    match body {
        &syn::Data::Struct(ref variant_data) => match variant_data.fields {
            syn::Fields::Unnamed(ref fields) => {
                if fields.unnamed.len() == 1 {
                    match fields.unnamed[0].ty {
                        syn::Type::Path(ref path) => {
                            if path.path.segments.len() == 1 {
                                return Some(path.path.segments[0].ident.clone());
                            }
                        }
                        _ => {}
                    }
                }
            }
            _ => {}
        },
        _ => {}
    };

    None
}

/// Fetch an attribute string from the derived struct.
fn fetch_attr(name: &str, attrs: &[syn::Attribute]) -> Option<String> {
    for attr in attrs {
        if let Ok(meta) = attr.parse_meta() {
            match meta {
                syn::Meta::NameValue(nv) => {
                    if nv.path.is_ident(name) {
                        match nv.lit {
                            syn::Lit::Str(ref s) => return Some(s.value()),
                            _ => {
                                panic!("attribute {} should be a string", name);
                            }
                        }
                    }
                }
                _ => {
                    panic!("attribute {} should be a string", name);
                }
            }
        }
    }

    None
}

// Implement PrimeFieldRepr for the wrapped ident `repr` with `limbs` limbs.
fn prime_field_repr_impl(repr: &syn::Ident, limbs: usize) -> proc_macro2::TokenStream {
    // Returns r{n} as an ident.
    fn get_temp(n: usize) -> syn::Ident {
        syn::Ident::new(&format!("r{}", n), proc_macro2::Span::call_site())
    }

    // The parameter list for the mont_reduce() internal method.
    // r0: u64, mut r1: u64, mut r2: u64, ...
    let mut mont_paramlist = proc_macro2::TokenStream::new();
    mont_paramlist.append_separated(
        (0..(limbs * 2)).map(|i| (i, get_temp(i))).map(|(i, x)| {
            if i != 0 {
                quote!{mut #x: u64}
            } else {
                quote!{#x: u64}
            }
        }),
        proc_macro2::Punct::new(',', proc_macro2::Spacing::Alone),
    );

    let mut into_normal_repr_params = proc_macro2::TokenStream::new();
    into_normal_repr_params.append_separated(
        (0..limbs)
            .map(|i| quote!{ self.0[#i] })
            .chain((0..limbs).map(|_| quote!{0})),
        proc_macro2::Punct::new(',', proc_macro2::Spacing::Alone),
    );

    // Implement montgomery reduction for some number of limbs
    fn mont_impl(limbs: usize) -> proc_macro2::TokenStream {
        let mut gen = proc_macro2::TokenStream::new();

        for i in 0..limbs {
            {
                let temp = get_temp(i);
                gen.extend(quote!{
                    let k = #temp.wrapping_mul(mont_inv);
                    let mut carry = 0;
                    crate::arithmetics::mac_with_carry(#temp, k, modulus.0[0], &mut carry);
                });
            }

            for j in 1..limbs {
                let temp = get_temp(i + j);
                gen.extend(quote!{
                    #temp = crate::arithmetics::mac_with_carry(#temp, k, modulus.0[#j], &mut carry);
                });
            }

            let temp = get_temp(i + limbs);

            if i == 0 {
                gen.extend(quote!{
                    #temp = crate::arithmetics::adc(#temp, 0, &mut carry);
                });
            } else {
                gen.extend(quote!{
                    #temp = crate::arithmetics::adc(#temp, carry2, &mut carry);
                });
            }

            if i != (limbs - 1) {
                gen.extend(quote!{
                    let carry2 = carry;
                });
            }
        }

        for i in 0..limbs {
            let temp = get_temp(limbs + i);

            gen.extend(quote!{
                self.0[#i] = #temp;
            });
        }

        gen
    }

    fn sqr_impl(a: proc_macro2::TokenStream, limbs: usize) -> proc_macro2::TokenStream {
        let mut gen = proc_macro2::TokenStream::new();

        for i in 0..(limbs - 1) {
            gen.extend(quote!{
                let mut carry = 0;
            });

            for j in (i + 1)..limbs {
                let temp = get_temp(i + j);
                if i == 0 {
                    gen.extend(quote!{
                        let #temp = crate::arithmetics::mac_with_carry(0, #a.0[#i], #a.0[#j], &mut carry);
                    });
                } else {
                    gen.extend(quote!{
                        let #temp = crate::arithmetics::mac_with_carry(#temp, #a.0[#i], #a.0[#j], &mut carry);
                    });
                }
            }

            let temp = get_temp(i + limbs);

            gen.extend(quote!{
                let #temp = carry;
            });
        }

        for i in 1..(limbs * 2) {
            let temp0 = get_temp(limbs * 2 - i);
            let temp1 = get_temp(limbs * 2 - i - 1);

            if i == 1 {
                gen.extend(quote!{
                    let #temp0 = #temp1 >> 63;
                });
            } else if i == (limbs * 2 - 1) {
                gen.extend(quote!{
                    let #temp0 = #temp0 << 1;
                });
            } else {
                gen.extend(quote!{
                    let #temp0 = (#temp0 << 1) | (#temp1 >> 63);
                });
            }
        }

        gen.extend(quote!{
            let mut carry = 0;
        });

        for i in 0..limbs {
            let temp0 = get_temp(i * 2);
            let temp1 = get_temp(i * 2 + 1);
            if i == 0 {
                gen.extend(quote!{
                    let #temp0 = crate::arithmetics::mac_with_carry(0, #a.0[#i], #a.0[#i], &mut carry);
                });
            } else {
                gen.extend(quote!{
                    let #temp0 = crate::arithmetics::mac_with_carry(#temp0, #a.0[#i], #a.0[#i], &mut carry);
                });
            }

            gen.extend(quote!{
                let #temp1 = crate::arithmetics::adc(#temp1, 0, &mut carry);
            });
        }

        let mut mont_calling = proc_macro2::TokenStream::new();
        mont_calling.append_separated(
            (0..(limbs * 2)).map(|i| get_temp(i)),
            proc_macro2::Punct::new(',', proc_macro2::Spacing::Alone),
        );

        gen.extend(quote!{
            self.mont_partial_reduce(modulus, mont_inv, #mont_calling);
        });

        gen
    }

    fn mul_impl(
        a: proc_macro2::TokenStream,
        b: proc_macro2::TokenStream,
        limbs: usize,
    ) -> proc_macro2::TokenStream {
        let mut gen = proc_macro2::TokenStream::new();

        for i in 0..limbs {
            gen.extend(quote!{
                let mut carry = 0;
            });

            for j in 0..limbs {
                let temp = get_temp(i + j);

                if i == 0 {
                    gen.extend(quote!{
                        let #temp = crate::arithmetics::mac_with_carry(0, #a.0[#i], #b.0[#j], &mut carry);
                    });
                } else {
                    gen.extend(quote!{
                        let #temp = crate::arithmetics::mac_with_carry(#temp, #a.0[#i], #b.0[#j], &mut carry);
                    });
                }
            }

            let temp = get_temp(i + limbs);

            gen.extend(quote!{
                let #temp = carry;
            });
        }

        let mut mont_calling = proc_macro2::TokenStream::new();
        mont_calling.append_separated(
            (0..(limbs * 2)).map(|i| get_temp(i)),
            proc_macro2::Punct::new(',', proc_macro2::Spacing::Alone),
        );

        gen.extend(quote!{
            self.mont_partial_reduce(modulus, mont_inv, #mont_calling);
        });

        gen
    }

    let squaring_impl = sqr_impl(quote!{self}, limbs);
    let multiply_impl = mul_impl(quote!{self}, quote!{other}, limbs);
    let montgomery_impl = mont_impl(limbs);

    quote! {

        #[derive(Copy, Clone, PartialEq, Eq, Default)]
        pub struct #repr(
            pub [u64; #limbs]
        );

        impl #repr {
            // #[inline(always)]
            // fn mont_reduce(
            //     &mut self,
            //     modulus: &#repr,
            //     mont_inv: u64,
            //     #mont_paramlist
            // )
            // {
            //     // The Montgomery reduction here is based on Algorithm 14.32 in
            //     // Handbook of Applied Cryptography
            //     // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

            //     #montgomery_impl
            //     self.reduce(modulus);
            // }

            #[inline(always)]
            fn mont_partial_reduce(
                &mut self,
                modulus: &#repr,
                mont_inv: u64,
                #mont_paramlist
            )
            {
                // The Montgomery reduction here is based on Algorithm 14.32 in
                // Handbook of Applied Cryptography
                // <http://cacr.uwaterloo.ca/hac/about/chap14.pdf>.

                #montgomery_impl
            }
        }

        impl ::std::fmt::Debug for #repr
        {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "0x")?;
                for i in self.0.iter().rev() {
                    write!(f, "{:016x}", *i)?;
                }

                Ok(())
            }
        }

        impl ::std::fmt::Display for #repr {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                write!(f, "0x")?;
                for i in self.0.iter().rev() {
                    write!(f, "{:016x}", *i)?;
                }

                Ok(())
            }
        }

        impl AsRef<[u64]> for #repr {
            #[inline(always)]
            fn as_ref(&self) -> &[u64] {
                &self.0
            }
        }

        impl AsMut<[u64]> for #repr {
            #[inline(always)]
            fn as_mut(&mut self) -> &mut [u64] {
                &mut self.0
            }
        }

        impl From<u64> for #repr {
            #[inline(always)]
            fn from(val: u64) -> #repr {
                use std::default::Default;

                let mut repr = Self::default();
                repr.0[0] = val;
                repr
            }
        }

        impl Ord for #repr {
            #[inline(always)]
            fn cmp(&self, other: &#repr) -> std::cmp::Ordering {
                for (a, b) in self.0.iter().rev().zip(other.0.iter().rev()) {
                    if a < b {
                        return std::cmp::Ordering::Less
                    } else if a > b {
                        return std::cmp::Ordering::Greater
                    }
                }

                std::cmp::Ordering::Equal
            }
        }

        impl PartialOrd for #repr {
            #[inline(always)]
            fn partial_cmp(&self, other: &#repr) -> Option<std::cmp::Ordering> {
                Some(self.cmp(other))
            }
        }

        impl crate::representation::IntoWnaf for #repr {
            #[inline]
            fn wnaf(&self, window: u32) -> Vec<i64> {
                let mut res = vec![];

                let mut e = *self;
                let max = (1 << window) as i64;
                let midpoint = (1 << (window-1)) as i64;
                let modulus_mask = ((1 << window) - 1) as u64;
                while !e.is_zero() {
                    let mut z: i64 = 0;
                    if e.is_odd() {
                        let masked_bits = (e.0[0] & modulus_mask) as i64;
                        // z = midpoint - (e.0[0] & modulus_mask) as i64;
                        if masked_bits > midpoint {
                            z = masked_bits - max;
                            e.add_nocarry(&Self::from((-z) as u64));
                        } else {
                            z = masked_bits;
                            e.sub_noborrow(&Self::from(z as u64));
                        }
                    } else {
                        z = 0;
                    }
                    res.push(z);
                    e.div2();
                }

                res
            }
        }

        impl crate::representation::ElementRepr for #repr {
            const NUM_LIMBS: usize = #limbs;

            #[inline(always)]
            fn is_odd(&self) -> bool {
                self.0[0] & 1 == 1
            }

            #[inline(always)]
            fn is_even(&self) -> bool {
                !self.is_odd()
            }

            #[inline(always)]
            fn is_zero(&self) -> bool {
                self.0.iter().all(|&e| e == 0)
            }

            #[inline(always)]
            fn shr(&mut self, mut n: u32) {
                if n as usize >= 64 * #limbs {
                    *self = Self::from(0);
                    return;
                }

                while n >= 64 {
                    let mut t = 0;
                    for i in self.0.iter_mut().rev() {
                        std::mem::swap(&mut t, i);
                    }
                    n -= 64;
                }

                if n > 0 {
                    let mut t = 0;
                    for i in self.0.iter_mut().rev() {
                        let t2 = *i << (64 - n);
                        *i >>= n;
                        *i |= t;
                        t = t2;
                    }
                }
            }

            #[inline(always)]
            fn div2(&mut self) {
                let mut t = 0;
                for i in self.0.iter_mut().rev() {
                    let t2 = *i << 63;
                    *i >>= 1;
                    *i |= t;
                    t = t2;
                }
            }

            #[inline(always)]
            fn mul2(&mut self) {
                let mut last = 0;
                for i in &mut self.0 {
                    let tmp = *i >> 63;
                    *i <<= 1;
                    *i |= last;
                    last = tmp;
                }
            }

            #[inline(always)]
            fn shl(&mut self, mut n: u32) {
                if n as usize >= 64 * #limbs {
                    *self = Self::from(0);
                    return;
                }

                while n >= 64 {
                    let mut t = 0;
                    for i in &mut self.0 {
                        std::mem::swap(&mut t, i);
                    }
                    n -= 64;
                }

                if n > 0 {
                    let mut t = 0;
                    for i in &mut self.0 {
                        let t2 = *i >> (64 - n);
                        *i <<= n;
                        *i |= t;
                        t = t2;
                    }
                }
            }

            #[inline(always)]
            fn num_bits(&self) -> u32 {
                let mut ret = (#limbs as u32) * 64;
                for i in self.0.iter().rev() {
                    let leading = i.leading_zeros();
                    ret -= leading;
                    if leading != 64 {
                        break;
                    }
                }

                ret
            }

            #[inline(always)]
            fn add_nocarry(&mut self, other: &#repr) {
                let mut carry = 0;

                for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
                    *a = crate::arithmetics::adc(*a, *b, &mut carry);
                }
            }

            #[inline(always)]
            fn sub_noborrow(&mut self, other: &#repr) {
                let mut borrow = 0;

                for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
                    *a = crate::arithmetics::sbb(*a, *b, &mut borrow);
                }
            }

            #[inline]
            fn mont_mul_assign(&mut self, other: &#repr, modulus: &#repr, mont_inv: u64)
            {
                #multiply_impl
                self.reduce(modulus);
            }

            #[inline]
            fn mont_square(&mut self, modulus: &#repr, mont_inv: u64)
            {
                #squaring_impl
                self.reduce(modulus);
            }

            #[inline]
            fn mont_mul_assign_with_partial_reduction(&mut self, other: &#repr, modulus: &#repr, mont_inv: u64)
            {
                #multiply_impl
            }

            #[inline]
            fn mont_square_with_partial_reduction(&mut self, modulus: &#repr, mont_inv: u64)
            {
                #squaring_impl
            }

            #[inline(always)]
            fn into_normal_repr(&self, modulus: &#repr, mont_inv: u64) -> #repr {
                let mut r = *self;
                r.mont_partial_reduce(
                    modulus,
                    mont_inv,
                    #into_normal_repr_params
                );

                r.reduce(modulus);

                r
            }

            #[inline(always)]
            fn reduce(
                &mut self,
                modulus: &#repr
            )
            {
                if &*self >= modulus {
                    self.sub_noborrow(&modulus);
                }
            }
        }
    }
}