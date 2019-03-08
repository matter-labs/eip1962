#![recursion_limit = "1024"]

extern crate byteorder;

extern crate proc_macro;
extern crate proc_macro2;
extern crate syn;
#[macro_use]
extern crate quote;

extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;

// use num_bigint::BigUint;
// use num_integer::Integer;
// use num_traits::{One, ToPrimitive, Zero};
// use quote::TokenStreamExt;
use std::str::FromStr;

#[proc_macro_derive(ElementRepresentation, attributes(NumberOfLimbs))]
pub fn element_repr(input: proc_macro::TokenStream) -> proc_macro::TokenStream {
    // Parse the type definition
    let ast: syn::DeriveInput = syn::parse(input).unwrap();

    // The struct we're deriving for is a wrapper around a "Repr" type we must construct.
    let repr_ident = fetch_wrapped_ident(&ast.data)
        .expect("PrimeField derive only operates over tuple structs of a single item");

    // // The arithmetic in this library only works if the modulus*2 is smaller than the backing
    // // representation. Compute the number of limbs we need.
    // let mut limbs = 1;
    // {
    //     let mod2 = (&modulus) << 1; // modulus * 2
    //     let mut cur = BigUint::one() << 64; // always 64-bit limbs for now
    //     while cur < mod2 {
    //         limbs += 1;
    //         cur = cur << 64;
    //     }
    // }

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
        if let Some(meta) = attr.interpret_meta() {
            match meta {
                syn::Meta::NameValue(nv) => {
                    if nv.ident.to_string() == name {
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
    quote! {

        #[derive(Copy, Clone, PartialEq, Eq, Default)]
        pub struct #repr(
            pub [u64; #limbs]
        );

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

        impl crate::representation::ElementRepr for #repr {
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
        }
    }
}