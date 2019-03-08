use byteorder;
use std::fmt;
use std::error::Error;
use std::io::{self, Read, Write};

use crate::traits::FieldElement;

/// This trait represents a wrapper around a biginteger which can encode any element of a particular
/// prime field. It is a smart wrapper around a sequence of `u64` limbs, least-significant digit
/// first.
pub trait ElementRepr:
    Sized
    + Copy
    + Clone
    + Eq
    + Ord
    + Send
    + Sync
    + Default
    + fmt::Debug
    + fmt::Display
    + 'static
    + AsRef<[u64]>
    + AsMut<[u64]>
    + From<u64>
{
    /// Subtract another represetation from this one.
    fn sub_noborrow(&mut self, other: &Self);

    /// Add another representation to this one.
    fn add_nocarry(&mut self, other: &Self);

    /// Compute the number of bits needed to encode this number. Always a
    /// multiple of 64.
    fn num_bits(&self) -> u32;

    /// Returns true iff this number is zero.
    fn is_zero(&self) -> bool;

    /// Returns true iff this number is odd.
    fn is_odd(&self) -> bool;

    /// Returns true iff this number is even.
    fn is_even(&self) -> bool;

    /// Performs a rightwise bitshift of this number, effectively dividing
    /// it by 2.
    fn div2(&mut self);

    /// Performs a rightwise bitshift of this number by some amount.
    fn shr(&mut self, amt: u32);

    /// Performs a leftwise bitshift of this number, effectively multiplying
    /// it by 2. Overflow is ignored.
    fn mul2(&mut self);

    /// Performs a leftwise bitshift of this number by some amount.
    fn shl(&mut self, amt: u32);

    /// Writes this `PrimeFieldRepr` as a big endian integer.
    fn write_be<W: Write>(&self, mut writer: W) -> io::Result<()> {
        use byteorder::{BigEndian, WriteBytesExt};

        for digit in self.as_ref().iter().rev() {
            writer.write_u64::<BigEndian>(*digit)?;
        }

        Ok(())
    }

    /// Reads a big endian integer into this representation.
    fn read_be<R: Read>(&mut self, mut reader: R) -> io::Result<()> {
        use byteorder::{BigEndian, ReadBytesExt};

        for digit in self.as_mut().iter_mut().rev() {
            *digit = reader.read_u64::<BigEndian>()?;
        }

        Ok(())
    }

    /// Writes this `PrimeFieldRepr` as a little endian integer.
    fn write_le<W: Write>(&self, mut writer: W) -> io::Result<()> {
        use byteorder::{LittleEndian, WriteBytesExt};

        for digit in self.as_ref().iter() {
            writer.write_u64::<LittleEndian>(*digit)?;
        }

        Ok(())
    }

    /// Reads a little endian integer into this representation.
    fn read_le<R: Read>(&mut self, mut reader: R) -> io::Result<()> {
        use byteorder::{LittleEndian, ReadBytesExt};

        for digit in self.as_mut().iter_mut() {
            *digit = reader.read_u64::<LittleEndian>()?;
        }

        Ok(())
    }
}

// /// This represents an element of a prime field.
// pub trait PrimeFieldElement: FieldElement {
//     /// The prime field can be converted back and forth into this biginteger
//     /// representation.
//     type Repr: ElementRepr + From<Self>;

//     /// Convert this prime field element into a biginteger representation.
//     fn from_repr(repr: Self::Repr) -> Result<Self, PrimeFieldDecodingError>;

//     /// Convert a biginteger representation into a prime field element, if
//     /// the number is an element of the field.
//     fn into_repr(&self) -> Self::Repr;
// }

/// An error that may occur when trying to interpret a `PrimeFieldRepr` as a
/// `PrimeField` element.
#[derive(Debug)]
pub enum RepresentationDecodingError {
    /// The encoded value is not in the field
    NotInField(String),
}

impl Error for RepresentationDecodingError {
    fn description(&self) -> &str {
        match *self {
            RepresentationDecodingError::NotInField(..) => "not an element of the field",
        }
    }
}

impl fmt::Display for RepresentationDecodingError {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        match *self {
            RepresentationDecodingError::NotInField(ref repr) => {
                write!(f, "{} is not an element of the field", repr)
            }
        }
    }
}