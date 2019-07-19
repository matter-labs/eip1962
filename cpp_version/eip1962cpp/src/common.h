#ifndef H_COMMON
#define H_COMMON

// Comment this define to turn on asserts and println
#define NDEBUG

#include <assert.h>
#include <optional>
#include <vector>
#include <cstdint>
#include <string>
#include <tuple>
#include <string>
#include <iostream>

// **************** RUST equivalent TYPES for ease of porting, and is also shorter to write ******************** //
using namespace std;

// C++ version of Rust u8
typedef std::uint8_t u8;

// C++ version of Rust u16
typedef std::uint16_t u16;

// C++ version of Rust u32
typedef std::uint32_t u32;

// C++ version of Rust i32
typedef std::int32_t i32;

// C++ version of Rust u64
typedef std::uint64_t u64;

// C++ version of Rust i64
typedef std::int64_t i64;

// C++ version of Rust usize
typedef std::size_t usize;

// C++ shortcut for Rust Option
template <class T>
using Option = std::optional<T>;

typedef const std::string str;

// **************** ERRORS as DEFINE so that the compiler can see that an exception is being thrown and stop with the warnings. ************************ //

#define unimplemented(ERR) throw std::runtime_error(err_concat("Unimplemented, ", ERR));

#define unreachable(ERR) throw std::runtime_error(err_concat("Unreachable reached", ERR));

#define api_err(ERR) throw std::domain_error(ERR);

#define input_err(ERR) throw std::domain_error(err_concat("Invalid input parameters, ", ERR));

#define unknown_parameter_err(ERR) throw std::domain_error(err_concat("parameter has value out of bounds, ", ERR));

#define unexpected_zero_err(ERR) throw std::domain_error(err_concat("parameter expected to be non-zero, ", ERR));

// **************** OTHER DEFINES ******************//

#define UNUSED(x) (void)(x)

// **************** FUNCTIONS ********************* //

std::string err_concat(std::string const &a, std::string const &b);

std::string stringf(const char *format, ...);

template <class... Args>
void println(const char *format, Args &&... args)
{
#ifndef NDEBUG
    std::cout << stringf(format, std::forward<Args>(args)...) << std::endl;
#endif
}

// ************** THINGS without a good place to put, so they are here ******************** //

enum TwistType
{
    D,
    M
};

#endif