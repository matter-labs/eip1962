#ifndef H_API
#define H_API

#include <variant>
#include <vector>
#include <string>

// Main API function for ABI.
std::variant<std::vector<std::uint8_t>, std::basic_string<char>> run(std::vector<std::uint8_t> const &input);

#endif
