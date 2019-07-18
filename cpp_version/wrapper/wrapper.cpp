#include "api.h"
#include "wrapper.h"

int run(const char *i, int i_len, char *o) {
    std::vector<std::uint8_t> input;
    input.resize(i_len);
    std::copy(i, i + i_len, input.begin());
    auto result = run(input);
    try {
        auto output = std::get<std::vector<std::uint8_t>>(result);
        std::copy(output.begin(), output.end(), o);
        return output.size();
    }
    catch (const std::bad_variant_access&) {
        return 0;
    }
}