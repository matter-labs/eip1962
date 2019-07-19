#include "api.h"
#include "wrapper.h"

int run(const char *i, uint32_t i_len, char *o, uint32_t *o_len, char *err, uint32_t *char_len) {
    std::vector<std::uint8_t> input;
    input.resize(i_len);
    std::copy(i, i + i_len, input.begin());
    auto result = run(input);
    if (auto answer = std::get_if<0>(&result))
    {
        std::copy(answer->begin(), answer->end(), o);
        *o_len = answer->size();
        return true;
    } else if (auto error_descr = std::get_if<1>(&result)) {
        auto str_len = error_descr->size();
        auto c_str = error_descr->c_str();
        std::copy(c_str, c_str + str_len + 1, err);
        *char_len = error_descr->size();
        return false;
    }

    return false;
}

int meter_gas(const char *i, uint32_t i_len, uint64_t *gas) {
    *gas = UINT64_MAX;
    return true;
}