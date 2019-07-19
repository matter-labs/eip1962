#include "common.h"
#include <cstdio>
#include <cstdarg>
#include <alloca.h>

#include <string>
#include <iostream>

std::string err_concat(std::string const &a, std::string const &b)
{
    return a + b;
}

std::string stringf(const char *format, ...)
{
    va_list arg_list;
    va_start(arg_list, format);

    // SUSv2 version doesn't work for buf NULL/size 0, so try printing
    // into a small buffer that avoids the double-rendering and alloca path too...
    char short_buf[256];
    const size_t needed = vsnprintf(short_buf, sizeof short_buf,
                                    format, arg_list) +
                          1;
    if (needed <= sizeof short_buf)
        return short_buf;

    // need more space...

    // OPTION 1
    std::string result(needed, ' ');
    vsnprintf(result.data(), needed, format, arg_list);
    return result; // RVO ensures this is cheap
}
