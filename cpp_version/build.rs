extern crate cc;

fn main() {
    cc::Build::new()
        .cpp(true) // Switch to C++ library compilation.
        .flag("-std=c++1z")
        .include("eip1962cpp/include")
        .file("eip1962cpp/src/api.cpp")
        .file("eip1962cpp/src/common.cpp")
        .file("eip1962cpp/src/wrapper.cpp")
        .file("eip1962cpp/src/repr.cpp")
        .warnings(false)
        // .static_flag(true)
        // .opt_level_str("3")
        .compile("eip1962cpp.a");
}