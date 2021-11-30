// exception.hpp
// wrapper for std::runtime_error

#include <stdexcept>
#include <string>

#define STD_RUNTIME_ERROR(error_str) \
    std::runtime_error(std::string(__FILE__) + ":" + std::to_string(__LINE__) + ": " + std::string(__PRETTY_FUNCTION__) + ": assertion failed: " + std::string(error_str))

#define runtime_assert(con, comment) \
    if(!(con)) { throw STD_RUNTIME_ERROR(std::string(comment)); }
