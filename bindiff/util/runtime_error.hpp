// exception.hpp
// wrapper for std::runtime_error

#include <stdexcept>
#include <string>

#define STD_RUNTIME_ERROR(error_str) \
    std::runtime_error(std::string(error_str) + " at " + std::string(__FILE__) + ":" + std::to_string(__LINE__) + ": " + std::string(__PRETTY_FUNCTION__))

#define RUNTIME_ASSERT(comment) \
    if(!(comment)) { throw STD_RUNTIME_ERROR(#comment); }
