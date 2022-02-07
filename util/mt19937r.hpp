#ifndef MT19937R_HPP_
#define MT19937R_HPP_
/// restartable mt19937 pseudo random engine
/// needs C++ >= 11

#include <random>
#include <string>
#include <cstdio>
#include <iostream>
#include "runtime_error.hpp"

namespace util {

class mt19937r {
public:
    using engine_type = std::mt19937;
    using result_type = engine_type::result_type;
    static constexpr auto default_seed = engine_type::default_seed;
private:
    std::mt19937 engine_;
    result_type seed_;
    unsigned long long count_;

public:
    explicit mt19937r(result_type value = engine_type::default_seed)
    : engine_(value),
      seed_(value),
      count_(0)
    {}

    mt19937r(const mt19937r&) = default;
    mt19937r(mt19937r&&) = default;

    void seed(result_type value = default_seed) { 
        engine_.seed(value);
        seed_ = value; 
        count_ = 0;
    }

    result_type operator()() { return count_++, engine_(); }
    void discard(unsigned long long z) { for(unsigned long long i=0; i<z; i++) { (*this)(); } }
    
    static result_type min() { return engine_type::min(); }
    static result_type max() { return engine_type::max(); }

public:
    void save(std::string filename, int verbose=0) const {
        FILE* fp = fopen(filename.c_str(), "wb");
        runtime_assert(fp != NULL, "IOError");
        fwrite(&seed_, sizeof(result_type), 1, fp);
        fwrite(&count_, sizeof(unsigned long long), 1, fp);
        fclose(fp);
        if(verbose >= 1) {
            std::cout << __PRETTY_FUNCTION__ << ": " << filename << ", seed = " << seed_ << ", count = " << count_ << std::endl;
        }
    }

    void restart(result_type seed_value, unsigned long long count) {
        seed(seed_value);
        discard(count);
    }

    void load(std::string filename, int verbose=0) {
        FILE* fp = fopen(filename.c_str(), "rb");
        runtime_assert(fp != NULL, "IOError");
        fread(&seed_, sizeof(result_type), 1, fp);
        fread(&count_, sizeof(unsigned long long), 1, fp);
        fclose(fp);
        restart(seed_, count_);
        if(verbose >= 1) {
            std::cout << __PRETTY_FUNCTION__ << ": " << filename << ", seed = " << seed_ << ", count = " << count_ << std::endl;
        }
    }
};

} // namespace

#endif
