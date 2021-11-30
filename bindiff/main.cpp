#define VERBOSE 0
#define EPSILON 1e-30
#define F_HALO 999999.9

#include <iostream>
#include <fstream>
#include <cstdint> // std::intptr_t
#include <cmath> // nan
#include <algorithm> // max,min,abs
#include <zlib.h>
#include <iomanip>

#include "util/argstr.hpp"

using real = float;

int main(int argc, char** argv) {

    // args: dir_left, dir_right, file1, file2, file3, ...
    const auto& args = util::argstr(argc, argv);
    const auto& prefix0 = args.at(0);
    const auto& prefix1 = args.at(1);
    std::cout << "comparing directories: " << prefix0 << " -- " << prefix1 << std::endl;

    real diffmin = 99999.9;
    real diffmax = 0.0;
    long double e2sum = 0.0;
    long double v2sum = 0.0;

    // for each file
    #pragma omp parallel for reduction(+:e2sum,v2sum) reduction(min:diffmin) reduction(max:diffmax)
    for(int i = 2; i < args.size(); i++) {
        const auto filename = args.at(i);
        /// do bin_diff of file
        auto file0 = gzopen((prefix0 + "/" + filename).c_str() , "rb");
        auto file1 = gzopen((prefix1 + "/" + filename).c_str() , "rb");
        constexpr auto bufsize = 4*4*4*8;
        auto buf0 = std::array<real, bufsize>();
        auto buf1 = std::array<real, bufsize>();
        auto diffmin_f = NAN;
        auto diffmax_f = NAN;
        auto count = 0l;

        while( bufsize*sizeof(real) ==  gzread(file0, buf0.data(), bufsize*sizeof(real))
            && bufsize*sizeof(real) ==  gzread(file1, buf1.data(), bufsize*sizeof(real))
        ) {
            for(auto j=0; j<bufsize; j++) {
                if(buf0.at(j) >= F_HALO/2 or buf1.at(j) >= F_HALO/2) { continue; }
                const auto diff = std::abs(buf0.at(j) - buf1.at(j));
                if(std::isnan(diff) || std::isnan(diffmin_f)) {
                    diffmin_f = diffmax_f = diff; 
                } else {
                    diffmin_f = std::min(diffmin_f, diff);
                    diffmax_f = std::max(diffmax_f, diff);
                }
                v2sum += buf0.at(j) * buf0.at(j);
                e2sum += diff*diff;
                count ++;
            }
        }

        gzclose(file0);
        gzclose(file1);

        if(VERBOSE >= 4) {
            std::cout << filename << ": " << diffmin_f << " --- " << diffmax_f
                << " with val_count = " << count << std::endl;
        } else {
            if(diffmax_f < EPSILON) {
                std::cout << '.' << std::flush;
            } else if(count == 0) {
                std::cout << '-' << std::flush;
            } else {
                std::cout << filename[0] << std::flush;
                if(VERBOSE >= 2) {
                    std::cout << std::endl;
                    std::cout << filename << ": " << diffmin_f << " --- " << diffmax_f
                        << " with val_count = " << count << std::endl;
                }
            }
        }
        if(!std::isnan(diffmin_f) && !std::isnan(diffmax_f)) {
            diffmin = std::min(diffmin, diffmin_f);
            diffmax = std::max(diffmax, diffmax_f);
        }
    }

    std::cout << std::endl;

    if(VERBOSE >= 1) {
        std::cout << "final result: " <<  diffmin << " --- " << diffmax << std::endl;
        std::cout << " with dirs: " << prefix0 << " -- " << prefix1 << std::endl;
    }
    if(diffmax >= EPSILON) {
        std::cerr << "result is not same; " << std::endl
           << "  max_diff = " << diffmax << ", " << std::endl
           << "  l2error = " << std::sqrt(e2sum/v2sum) << std::endl;
        std::cout << std::endl;
        return 2;
    }
} // main()

