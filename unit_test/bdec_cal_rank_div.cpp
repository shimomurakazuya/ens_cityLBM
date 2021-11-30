
#include <string>
#include <iostream>
#include <functional>
#include <map>
#include <tuple>

#include "../util/divisor.hpp"
#include "../util/argstr.hpp"

int main(int argc, char** argv) {
    auto args = util::argstr(argc, argv);

    for(auto a: args) {
        auto ncpu = std::stoi(a);

        //// a brute-force map for BDEC (2021/07/16)
        const std::map<int, std::tuple<int, int>> ncpu_xy_map {
            // ss 40 -- 360
            {40, {10, 4}},
            {80, {10, 8}},
            {120, {12, 10}},
            {240, {20, 12}},
            {360, {20, 18}},
            // ss 36 -- 324
            {36, {6, 6}},
            {64, {8, 8}},
            {144, {12, 12}},
            {256, {16, 16}},
            {324, {18, 18}}
        };
        const auto [ncpu_x, ncpu_y] = ncpu_xy_map.at(ncpu);

        std::cout << " ncpux,ncpuy(ncpu) = " << ncpu_x << ", " << ncpu_y << " (" << ncpu << ")" << std::endl;
    }
}
