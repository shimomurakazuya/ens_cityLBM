#pragma once
#ifndef DEFINEOBJBC_H_
#define DEFINEOBJBC_H_


#include <iostream>
#include <cstdlib>

namespace BCTypes {
    const int CalRegion   = 0;
    const int BCRegion    = 1;

    const int BCNeumann   = 2;
    const int BCDirichlet = 3;
};


#endif
