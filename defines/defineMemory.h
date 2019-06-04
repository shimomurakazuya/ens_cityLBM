#pragma once
#ifndef DEFINEMEMORY_H_
#define DEFINEMEMORY_H_


#include <iostream>
#include <cstdlib>
#include "defineCal.h"


enum struct MemType {
    Host,
    Device,
    Pinned,
    Default,
    Managed,
    NullPtr
};


#endif
