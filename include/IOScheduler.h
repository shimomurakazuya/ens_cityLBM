#pragma once
#ifndef IOSCHEDULER_H_
#define IOSCHEDULER_H_


#include <iostream>
#include <thread>


class  ThreadScheduler {
private:
    int   nstep_ = 1;

public:
    ThreadScheduler () {}
    ~ThreadScheduler () {}

public:
    bool  is_start(const int  i)  const { return  ( i%nstep_ == 0        ) ? true : false; }
    bool  is_end  (const int  i)  const { return  ( i%nstep_ == nstep_-1 ) ? true : false; }

public:
    void  set_nstep(const int  nstep) { nstep_ = nstep; }


};


#endif
