#pragma once
#ifndef WORKERTHREAD_H_
#define WORKERTHREAD_H_


#include <iostream>
#include <thread>


class  WorkerThread {
private:
    std::thread* thread_;
    bool  is_allocated_ = false;

    std::string  wname_ = "none";

public:
     WorkerThread () {}
    ~WorkerThread () {}

public:
    void check_working_func() {
        std::cout << "working func : " << wname_ << std::endl;
    }

    void work(const bool  is_work, std::function<void()>  func)
    {
        if ( is_work ) {
            if ( is_allocated_ ) {
                std::cout << __PRETTY_FUNCTION__ << " : "
                          << "error worker thread is already allocated!!" << std::endl;
                exit(-1);
            }

            wname_ = __PRETTY_FUNCTION__;

            is_allocated_ = true;
            thread_ = new std::thread( func );
        }
    }

    void join(const bool  is_join)
    {
        if ( is_join ) {
            if ( is_allocated_ ) {
                thread_->join();

                delete thread_;
                is_allocated_ = false;
                wname_ = "none";
            }
        }
    }

};


#endif
