#pragma once
#ifndef TIME_H_
#define TIME_H_

#include <iostream>
#include <sys/time.h>
#include <string>
#include <functional>


class Timer {
private:
    struct timeval  t_begin_;
    struct timeval  t_end_;
    float  total_elapsed_time_;

    std::string message_;

public :
    Timer() { Timer("NO MESSAGE"); }

    Timer(const std::string& message) :
        message_(message)
    {
        gettimeofday(&t_begin_, NULL);
        gettimeofday(&t_end_,   NULL);

        total_elapsed_time_ = 0.0;
    }

    ~Timer() {}

public:
    float  measure_func(const std::function<void()> func);
    float  measure_func(const std::function<void()> func, const std::string  message);

    void start();
    void stop ();

    float clean();

    float get_total_elapsed_time_msec() const;
    float get_elapsed_time_msec() const;

    void  print_total_elapsed_time_msec() const;
    void  print_elapsed_time_msec() const;


    float  cal_elapsed_time_msec (
        const struct timeval&   begin,
        const struct timeval&   end
        )
    const ;

};


#endif
