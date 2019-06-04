#include "Timer.h"


float  Timer::measure_func(const std::function<void()> func)
{
    start();
    func();
    stop();

    const float  elapsed_time = get_elapsed_time_msec();
    total_elapsed_time_ += elapsed_time;

    return  elapsed_time;
}


float  Timer::measure_func(const std::function<void()> func, const std::string  message)
{
    message_ = message;
    return  measure_func( func );
}


void  Timer::start()
{
    gettimeofday(&t_begin_, NULL);
    gettimeofday(&t_end_,   NULL);
}


void  Timer::stop()
{
    gettimeofday(&t_end_, NULL);
}


float  Timer::clean()
{
    total_elapsed_time_ = 0.0;
    return  total_elapsed_time_;
}


float  Timer::get_total_elapsed_time_msec ()
const
{
    return  total_elapsed_time_;
}


float  Timer::get_elapsed_time_msec ()
const
{
    return  cal_elapsed_time_msec(t_begin_, t_end_);
}


void  Timer::print_total_elapsed_time_msec()
const
{
    printf("%s = %f sec.\n", message_.c_str(), get_total_elapsed_time_msec() / 1000.0);
}


void  Timer::print_elapsed_time_msec()
const
{
    printf("%s = %f sec.\n", message_.c_str(), get_elapsed_time_msec() / 1000.0);
}


float  Timer::cal_elapsed_time_msec (
    const struct timeval&   begin,
    const struct timeval&   end
    )
const
{
    // msec //
    return    (end.tv_sec  - begin.tv_sec) * 1000
            + (end.tv_usec - begin.tv_usec) / 1000.0;
}
