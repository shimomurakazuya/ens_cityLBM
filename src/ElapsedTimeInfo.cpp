#include "ElapsedTimeInfo.h"


void  ElapsedTimeInfo::
InitElapsedTimeInfo(const std::string& fname, const Parameters& parameters)
{
    is_ios_base_app_ = false;
    fname_ = fname;

    setOutputPeriod(parameters);
}


void  ElapsedTimeInfo::
ReadElapsedTimeInfo(const std::string& fname, const Parameters& parameters)
{
    is_ios_base_app_ = true;
    fname_ = fname;

    setOutputPeriod(parameters);
}


void  ElapsedTimeInfo::
OutputElapsedTimeInfo(int t)
{
    if ( !flag_measure_ ) { return; }
    if ( t%output_period_ != 0 ) { return; }


    std::ofstream  fout;
    auto&& fname = filename(comm_.world().rank());

    if ( is_ios_base_app_ ) {
        fout.open(fname, std::ios::app);
    }
    else                    {
        fout.open(fname);

        if (true /* rank_ == 0 */) {
            fout << "func_name" << "," << "rank" << ","
                 << "step" << ","
                 << "nn_leaf" << "," << "lv" << "," << "count" << ","
                 << "elapsed_time_msec" << ","
                 << "elapsed_time_msec_per_leaf"
                 << std::endl;
        }
    }

    for (auto elem : st_elapsedTimeInfo_) {
        int  count = (elem.count != 0) ? elem.count : 1;

        fout << elem.func_name << "," << elem.rank << ","
             << elem.step << ","
             << elem.nn_leaf << "," << elem.lv << "," << elem.count << ","
             << elem.elapsed_time_msec << ","
             << elem.elapsed_time_msec/count
             << std::endl;
    }

    fout.close();


    // flags&reset //
    is_ios_base_app_ = true;
    st_elapsedTimeInfo_.clear ();
    st_elapsedTimeInfo_.shrink_to_fit();
}


void  ElapsedTimeInfo::
StartTimer(const std::string& func_name)
{
    if ( !flag_measure_ ) { return; }
    check_namelist(func_name);

    map_timer_[func_name].start();
}


void  ElapsedTimeInfo::
StopTimer(const std::string& func_name)
{
    if ( !flag_measure_ ) { return; }

    map_timer_[func_name].stop();
}


void  ElapsedTimeInfo::
SubmitElapsedTimeInfo(const std::string& func_name, int step)
{
    const int lv = 0;
    const int count = 1;

    SubmitElapsedTimeInfo(func_name, lv, count, step);
}


void  ElapsedTimeInfo::
SubmitElapsedTimeInfo(const std::string& func_name, int lv, int count, int step)
{
    stElapsedTimeInfo  elem;

    // submit //
    elem.func_name = func_name;

    elem.rank      = comm_.world().rank();
    elem.step      = step;
    elem.nn_leaf   = DefAMR::NN_LEAF;
    elem.lv        = lv;
    elem.count     = count;

    elem.elapsed_time_msec = map_timer_[func_name].get_elapsed_time_msec();

    // emplace_back //
    st_elapsedTimeInfo_.emplace_back( elem );
}


void  ElapsedTimeInfo::
setOutputPeriod(const Parameters& parameters)
{
    output_period_ = parameters.coefCFROutput().cout_step;
}


void  ElapsedTimeInfo::
check_namelist(const std::string& func_name)
{
    if ( map_timer_.count(func_name.c_str()) == 0 ) {
        timer_.emplace_back();
        map_timer_[ func_name.c_str() ] = timer_.back();
    }
}


std::string  ElapsedTimeInfo::
filename(const int rank)
const
{
    return    Foldernames::io_folder + "/"
            + fname_ + "-"
            + Filenames::ElapsedTimeInfo_name0
            + "-rank" + std::to_string(rank)
            + ".dat";
}
