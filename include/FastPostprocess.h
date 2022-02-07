#ifndef FASTPOSTPROCESS_H_
#define FASTPOSTPROCESS_H_

#include "defineAMR.h"
#include "Tree.h"
#include "MeshValue.h"
#include "ValueTimeAverage.h"
#include "Parameters.h"
#include "defineFilenames.h"

#include "managed_vector.hpp"
#include "range.hpp"

#include <string>
#include <iostream>
#include <sys/stat.h>

class FastPostprocess {
private:
    const MPICommEnsemble comm_;
    static constexpr real epsilon = 0.01;

    const std::string in_filename_ = "input-";
    std::string out_prefix0_() const { 
        return Foldernames::output_folder 
            + "/" + std::to_string(comm_.ensemble_id());
    }
    std::string out_prefix_ens_() const { 
        return out_prefix0_()
           + "/" + std::to_string(comm_.row_id()); 
    }
    std::string out_filename_(int minute, int fileno) const {
        std::string ret = out_prefix_ens_() + "/output";
        if(minute >= 0) { ret += "-minute" + std::to_string(minute); }
        ret += "-station" + std::to_string(fileno);
        ret += ".csv";
        return ret;
    }

public:
    struct st_ids { ssize_t i, j, k, l, lv; float xx, yy, zz; };
    template<class T> using vector_t = util::managed_vector<T>;

private: // data to output
    int time_minutes_now_ = -1;
    vector_t<st_ids> ids;
    vector_t<char> calcflg;
    vector_t<int> d_rank;
    vector_t<float> u, v, w;
    vector_t<float> levelset_obj;
    vector_t<float> scalars, T;
    vector_t<float> x, y, z;
    vector_t<float> xitp, yitp, zitp;

    vector_t<float> u_ave, v_ave, w_ave;
    vector_t<float> vel2_fluc_ave, T_ave;
    // ValueTimeAverage

    vector_t<size_t> heads_input;

public:
    FastPostprocess() = delete;
    FastPostprocess(const MPICommEnsemble comm): comm_(comm) {
        const std::string prefix = Foldernames::output_folder + "/" + out_prefix_ens_();
        ::mkdir(prefix.c_str(), 0755);
    }

    ~FastPostprocess() { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    size_t monitorcount() const noexcept { return x.size(); }

    void setupdata(const Tree&, const MeshValue* meshValues_ptr[]);
    void setupdata(const Tree& tree, const MeshValue (&meshValues)[DefAMR::LV_MAX]) { // compat
        const MeshValue* meshValues_ptr[DefAMR::LV_MAX];
        for(auto&& lv: util::irange(DefAMR::LV_MAX)) {
            meshValues_ptr[lv] = &meshValues[lv]; // &(meshValues_new(lv))
        }
        setupdata(tree, meshValues_ptr);
    }

    void OutputMonitorData(int step, int time_minutes, const  Tree&, const Parameters&, const MeshValue* meshValues_ptr[]
        , const ValueTimeAverage* valueTimeAverage1min 
        );

    void OutputMonitorData(int step, int time_minutes, const  Tree& tree, const Parameters& parameters, const MeshValue (&meshValues)[DefAMR::LV_MAX]
        , const ValueTimeAverage* valueTimeAverage1min
        ) { // compat
        const MeshValue* meshValues_ptr[DefAMR::LV_MAX];
        for(auto&& lv: util::irange(DefAMR::LV_MAX)) {
            meshValues_ptr[lv] = &meshValues[lv]; // &(meshValues_new(lv))
        }
        OutputMonitorData(step, time_minutes, tree, parameters, meshValues_ptr
            , valueTimeAverage1min
            );
    }



public: // trick for __device__ lambda
    void updateMonitorData(const Tree& tree, const Parameters& parameters, const MeshValue* meshValues_ptr[]
        , const ValueTimeAverage* valueTimeAverage1min
        );

    void writeMonitorData(int fileno, std::string filename, int step, bool trancate, const MeshValue* meshValues_ptr[]);

private:
    void readdata();
    void mallocdata(int n_scalars);

public: // to be private, but public is needed for cuda lambda...
    void setup_coordinate_to_cell(const Tree&, const MeshValue* meshValues_ptr[]);

};

#endif

