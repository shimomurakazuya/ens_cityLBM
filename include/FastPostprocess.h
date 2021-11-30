#ifndef FASTPOSTPROCESS_H_
#define FASTPOSTPROCESS_H_

#include "defineAMR.h"
#include "Tree.h"
#include "MeshValue.h"
#include "ValueStat.h"
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
    std::string out_filename_(int step, int fileno) const {
        std::string ret = out_prefix_ens_() + "/output";
        if(step >= 0) { ret += "-timestep" + std::to_string(step); }
        ret += "-station" + std::to_string(fileno);
        ret += ".csv";
        return ret;
    }

public:
    struct st_ids { ssize_t i, j, k, l, lv; float xx, yy, zz; };
    template<class T> using vector_t = util::managed_vector<T>;

private: // data to output
    vector_t<st_ids> ids;
    vector_t<char> calcflg;
    vector_t<int> d_rank;
    vector_t<float> u, v, w;
    vector_t<float> levelset_obj;
    vector_t<float> scalars, T;
    vector_t<float> x, y, z;
    vector_t<float> xitp, yitp, zitp;

#ifdef USE_VALUE_STAT
    vector_t<float> u_stat, v_stat, w_stat;
    vector_t<float> sc_stat, T_stat;
    vector_t<float> uu_stat, vv_stat, ww_stat, TT_stat;
#endif // USE_VALUE_STAT

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

    void OutputMonitorData(int step, const  Tree&, const Parameters&, const MeshValue* meshValues_ptr[]
        #ifdef USE_VALUE_STAT
        , const ValueStat* valueStat
        #endif // USE_VALUE_STAT
        );

    void OutputMonitorData(int step, const  Tree& tree, const Parameters& parameters, const MeshValue (&meshValues)[DefAMR::LV_MAX]
        #ifdef USE_VALUE_STAT
        , const ValueStat* valueStat
        #endif // USE_VALUE_STAT
        ) { // compat
        const MeshValue* meshValues_ptr[DefAMR::LV_MAX];
        for(auto&& lv: util::irange(DefAMR::LV_MAX)) {
            meshValues_ptr[lv] = &meshValues[lv]; // &(meshValues_new(lv))
        }
        OutputMonitorData(step, tree, parameters, meshValues_ptr
            #ifdef USE_VALUE_STAT
            , valueStat
            #endif // USE_VALUE_STAT
            );
    }


private:
    void readdata();
    void mallocdata(int n_scalars);

public: // to be private, but public is needed for cuda lambda...
    void setup_coordinate_to_cell(const Tree&, const MeshValue* meshValues_ptr[]);

};

#endif

