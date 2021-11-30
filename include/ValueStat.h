#pragma once
#ifndef VALUESTAT_H_
#define VALUESTAT_H_

#ifdef USE_VALUE_STAT

#include <iostream>
#include <cstdlib>
#include <cstdint>

#include "definePrecision.h"
#include "defineMemory.h"
#include "Grid.h"
#include "MeshValue.h"


class  ValueStat {
public: 
    using this_type = ValueStat;
    using stat_type = double;
private:
    MemType memType_;

    std::intptr_t nn_max_;
    std::intptr_t t_count_;

    // vel
    stat_type*   u_sum_ = nullptr;
    stat_type*   v_sum_ = nullptr;
    stat_type*   w_sum_ = nullptr;
    // temp
    stat_type*   T_sum_ = nullptr;
    // en
    stat_type*   uu_sum_ = nullptr;
    stat_type*   vv_sum_ = nullptr;
    stat_type*   ww_sum_ = nullptr;
    stat_type*   uv_sum_ = nullptr;
    stat_type*   vw_sum_ = nullptr;
    stat_type*   wu_sum_ = nullptr;
    // TT
    stat_type*   TT_sum_ = nullptr;
    // T flux
    stat_type*   Tu_sum_ = nullptr;
    stat_type*   Tv_sum_ = nullptr;
    stat_type*   Tw_sum_ = nullptr;
    // pow3
    stat_type*   uuu_sum_ = nullptr;
    stat_type*   vvv_sum_ = nullptr;
    stat_type*   www_sum_ = nullptr;
    stat_type*   TTT_sum_ = nullptr;
    // en flux
    stat_type*   uuw_sum_ = nullptr;
    stat_type*   vvw_sum_ = nullptr;
    // scalar
    stat_type*   sc_sum_  = nullptr;
    stat_type*   sc_ssum_ = nullptr;

public:
    ValueStat(const MemType& memType = MemType::Managed): memType_(memType), nn_max_(0), t_count_(0) {}

    ~ValueStat () { release(); }

public:
    MemType  memType() const { return memType_; }

    const stat_type* u_sum()   const { return u_sum_  ; }
    const stat_type* v_sum()   const { return v_sum_  ; }
    const stat_type* w_sum()   const { return w_sum_  ; }
    const stat_type* T_sum()   const { return T_sum_  ; }
    const stat_type* uu_sum()  const { return uu_sum_ ; }
    const stat_type* vv_sum()  const { return vv_sum_ ; }
    const stat_type* ww_sum()  const { return ww_sum_ ; }
    const stat_type* uv_sum()  const { return uv_sum_ ; }
    const stat_type* vw_sum()  const { return vw_sum_ ; }
    const stat_type* wu_sum()  const { return wu_sum_ ; }
    const stat_type* TT_sum()  const { return TT_sum_ ; }
    const stat_type* Tu_sum()  const { return Tu_sum_ ; }
    const stat_type* Tv_sum()  const { return Tv_sum_ ; }
    const stat_type* Tw_sum()  const { return Tw_sum_ ; }
    const stat_type* uuu_sum() const { return uuu_sum_; }
    const stat_type* vvv_sum() const { return vvv_sum_; }
    const stat_type* www_sum() const { return www_sum_; }
    const stat_type* TTT_sum() const { return TTT_sum_; }
    const stat_type* uuw_sum() const { return uuw_sum_; }
    const stat_type* vvw_sum() const { return vvw_sum_; }
    const stat_type* sc_sum()  const { return sc_sum_; }
    const stat_type* sc_ssum() const { return sc_ssum_; }

    const std::intptr_t& nn_max() const  { return nn_max_; }
    const std::intptr_t& t_count() const { return t_count_; }

public:
    void  init(const std::intptr_t&  nn_max);
    void  zeroset();
    void  sumup(const MeshValue& other);
    void  clone_from(const this_type& other);

private:
    void  allocate();
    void  release();
};

#else /// with nvcc -U USE_VALUE_STAT

class ValueStat { 
    public: // do nothing, but no problem to compile
    using this_type = ValueStat;
    using stat_type = double;
    ValueStat(const MemType& memType = MemType::Managed) {}
    ~ValueStat () {}
    void  init(const std::intptr_t&  nn_max) {}
    void  zeroset() {}
    void  sumup(const MeshValue& other) {}
    void  clone_from(const this_type& other) {}
    std::intptr_t nn_max() const  { return -1; }
    std::intptr_t t_count() const { return -1; }
    public: // compilation error if they are called within -U USE_VALUE_STAT defined
    const stat_type* u_sum()   const = delete;
    const stat_type* v_sum()   const = delete;
    const stat_type* w_sum()   const = delete;
    const stat_type* T_sum()   const = delete;
    const stat_type* uu_sum()  const = delete;
    const stat_type* vv_sum()  const = delete;
    const stat_type* ww_sum()  const = delete;
    const stat_type* uv_sum()  const = delete;
    const stat_type* vw_sum()  const = delete;
    const stat_type* wu_sum()  const = delete;
    const stat_type* TT_sum()  const = delete;
    const stat_type* Tu_sum()  const = delete;
    const stat_type* Tv_sum()  const = delete;
    const stat_type* Tw_sum()  const = delete;
    const stat_type* uuu_sum() const = delete;
    const stat_type* vvv_sum() const = delete;
    const stat_type* www_sum() const = delete;
    const stat_type* TTT_sum() const = delete;
    const stat_type* uuw_sum() const = delete;
    const stat_type* vvw_sum() const = delete;
    const stat_type* sc_sum()  const = delete;
    const stat_type* sc_ssum() const = delete;
};

#endif

#endif
