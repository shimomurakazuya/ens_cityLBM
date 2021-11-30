#ifdef USE_VALUE_STAT

#include "ValueStat.h"
#include "FuncAllocate.h"
#include "foreach.h"

// public //
void  ValueStat::
init(const std::intptr_t&  nn_max)
{
    if(nn_max_ <= 0) { release(); }
    nn_max_ = nn_max;
    allocate();
    zeroset();
}

void ValueStat::
zeroset() {
    t_count_ = 0;
    const auto nn_max = nn_max_;

    auto* u_sum_   = this->u_sum_  ;
    auto* v_sum_   = this->v_sum_  ;
    auto* w_sum_   = this->w_sum_  ;
    auto* T_sum_   = this->T_sum_  ;
    auto* uu_sum_  = this->uu_sum_ ;
    auto* vv_sum_  = this->vv_sum_ ;
    auto* ww_sum_  = this->ww_sum_ ;
    auto* uv_sum_  = this->uv_sum_ ;
    auto* vw_sum_  = this->vw_sum_ ;
    auto* wu_sum_  = this->wu_sum_ ;
    auto* TT_sum_  = this->TT_sum_ ;
    auto* Tu_sum_  = this->Tu_sum_ ;
    auto* Tv_sum_  = this->Tv_sum_ ;
    auto* Tw_sum_  = this->Tw_sum_ ;
    auto* uuu_sum_ = this->uuu_sum_;
    auto* vvv_sum_ = this->vvv_sum_;
    auto* www_sum_ = this->www_sum_;
    auto* TTT_sum_ = this->TTT_sum_;
    auto* uuw_sum_ = this->uuw_sum_;
    auto* vvw_sum_ = this->vvw_sum_;
    auto* sc_sum_ = this->sc_sum_;
    auto* sc_ssum_ = this->sc_sum_;

    foreach::exec_1d<foreach::opti>(
        nn_max,
        [=] __HD__ () {
            FOR_EACH1D_XX(i, nn_max) {
                u_sum_  [i] = static_cast<stat_type>(0);
                v_sum_  [i] = static_cast<stat_type>(0);
                w_sum_  [i] = static_cast<stat_type>(0);
                T_sum_  [i] = static_cast<stat_type>(0);
                uu_sum_ [i] = static_cast<stat_type>(0);
                vv_sum_ [i] = static_cast<stat_type>(0);
                ww_sum_ [i] = static_cast<stat_type>(0);
                uv_sum_ [i] = static_cast<stat_type>(0);
                vw_sum_ [i] = static_cast<stat_type>(0);
                wu_sum_ [i] = static_cast<stat_type>(0);
                TT_sum_ [i] = static_cast<stat_type>(0);
                Tu_sum_ [i] = static_cast<stat_type>(0);
                Tv_sum_ [i] = static_cast<stat_type>(0);
                Tw_sum_ [i] = static_cast<stat_type>(0);
                uuu_sum_[i] = static_cast<stat_type>(0);
                vvv_sum_[i] = static_cast<stat_type>(0);
                www_sum_[i] = static_cast<stat_type>(0);
                TTT_sum_[i] = static_cast<stat_type>(0);
                uuw_sum_[i] = static_cast<stat_type>(0);
                vvw_sum_[i] = static_cast<stat_type>(0);
                sc_sum_ [i] = static_cast<stat_type>(0);
                sc_ssum_[i] = static_cast<stat_type>(0);
            }
        }
    );
}

void ValueStat::
sumup(const MeshValue& other) 
{
    ++ t_count_;

    const auto nn_max = other.nn_max();
    if(nn_max != nn_max_) { init(nn_max); }

    const auto* u_ = other.valueNS().u();
    const auto* v_ = other.valueNS().v();
    const auto* w_ = other.valueNS().w();
    const auto* T_ = other.valueNS().T();

    auto* u_sum_   = this->u_sum_  ;
    auto* v_sum_   = this->v_sum_  ;
    auto* w_sum_   = this->w_sum_  ;
    auto* T_sum_   = this->T_sum_  ;
    auto* uu_sum_  = this->uu_sum_ ;
    auto* vv_sum_  = this->vv_sum_ ;
    auto* ww_sum_  = this->ww_sum_ ;
    auto* uv_sum_  = this->uv_sum_ ;
    auto* vw_sum_  = this->vw_sum_ ;
    auto* wu_sum_  = this->wu_sum_ ;
    auto* TT_sum_  = this->TT_sum_ ;
    auto* Tu_sum_  = this->Tu_sum_ ;
    auto* Tv_sum_  = this->Tv_sum_ ;
    auto* Tw_sum_  = this->Tw_sum_ ;
    auto* uuu_sum_ = this->uuu_sum_;
    auto* vvv_sum_ = this->vvv_sum_;
    auto* www_sum_ = this->www_sum_;
    auto* TTT_sum_ = this->TTT_sum_;
    auto* uuw_sum_ = this->uuw_sum_;
    auto* vvw_sum_ = this->vvw_sum_;

    const auto* sc_ = other.valueNS().scalar();
    auto* sc_sum_ = this->sc_sum_;
    auto* sc_ssum_ = this->sc_sum_;

    foreach::exec_1d<foreach::opti>(
        nn_max,
        [=] __HD__ () {
            FOR_EACH1D_XX(i, nn_max) {
                const auto u = u_[i];
                const auto v = v_[i];
                const auto w = w_[i];
                const auto T = T_[i];
                u_sum_  [i] += u;
                v_sum_  [i] += v;
                w_sum_  [i] += w;
                T_sum_  [i] += T;
                uu_sum_ [i] += u*u;
                vv_sum_ [i] += v*v;
                ww_sum_ [i] += w*w;
                uv_sum_ [i] += u*v;
                vw_sum_ [i] += v*w;
                wu_sum_ [i] += w*u;
                TT_sum_ [i] += T*T;
                Tu_sum_ [i] += T*u;
                Tv_sum_ [i] += T*v;
                Tw_sum_ [i] += T*w;
                uuu_sum_[i] += u*u*u;
                vvv_sum_[i] += v*v*v;
                www_sum_[i] += w*w*w;
                TTT_sum_[i] += T*T*T;
                uuw_sum_[i] += u*u*w;
                vvw_sum_[i] += v*v*w;
                const auto sc = sc_[i];
                sc_sum_[i] += sc;
                sc_ssum_[i] += sc*sc;
            }
        }
    );
}


void  ValueStat::
clone_from(const this_type&  other)
{
    if(other.nn_max_ != nn_max_) { init(other.nn_max_); }
    t_count_ = other.t_count_;
    memType_ = other.memType_;
    FuncAllocate::copy_values(u_sum_  , other.u_sum_  , nn_max_, memType_);
    FuncAllocate::copy_values(v_sum_  , other.v_sum_  , nn_max_, memType_);
    FuncAllocate::copy_values(w_sum_  , other.w_sum_  , nn_max_, memType_);
    FuncAllocate::copy_values(T_sum_  , other.T_sum_  , nn_max_, memType_);
    FuncAllocate::copy_values(uu_sum_ , other.uu_sum_ , nn_max_, memType_);
    FuncAllocate::copy_values(vv_sum_ , other.vv_sum_ , nn_max_, memType_);
    FuncAllocate::copy_values(ww_sum_ , other.ww_sum_ , nn_max_, memType_);
    FuncAllocate::copy_values(uv_sum_ , other.uv_sum_ , nn_max_, memType_);
    FuncAllocate::copy_values(vw_sum_ , other.vw_sum_ , nn_max_, memType_);
    FuncAllocate::copy_values(wu_sum_ , other.wu_sum_ , nn_max_, memType_);
    FuncAllocate::copy_values(TT_sum_ , other.TT_sum_ , nn_max_, memType_);
    FuncAllocate::copy_values(Tu_sum_ , other.Tu_sum_ , nn_max_, memType_);
    FuncAllocate::copy_values(Tv_sum_ , other.Tv_sum_ , nn_max_, memType_);
    FuncAllocate::copy_values(Tw_sum_ , other.Tw_sum_ , nn_max_, memType_);
    FuncAllocate::copy_values(uuu_sum_, other.uuu_sum_, nn_max_, memType_);
    FuncAllocate::copy_values(vvv_sum_, other.vvv_sum_, nn_max_, memType_);
    FuncAllocate::copy_values(www_sum_, other.www_sum_, nn_max_, memType_);
    FuncAllocate::copy_values(TTT_sum_, other.TTT_sum_, nn_max_, memType_);
    FuncAllocate::copy_values(uuw_sum_, other.uuw_sum_, nn_max_, memType_);
    FuncAllocate::copy_values(vvw_sum_, other.vvw_sum_, nn_max_, memType_);
    FuncAllocate::copy_values(sc_sum_ , other.sc_sum_ , nn_max_, memType_);
    FuncAllocate::copy_values(sc_ssum_, other.sc_ssum_, nn_max_, memType_);
}


// private //
void  ValueStat::
allocate()
{
    FuncAllocate::allocate_value<stat_type>(&u_sum_  ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&v_sum_  ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&w_sum_  ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&T_sum_  ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&uu_sum_ ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&vv_sum_ ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&ww_sum_ ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&uv_sum_ ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&vw_sum_ ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&wu_sum_ ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&TT_sum_ ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&Tu_sum_ ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&Tv_sum_ ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&Tw_sum_ ,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&uuu_sum_,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&vvv_sum_,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&www_sum_,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&TTT_sum_,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&uuw_sum_,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&vvw_sum_,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&sc_sum_,  nn_max_, memType_);
    FuncAllocate::allocate_value<stat_type>(&sc_ssum_,  nn_max_, memType_);
}


void  ValueStat::
release()
{
    FuncAllocate::release_value<stat_type>(u_sum_  ,  memType_);
    FuncAllocate::release_value<stat_type>(v_sum_  ,  memType_);
    FuncAllocate::release_value<stat_type>(w_sum_  ,  memType_);
    FuncAllocate::release_value<stat_type>(T_sum_  ,  memType_);
    FuncAllocate::release_value<stat_type>(uu_sum_ ,  memType_);
    FuncAllocate::release_value<stat_type>(vv_sum_ ,  memType_);
    FuncAllocate::release_value<stat_type>(ww_sum_ ,  memType_);
    FuncAllocate::release_value<stat_type>(uv_sum_ ,  memType_);
    FuncAllocate::release_value<stat_type>(vw_sum_ ,  memType_);
    FuncAllocate::release_value<stat_type>(wu_sum_ ,  memType_);
    FuncAllocate::release_value<stat_type>(TT_sum_ ,  memType_);
    FuncAllocate::release_value<stat_type>(Tu_sum_ ,  memType_);
    FuncAllocate::release_value<stat_type>(Tv_sum_ ,  memType_);
    FuncAllocate::release_value<stat_type>(Tw_sum_ ,  memType_);
    FuncAllocate::release_value<stat_type>(uuu_sum_,  memType_);
    FuncAllocate::release_value<stat_type>(vvv_sum_,  memType_);
    FuncAllocate::release_value<stat_type>(www_sum_,  memType_);
    FuncAllocate::release_value<stat_type>(TTT_sum_,  memType_);
    FuncAllocate::release_value<stat_type>(uuw_sum_,  memType_);
    FuncAllocate::release_value<stat_type>(vvw_sum_,  memType_);
    FuncAllocate::release_value<stat_type>(sc_sum_,  memType_);
    FuncAllocate::release_value<stat_type>(sc_ssum_,  memType_);
    nn_max_ = 0;
}

#endif

