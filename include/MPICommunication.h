#pragma once
#ifndef MPICOMMUNICATION_H_
#define MPICOMMUNICATION_H_


#include <iostream>
#include <vector>
#include <mpi.h>
#include "defineAMR.h"
#include "defineLBM.h"
#include "MPIDatatypeInfo.h"
#include "Field.h"
#include "MeshValue.h"


class  MPICommunication {
private:
    int  rank_;
    int  nprocs_;

    int  halo_;
    int  val_outdated_[DefAMR::LV_MAX];
    int  lbm_outdated_[DefAMR::LV_MAX];

public:
    MPICommunication(){
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs_);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

        halo_ = DefAMR::NX_LEAF;
        for (int i=0; i<DefAMR::LV_MAX; i++) {
            val_outdated_[i] = halo_;
            lbm_outdated_[i] = halo_;
        }
    }
    ~MPICommunication(){}

public:
    bool ValCommPackUnpack(
        const int                       lv,
              Field&                    field, 
        const MPIPackUnpackInfo&        mpiPackUnpackInfo,
        const Tree&                     tree,
              ValueNS&                  valueNS,
              ValueObjLS&               valueObjLS,
              ValueBuff&                valueBuff
        );


    bool LBMCommPackUnpack(
        const int                       lv,
              Field&                    field, 
        const std::vector<MPIDatatypeInfoSrcDst>& lbmPutInfoSrcDst,
        const MPIPackUnpackInfo&        mpiPackUnpackInfo,
        const int                       nmax,
        const int                       nmax_global,
        const Tree&                     tree,
              ValueLBM&                 valueLBM,
              ValueBuff&                valueBuff
        );

private:
    void ValComm_pack_unpack(
        const int                       lv,
        const MPIPackUnpackInfo&        mpiPackUnpackInfo,
        const Tree&                     tree,
              real*                     val,
              ValueBuff&                valueBuff
        );


    void ValsComm_pack_unpack(
        const int                       lv,
        const MPIPackUnpackInfo&        mpiPackUnpackInfo,
        const Tree&                     tree,
              std::vector<real*>        vals,
              ValueBuff&                valueBuff
        );


    void LBMComm_pack_unpack(
        const int                       lv,
              Field&                    field, 
        const MPIPackUnpackInfo&        mpiPackUnpackInfo,
        const Tree&                     tree,
              ValueLBM&                 valueLBM,
              ValueBuff&                valueBuff
        );


    void LBMComm_pack_unpack_host(
        const int                       lv,
              Field&                    field, 
        const MPIPackUnpackInfo&        mpiPackUnpackInfo,
        const Tree&                     tree,
              ValueLBM&                 valueLBM,
              ValueBuff&                valueBuff
        );


    void clear_outdated(int& val);
    bool require_mpicomm(int lv, int val);


    // prefetch //
//    void cudaMemAdvise_ValComm(
//        const MPIDatatypeInfoSrcDst&  elem,
//        const Tree&                   tree,
//              ValueNS&                valueNS
//        );
//
//    void cudaMemAdvise_LBMComm(
//        const MPIDatatypeInfoSrcDst&  elem,
//        const Tree&                   tree,
//              ValueLBM&               valueLBM
//        );

    void
    cudaMemAdvise_val(
              real* val,
        const int n
        );

    void get_communication_time(
        const struct timeval& t_begin,
        const struct timeval& t_end,
        const struct timeval& t_begin_mpi,
        const struct timeval& t_end_mpi,
        const struct timeval& t_begin_pack,
        const struct timeval& t_end_pack,
        const struct timeval& t_begin_unpack,
        const struct timeval& t_end_unpack,
        const int   lv,
        const MPIPackUnpackInfo&    mpiPackUnpackInfo,
        const int   nQ
        );

};


#endif
