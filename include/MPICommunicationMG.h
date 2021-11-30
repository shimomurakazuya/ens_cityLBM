#pragma once
#ifndef MPICOMMUNICATIONMG_H_
#define MPICOMMUNICATIONMG_H_


#include <iostream>
#include <vector>
#include <mpi.h>
#include "defineAMR.h"
#include "MPIDatatypeInfo.h"
#include "Field.h"
#include "PackUnpackCPUGPU.h"
#include "MPICommEnsemble.h"
#include "mpi_wrapper.hpp"


class  MPICommunicationMG {
private:
    const MPICommEnsemble comm_;

    int  halo_;

public:
    MPICommunicationMG() = delete;

    MPICommunicationMG(const MPICommEnsemble comm): comm_(comm) {
        halo_ = DefAMR::NX_LEAF;
    }

    ~MPICommunicationMG(){}

private:
    void clear_outdated(int& val) { val = 0; }

    bool require_mpicomm(int lv, int val)
    {
        return  (val >= 1);
//    #ifdef TEMPORAL_BLOCKING_
//        return  (val >= halo_-1);
//    #else
//        return  (val >= 1);
//    #endif
    }

public:
    template<int NX_LEAF, typename T>
    void ValCommPackUnpackMG(
        const std::vector<MPIPackUnpackCommInfo>&   sendPackUnpackCommInfo,
        const std::vector<MPIPackUnpackCommInfo>&   recvPackUnpackCommInfo,
        const PackUnpackVal&                        packUnpackVal,
              std::vector<T*>        vals,
              ValueBuff&             valueBuff
        )
    {
        auto func_ovlp = [](){};
        ValCommPackUnpackMG<NX_LEAF>(
            func_ovlp,
            sendPackUnpackCommInfo, recvPackUnpackCommInfo,
            packUnpackVal,
            vals, valueBuff
            );
    }


    template<int NX_LEAF, typename T>
    void ValCommPackUnpackMG(
        const std::function<void()>                 func_ovlp,
        const std::vector<MPIPackUnpackCommInfo>&   sendPackUnpackCommInfo,
        const std::vector<MPIPackUnpackCommInfo>&   recvPackUnpackCommInfo,
        const PackUnpackVal&                        packUnpackVal,
              std::vector<T*>     vals,
              ValueBuff&          valueBuff
        )
    {
        const auto comm_target = comm_.col_vector().comm();

        T* sbuff  = valueBuff.sbuff_cal<T>();
        T* rbuff  = valueBuff.rbuff_cal<T>();
        const int  num_slist = packUnpackVal.num_slist_val;
        const int  num_rlist = packUnpackVal.num_rlist_val;
        const int* slist     = packUnpackVal.slist_val;
        const int* rlist     = packUnpackVal.rlist_val;

        // packing //
        int offset_slist = 0;
        for (auto& val : vals) {
#ifdef USE_NVCC
            if      (valueBuff.memType()==MemType::Device || valueBuff.memType()==MemType::Managed) { PackUnpack::pack_device<NX_LEAF>(slist, val, &sbuff[offset_slist], num_slist); }
            else if (valueBuff.memType()==MemType::Host)                                            { PackUnpack::pack_host           (slist, val, &sbuff[offset_slist], num_slist); }
#else
            PackUnpack::pack_host(slist, val, &sbuff[offset_slist], num_slist);
#endif
            offset_slist += num_slist;
        }

        // mpi //
        const int  num_rank_send  = sendPackUnpackCommInfo.size();
        const int  num_rank_recv  = recvPackUnpackCommInfo.size();

        // check
        const int _num_rank_send  = packUnpackVal.num_slist_val_rank.size();
        const int _num_rank_recv  = packUnpackVal.num_rlist_val_rank.size();
        if ( num_rank_send != _num_rank_send ) { std::cout << "num_rank_send != num_slist.size : " << num_rank_send << "!=" << _num_rank_send << std::endl; }
        if ( num_rank_recv != _num_rank_recv ) { std::cout << "num_rank_recv != num_rlist.size : " << num_rank_recv << "!=" << _num_rank_recv << std::endl; }

        const int  num_rank_send_recv = num_rank_send + num_rank_recv;
        MPI_Request request[num_rank_send_recv * vals.size()];
        MPI_Status  status [num_rank_send_recv * vals.size()];

        int irequ = 0;
        int offset_s = 0;
        int offset_r = 0;
        for (auto& val : vals) {
            // recv //
            for (int i=0; i<num_rank_recv; i++) {
                const int  rank_recv   = recvPackUnpackCommInfo[i].rank;
                const int  tag_recv    = recvPackUnpackCommInfo[i].tag;
                const int  num_recv    = packUnpackVal.num_rlist_val_rank[i];
                const int  offset_recv = packUnpackVal.offset_rlist_val_rank[i] + offset_r;

                MPI_Irecv( &rbuff[offset_recv], num_recv, util::mpi::MPItypename<T>::name(), rank_recv, tag_recv, comm_target, &request[irequ]);
                irequ++;
            }
            offset_r += num_rlist;

            // send //
            for (int i=0; i<num_rank_send; i++) {
                const int  rank_send   = sendPackUnpackCommInfo[i].rank;
                const int  tag_send    = sendPackUnpackCommInfo[i].tag;
                const int  num_send    = packUnpackVal.num_slist_val_rank[i];
                const int  offset_send = packUnpackVal.offset_slist_val_rank[i] + offset_s;

                MPI_Isend( &sbuff[offset_send], num_send, util::mpi::MPItypename<T>::name(), rank_send, tag_send, comm_target, &request[irequ]);
                irequ++;
            }
            offset_s += num_slist;
        }

        ///// ovlp /////
        func_ovlp();
        ///// ovlp /////

        MPI_Waitall(num_rank_send_recv * vals.size(), request, status);

        // unpacking //
        int offset_rlist = 0;
        for (auto& val : vals) {
#ifdef USE_NVCC
            if      (valueBuff.memType()==MemType::Device || valueBuff.memType()==MemType::Managed) { PackUnpack::unpack_device<NX_LEAF>(rlist, val, &rbuff[offset_rlist], num_rlist); }
            else if (valueBuff.memType()==MemType::Host)                                            { PackUnpack::unpack_host           (rlist, val, &rbuff[offset_rlist], num_rlist); }
#else
            PackUnpack::unpack_host(rlist, val, &rbuff[offset_rlist], num_rlist);
#endif
            offset_rlist += num_rlist;
        }
    }

};


#endif
