#include "MPICommunication.h"
#include "defineMemory.h"
#include "definePrecision.h"
#include "PackUnpackCPUGPU.h"
#include "FuncMath.h"
#include <sys/time.h>


bool MPICommunication::TCommPackUnpack(
    const int                       lv,
          Field&                    field,
    const MPIPackUnpackInfo&        mpiPackUnpackInfo,
    const Tree&                     tree,
          ValueNS&                  valueNS,
          ValueObjLS&               valueObjLS,
          ValueBuff&                valueBuff
    )
{
    // update //
    T_outdated_[lv] = T_outdated_[lv] + n_T_outdated_;
    if ( !require_mpicomm(lv, T_outdated_[lv]) ) { return false; } // tb


    // MPI Comm //
    const int nn_max = field.meshValue(lv).nn_max();

    std::vector<real*>  vals;
    vals.push_back(valueNS.T());
    ValsComm_pack_unpack(lv, mpiPackUnpackInfo, tree, vals, valueBuff);


    // clear //
    clear_outdated(T_outdated_[lv]);

    return true;
}


bool MPICommunication::ScalarCommPackUnpack(
    const int                       lv,
          Field&                    field,
    const MPIPackUnpackInfo&        mpiPackUnpackInfo,
    const Tree&                     tree,
          ValueNS&                  valueNS,
          ValueObjLS&               valueObjLS,
          ValueBuff&                valueBuff
    )
{
    // update //
    scalar_outdated_[lv] = scalar_outdated_[lv] + n_scalar_outdated_;
    if ( !require_mpicomm(lv, scalar_outdated_[lv]) ) { return false; } // tb


    // MPI Comm //
    const int nn_max = field.meshValue(lv).nn_max();
    const int n_scalars = field.meshValue(lv).n_scalars();

    std::vector<real*>  vals;
    for(int n=0; n<n_scalars; n++) {
        int idx = n * nn_max;
        vals.push_back(&(valueNS.scalar()[idx]));
    }
    ValsComm_pack_unpack(lv, mpiPackUnpackInfo, tree, vals, valueBuff);


    // clear //
    clear_outdated(scalar_outdated_[lv]);

    return true;
}


bool MPICommunication::LBMCommPackUnpack(
    const int                       lv,
          Field&                    field,
    const std::vector<MPIDatatypeInfoSrcDst>& lbmPutInfoSrcDst,
    const MPIPackUnpackInfo&        mpiPackUnpackInfo,
    const int                       nmax,
    const int                       nmax_global,
    const Tree&                     tree,
          ValueNS&                  valueNS,
          ValueLBM&                 valueLBM,
          ValueBuff&                valueBuff
    )
{
    // update //
    lbm_outdated_[lv] = lbm_outdated_[lv] + n_lbm_outdated_;
    if ( !require_mpicomm(lv, lbm_outdated_[lv]) ) { return false; } // tb


    // MPI Comm //
    LBMComm_pack_unpack(lv, field, mpiPackUnpackInfo, tree, valueLBM, valueBuff);

    #ifdef USE_TEMPORAL_BLOCKING
    ValsComm_pack_unpack(lv, mpiPackUnpackInfo, tree, { valueNS.sgs_vis() }, valueBuff);
    #endif

    // clear //
    clear_outdated(lbm_outdated_[lv]);

    return true;
}


void MPICommunication::ValComm_pack_unpack(
    const int                       lv,
    const MPIPackUnpackInfo&        mpiPackUnpackInfo,
    const Tree&                     tree,
          real*                     val,
          ValueBuff&                valueBuff
    )
{
    const auto comm_target = comm_.col_vector().comm();
    real* sbuff  = valueBuff.sbuff_cal<real>();
    real* rbuff  = valueBuff.rbuff_cal<real>();

    const int  num_slist = mpiPackUnpackInfo.num_slist_val;
    const int  num_rlist = mpiPackUnpackInfo.num_rlist_val;
    const int* slist     = mpiPackUnpackInfo.slist_val;
    const int* rlist     = mpiPackUnpackInfo.rlist_val;


    // packing //
//    PackUnpack::pack<real>(slist, val, sbuff, num_slist);
    if (valueBuff.memType() == MemType::Device || valueBuff.memType() == MemType::Managed) {
#ifdef USE_NVCC
        PackUnpack::pack_device<DefAMR::NX_LEAF>(slist, val, sbuff, num_slist);
#endif
    }
    else if (valueBuff.memType() == MemType::Host) {
        PackUnpack::pack_host(slist, val, sbuff, num_slist);
    }


    // mpi //
    const int  num_rank_send  = mpiPackUnpackInfo.sendPackUnpackCommInfo.size();
    const int  num_rank_recv  = mpiPackUnpackInfo.recvPackUnpackCommInfo.size();

    const int  num_rank_send_recv = num_rank_send + num_rank_recv;
    MPI_Request request[num_rank_send_recv];
    MPI_Status  status [num_rank_send_recv];

    // recv //
    for (int i=0; i<num_rank_recv; i++) {
        const int  rank_recv   = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].rank;
        const int  tag_recv    = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].tag;
        const int  num_recv    = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].nleaf    * DefAMR::NN_LEAF;
        const int  offset_recv = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].sum_leaf * DefAMR::NN_LEAF;

        MPI_Irecv( &rbuff[offset_recv], num_recv, MFLOAT, rank_recv, tag_recv, comm_target, &request[i]);
    }

    // send //
    for (int i=0; i<num_rank_send; i++) {
        const int  rank_send   = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].rank;
        const int  tag_send    = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].tag;
        const int  num_send    = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf    * DefAMR::NN_LEAF;
        const int  offset_send = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].sum_leaf * DefAMR::NN_LEAF;

        MPI_Isend( &sbuff[offset_send], num_send, MFLOAT, rank_send, tag_send, comm_target, &request[num_rank_recv + i]);
    }

    // wait //
    MPI_Waitall(num_rank_send_recv, request, status);


    // unpacking //
//    PackUnpack::unpack<real>(rlist, val, rbuff, num_rlist);
    if (valueBuff.memType() == MemType::Device || valueBuff.memType() == MemType::Managed) {
#ifdef USE_NVCC
        PackUnpack::unpack_device<DefAMR::NX_LEAF>(rlist, val, rbuff, num_rlist);
#endif
    }
    else if (valueBuff.memType() == MemType::Host) {
        PackUnpack::unpack_host(rlist, val, rbuff, num_rlist);
    }
}


void MPICommunication::ValsComm_pack_unpack(
    const int                       lv,
    const MPIPackUnpackInfo&        mpiPackUnpackInfo,
    const Tree&                     tree,
          std::vector<real*>        vals,
          ValueBuff&                valueBuff
    )
{
    const auto comm_target = comm_.col_vector().comm();
    real* sbuff  = valueBuff.sbuff_cal<real>();
    real* rbuff  = valueBuff.rbuff_cal<real>();

    const int  num_slist = mpiPackUnpackInfo.num_slist_val;
    const int  num_rlist = mpiPackUnpackInfo.num_rlist_val;
    const int* slist     = mpiPackUnpackInfo.slist_val;
    const int* rlist     = mpiPackUnpackInfo.rlist_val;

    // packing //
    int offset_slist = 0;
    for (auto& val : vals) {
        if (valueBuff.memType() == MemType::Device || valueBuff.memType() == MemType::Managed) {
#ifdef USE_NVCC
            PackUnpack::pack_device<DefAMR::NX_LEAF>(slist, val, &sbuff[offset_slist], num_slist);
#endif
        }
        else if (valueBuff.memType() == MemType::Host) {
            PackUnpack::pack_host(slist, val, &sbuff[offset_slist], num_slist);
        }
        offset_slist += num_slist;
    }


    // mpi //
    const int  num_rank_send  = mpiPackUnpackInfo.sendPackUnpackCommInfo.size();
    const int  num_rank_recv  = mpiPackUnpackInfo.recvPackUnpackCommInfo.size();

    const int  num_rank_send_recv = num_rank_send + num_rank_recv;
    MPI_Request request[num_rank_send_recv * vals.size()];
    MPI_Status  status [num_rank_send_recv * vals.size()];

    int irequ = 0;
    int offset_s = 0;
    int offset_r = 0;
    for (auto& val : vals) {
        // recv //
        for (int i=0; i<num_rank_recv; i++) {
            const int  rank_recv   = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].rank;
            const int  tag_recv    = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].tag;
            const int  num_recv    = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].nleaf    * DefAMR::NN_LEAF;
            const int  offset_recv = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].sum_leaf * DefAMR::NN_LEAF + offset_r;

            MPI_Irecv( &rbuff[offset_recv], num_recv, MFLOAT, rank_recv, tag_recv, comm_target, &request[irequ]);
            irequ++;
        }
        offset_r += num_rlist;

        // send //
        for (int i=0; i<num_rank_send; i++) {
            const int  rank_send   = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].rank;
            const int  tag_send    = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].tag;
            const int  num_send    = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf    * DefAMR::NN_LEAF;
            const int  offset_send = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].sum_leaf * DefAMR::NN_LEAF + offset_s;

            MPI_Isend( &sbuff[offset_send], num_send, MFLOAT, rank_send, tag_send, comm_target, &request[irequ]);
            irequ++;
        }
        offset_s += num_slist;
    }

    // wait //
    MPI_Waitall(num_rank_send_recv * vals.size(), request, status);


    // unpacking //
    int offset_rlist = 0;
    for (auto& val : vals) {
        if (valueBuff.memType() == MemType::Device || valueBuff.memType() == MemType::Managed) {
#ifdef USE_NVCC
            PackUnpack::unpack_device<DefAMR::NX_LEAF>(rlist, val, &rbuff[offset_rlist], num_rlist);
#endif
        }
        else if (valueBuff.memType() == MemType::Host) {
            PackUnpack::unpack_host(rlist, val, &rbuff[offset_rlist], num_rlist);
        }
        offset_rlist += num_rlist;
    }
}


void MPICommunication::LBMComm_pack_unpack(
    const int                       lv,
          Field&                    field,
    const MPIPackUnpackInfo&        mpiPackUnpackInfo,
    const Tree&                     tree,
          ValueLBM&                 valueLBM,
          ValueBuff&                valueBuff
    )
{
    const auto comm_target = comm_.col_vector().comm();
    real* f_lbm = valueLBM.f_lbm();
    real* sbuff = valueBuff.sbuff_cal<real>();
    real* rbuff = valueBuff.rbuff_cal<real>();

    const int  num_slist = mpiPackUnpackInfo.num_slist_lbm;
    const int  num_rlist = mpiPackUnpackInfo.num_rlist_lbm;
    const int* slist_lbm = mpiPackUnpackInfo.slist_lbm;
    const int* rlist_lbm = mpiPackUnpackInfo.rlist_lbm;

//#define CHECK_MPI_TIME

#ifdef CHECK_MPI_TIME
    struct timeval  t_begin, t_begin_mpi, t_begin_pack, t_begin_unpack;
    struct timeval  t_end,   t_end_mpi,   t_end_pack,   t_end_unpack;
    constexpr int ncout = 500;


    gettimeofday(&t_begin, NULL);
    gettimeofday(&t_begin_pack, NULL);
#endif

    field.mpiTimeInfo().StartTimer("LBM_pack"); // timer //

    // packing //
//    PackUnpack::pack<real>(slist_lbm, f_lbm, sbuff, num_slist);
    if (valueBuff.memType() == MemType::Device || valueBuff.memType() == MemType::Managed) {
#ifdef USE_NVCC
        PackUnpack::pack_device<DefAMR::NX_LEAF>(slist_lbm, f_lbm, sbuff, num_slist);
#endif
    }
    else if (valueBuff.memType() == MemType::Host) {
        PackUnpack::pack_host(slist_lbm, f_lbm, sbuff, num_slist);
    }

    field.mpiTimeInfo().StopTimer("LBM_pack"); // timer //
    field.mpiTimeInfo().SubmitElapsedTimeInfo("LBM_pack", lv, num_slist, field.parameters().step_now()); // timer //

#ifdef CHECK_MPI_TIME
    gettimeofday(&t_end_pack, NULL);
    gettimeofday(&t_begin_mpi, NULL);
#endif

    field.mpiTimeInfo().StartTimer("LBM_mpi"); // timer //

    // mpi //
    constexpr int nQ = LBM_velocity_model::nQ;

    const int  num_rank_send  = mpiPackUnpackInfo.sendPackUnpackCommInfo.size();
    const int  num_rank_recv  = mpiPackUnpackInfo.recvPackUnpackCommInfo.size();

    const int  num_rank_send_recv = num_rank_send + num_rank_recv;
    MPI_Request request[num_rank_send_recv];
    MPI_Status  status [num_rank_send_recv];

    // recv //
    for (int i=0; i<num_rank_recv; i++) {
        const int  rank_recv   = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].rank;
        const int  tag_recv    = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].tag;
        const int  num_recv    = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].nleaf    * DefAMR::NN_LEAF*nQ;
        const int  offset_recv = mpiPackUnpackInfo.recvPackUnpackCommInfo[i].sum_leaf * DefAMR::NN_LEAF*nQ;

        MPI_Irecv( &rbuff[offset_recv], num_recv, MFLOAT, rank_recv, tag_recv, comm_target, &request[i]);
    }

    // send //
    for (int i=0; i<num_rank_send; i++) {
        const int  rank_send   = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].rank;
        const int  tag_send    = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].tag;
        const int  num_send    = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf    * DefAMR::NN_LEAF*nQ;
        const int  offset_send = mpiPackUnpackInfo.sendPackUnpackCommInfo[i].sum_leaf * DefAMR::NN_LEAF*nQ;

        MPI_Isend( &sbuff[offset_send], num_send, MFLOAT, rank_send, tag_send, comm_target, &request[num_rank_recv+i]);
    }


    // wait //
    MPI_Waitall(num_rank_send_recv, request, status);

    field.mpiTimeInfo().StopTimer("LBM_mpi"); // timer //
    field.mpiTimeInfo().SubmitElapsedTimeInfo("LBM_mpi", lv, num_slist, field.parameters().step_now()); // timer //

#ifdef CHECK_MPI_TIME
    gettimeofday(&t_end_mpi, NULL);
    gettimeofday(&t_begin_unpack, NULL);
#endif

    field.mpiTimeInfo().StartTimer("LBM_unpack"); // timer //

    // unpacking //
//    PackUnpack::unpack<real>(rlist_lbm, f_lbm, rbuff, num_rlist);
    if (valueBuff.memType() == MemType::Device || valueBuff.memType() == MemType::Managed) {
#ifdef USE_NVCC
        PackUnpack::unpack_device<DefAMR::NX_LEAF>(rlist_lbm, f_lbm, rbuff, num_rlist);
#endif
    }
    else if (valueBuff.memType() == MemType::Host) {
        PackUnpack::unpack_host(rlist_lbm, f_lbm, rbuff, num_rlist);
    }

    field.mpiTimeInfo().StopTimer("LBM_unpack"); // timer //
    field.mpiTimeInfo().SubmitElapsedTimeInfo("LBM_unpack", lv, num_slist, field.parameters().step_now()); // timer //

#ifdef CHECK_MPI_TIME
    gettimeofday(&t_end_unpack, NULL);
    gettimeofday(&t_end, NULL);
#endif

#ifdef CHECK_MPI_TIME
    static int  count;
    if (count%ncout == 0) {
        get_communication_time(
            t_begin,        t_end,
            t_begin_mpi,    t_end_mpi,
            t_begin_pack,   t_end_pack,
            t_begin_unpack, t_end_unpack,
            lv, mpiPackUnpackInfo, nQ );
    }
    count++;
#endif

#ifdef  CHECK_MPI_TIME
#undef  CHECK_MPI_TIME
#endif
}


void MPICommunication::clear_outdated(int& val)
{
    val = 0;
}


bool MPICommunication::require_mpicomm(int lv, int val)
{
#ifdef USE_TEMPORAL_BLOCKING
//    return  (DefAMR::is_level_max) ? (val >= halo_) : (val >= halo_-1);
    return  (val >= halo_-1); // tb : interpolation + sgs

//    return  (val >= halo_); // tb
//    return  (val >= halo_-1); // tb : interpolation
#else
//    return  (val >= halo_-1); // tb : interpolation
//    return  (val >= halo_-2);
//    return  (val >= 2); // normal
    return  (val >= 1); // normal
#endif
}


// prefetch //
//void MPICommunication::
//cudaMemAdvise_ValComm(
//    const MPIDatatypeInfoSrcDst&  elem,
//    const Tree&                   tree,
//          ValueNS&                valueNS
//    )
//{
//    real* scalar = valueNS.scalar();
//
//    int device_id; cudaGetDevice(&device_id);
//    constexpr int nput = DefAMR::NN_LEAF;
//    for (int i=0; i<(int)elem.mpiDatatypeInfo_dst.id.size(); i++) {
//        const int  id_src   = elem.mpiDatatypeInfo_src.id[i];
//
//        cudaMemAdvise(&scalar[id_src], nput, cudaMemAdviseSetReadMostly, device_id);
//    }
//}
//
//
//void MPICommunication::
//cudaMemAdvise_LBMComm(
//    const MPIDatatypeInfoSrcDst&  elem,
//    const Tree&                   tree,
//          ValueLBM&               valueLBM
//    )
//{
//    constexpr int nQ = LBM_velocity_model::nQ;
//
//    real* f_lbm = valueLBM.f_lbm();
//
//    int device_id; cudaGetDevice(&device_id);
//    constexpr int nput = DefAMR::NN_LEAF * nQ;
//    for (int i=0; i<(int)elem.mpiDatatypeInfo_dst.id.size(); i++) {
//        const int  id_src   = elem.mpiDatatypeInfo_src.id[i];
//
//        cudaMemAdvise(&f_lbm[id_src], nput, cudaMemAdviseSetReadMostly, device_id);
//    }
//}


//void MPICommunication::
//cudaMemAdvise_LBMComm(
//    const MPIDatatypeInfoSrcDst&  elem,
//    const Tree&                   tree,
//          ValueLBM&               valueLBM
//    )
//{
//    constexpr int nQ = LBM_velocity_model::nQ;
//
//    real* f_lbm = valueLBM.f_lbm();
//
//    int device_id; cudaGetDevice(&device_id);
//    constexpr int nput = DefAMR::NN_LEAF * nQ;
//    for (int i=0; i<(int)elem.mpiDatatypeInfo_dst.id.size(); i++) {
//        const int  id_src   = elem.mpiDatatypeInfo_src.id[i];
//
//        cudaMemAdvise(&f_lbm[id_src], nput, cudaMemAdviseSetReadMostly, device_id);
//    }
//}


void MPICommunication::
cudaMemAdvise_val(
          real* val,
    const int n
    )
{
    #ifdef GPU_CALCULATION__
      #if defined(USE_NVCC)
        int device_id; cudaGetDevice(&device_id);
        cudaMemAdvise(val, sizeof(real)*n, cudaMemAdviseSetReadMostly, device_id);
      #endif
    #endif
}


void MPICommunication::get_communication_time(
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
    )
{
    const auto comm_target = comm_.col_vector().comm();
    const int  num_rank_send  = mpiPackUnpackInfo.sendPackUnpackCommInfo.size();
    const int  num_rank_recv  = mpiPackUnpackInfo.recvPackUnpackCommInfo.size();

    const float msec_total  = (t_end       .tv_sec  - t_begin       .tv_sec) * 1000 + (t_end       .tv_usec - t_begin       .tv_usec) / 1000.0;
    const float msec_mpi    = (t_end_mpi   .tv_sec  - t_begin_mpi   .tv_sec) * 1000 + (t_end_mpi   .tv_usec - t_begin_mpi   .tv_usec) / 1000.0;
    const float msec_pack   = (t_end_pack  .tv_sec  - t_begin_pack  .tv_sec) * 1000 + (t_end_pack  .tv_usec - t_begin_pack  .tv_usec) / 1000.0;
    const float msec_unpack = (t_end_unpack.tv_sec  - t_begin_unpack.tv_sec) * 1000 + (t_end_unpack.tv_usec - t_begin_unpack.tv_usec) / 1000.0;

    int  nleaf_s = 0;
    int  nleaf_r = 0;
    for (int i=0; i<num_rank_send; i++) { nleaf_s += mpiPackUnpackInfo.sendPackUnpackCommInfo[i].nleaf    * DefAMR::NN_LEAF*nQ; }
    for (int i=0; i<num_rank_recv; i++) { nleaf_r += mpiPackUnpackInfo.recvPackUnpackCommInfo[i].nleaf    * DefAMR::NN_LEAF*nQ; }

    std::cout << std::endl << std::flush;
    MPI_Barrier(MPI_COMM_WORLD);
    const auto ncpu = comm_.col_vector().size(); 

    int nleaf_s_max = 0;          int nleaf_r_max = 0;
    int nleaf_s_min = 100000000;  int nleaf_r_min = 100000000;
    int nleaf_s_sum = 0;          int nleaf_r_sum = 0;

    float msec_total_max = 0.0;     float msec_mpi_max = 0.0;     float msec_pack_max = 0.0;     float msec_unpack_max = 0.0;
    float msec_total_min = 1.0e16;  float msec_mpi_min = 1.0e16;  float msec_pack_min = 1.0e16;  float msec_unpack_min = 1.0e16;
    float msec_total_sum = 0.0;     float msec_mpi_sum = 0.0;     float msec_pack_sum = 0.0;     float msec_unpack_sum = 0.0;

    MPI_Request requ_cout[18];
    MPI_Iallreduce(&nleaf_s,  &nleaf_s_max,  1, MPI_INT, MPI_MAX, comm_target, &requ_cout[0]);
    MPI_Iallreduce(&nleaf_r,  &nleaf_r_max,  1, MPI_INT, MPI_MAX, comm_target, &requ_cout[1]);
    MPI_Iallreduce(&nleaf_s,  &nleaf_s_min,  1, MPI_INT, MPI_MIN, comm_target, &requ_cout[2]);
    MPI_Iallreduce(&nleaf_r,  &nleaf_r_min,  1, MPI_INT, MPI_MIN, comm_target, &requ_cout[3]);
    MPI_Iallreduce(&nleaf_s,  &nleaf_s_sum,  1, MPI_INT, MPI_SUM, comm_target, &requ_cout[4]);
    MPI_Iallreduce(&nleaf_r,  &nleaf_r_sum,  1, MPI_INT, MPI_SUM, comm_target, &requ_cout[5]);

    MPI_Iallreduce(&msec_total,  &msec_total_max,  1, MPI_FLOAT, MPI_MAX, comm_target, &requ_cout[6]);
    MPI_Iallreduce(&msec_total,  &msec_total_min,  1, MPI_FLOAT, MPI_MIN, comm_target, &requ_cout[7]);
    MPI_Iallreduce(&msec_total,  &msec_total_sum,  1, MPI_FLOAT, MPI_SUM, comm_target, &requ_cout[8]);

    MPI_Iallreduce(&msec_mpi,  &msec_mpi_max,  1, MPI_FLOAT, MPI_MAX, comm_target, &requ_cout[9]);
    MPI_Iallreduce(&msec_mpi,  &msec_mpi_min,  1, MPI_FLOAT, MPI_MIN, comm_target, &requ_cout[10]);
    MPI_Iallreduce(&msec_mpi,  &msec_mpi_sum,  1, MPI_FLOAT, MPI_SUM, comm_target, &requ_cout[11]);

    MPI_Iallreduce(&msec_pack,  &msec_pack_max,  1, MPI_FLOAT, MPI_MAX, comm_target, &requ_cout[12]);
    MPI_Iallreduce(&msec_pack,  &msec_pack_min,  1, MPI_FLOAT, MPI_MIN, comm_target, &requ_cout[13]);
    MPI_Iallreduce(&msec_pack,  &msec_pack_sum,  1, MPI_FLOAT, MPI_SUM, comm_target, &requ_cout[14]);

    MPI_Iallreduce(&msec_unpack,  &msec_unpack_max,  1, MPI_FLOAT, MPI_MAX, comm_target, &requ_cout[15]);
    MPI_Iallreduce(&msec_unpack,  &msec_unpack_min,  1, MPI_FLOAT, MPI_MIN, comm_target, &requ_cout[16]);
    MPI_Iallreduce(&msec_unpack,  &msec_unpack_sum,  1, MPI_FLOAT, MPI_SUM, comm_target, &requ_cout[17]);
    MPI_Waitall(18, requ_cout, MPI_STATUSES_IGNORE);

    const int  nleaf_s_mean = nleaf_s_sum/ncpu;
    const int  nleaf_r_mean = nleaf_r_sum/ncpu;
//    const float  msec_total_mean    = msec_total_sum   /ncpu;
    const float  msec_mpi_mean      = msec_mpi_sum     /ncpu;
//    const float  msec_pack_mean     = msec_pack_sum    /ncpu;
//    const float  msec_unpack_mean   = msec_unpack_sum  /ncpu;


    if (comm_.is_rank0()) {
//        std::cout  << "time (total, mpi, pack, unpack) = "
//                   << msec_total_mean << ", " << msec_mpi_mean << ", " << msec_pack_mean << ", " << msec_unpack_mean
//        << " : ( " << msec_total_mean/msec_total_mean*100.0 << ", " << msec_mpi_mean/msec_total_mean*100.0 << ", " << msec_pack_mean/msec_total_mean*100.0 << ", " << msec_unpack_mean/msec_total_mean*100.0 << ") % " << std::endl;
//        std::cout << std::endl;

        std::cout  << "time (total, mpi, pack, unpack) = "
                   << msec_total_max << ", " << msec_mpi_max << ", " << msec_pack_max << ", " << msec_unpack_max
        << " : ( " << msec_total_max/msec_total_max*100.0 << ", " << msec_mpi_max/msec_total_max*100.0 << ", " << msec_pack_max/msec_total_max*100.0 << ", " << msec_unpack_max/msec_total_max*100.0 << " % ) " << std::endl;
        std::cout << std::endl;

        std::cout << "mpi bandwidth( w/o pack/unpack )" << std::endl;
        std::cout << "lv, msec(min, max) = "
                  <<  lv << ", " << msec_mpi_mean << " ( " << msec_mpi_min << ", " << msec_mpi_max << " ) " << std::endl
                  << "mean_bandwidth : GB/s(send), GB/s(recv),  ( GB(send), GB(recv) ) = "
                  << nleaf_s_mean*sizeof(real)/(msec_mpi_mean/1000.0)*1.0e-9 << " , " << nleaf_r_mean*sizeof(real)/(msec_mpi_mean/1000.0)*1.0e-9 << " ( " << nleaf_s_mean*sizeof(real)*1.0e-9 << " , " << nleaf_r_mean*sizeof(real)*1.0e-9 << " ) " << std::endl;
        std::cout << "max_bandwidth  : GB/s(send), GB/s(recv),  ( GB(send), GB(recv) ) = "
                  << nleaf_s_min*sizeof(real)/(msec_mpi_min/1000.0)*1.0e-9 << " , " << nleaf_r_min*sizeof(real)/(msec_mpi_min/1000.0)*1.0e-9 << " ( " << nleaf_s_min*sizeof(real)*1.0e-9 << " , " << nleaf_r_min*sizeof(real)*1.0e-9 << " ) " << std::endl;
        std::cout << "min_bandwidth  : GB/s(send), GB/s(recv),  ( GB(send), GB(recv) ) = "
                  << nleaf_s_max*sizeof(real)/(msec_mpi_max/1000.0)*1.0e-9 << " , " << nleaf_r_max*sizeof(real)/(msec_mpi_max/1000.0)*1.0e-9 << " ( " << nleaf_s_max*sizeof(real)*1.0e-9 << " , " << nleaf_r_max*sizeof(real)*1.0e-9 << " ) " << std::endl;
    }

}
