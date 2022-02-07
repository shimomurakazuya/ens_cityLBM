#pragma once
#ifndef LBMCALCULATION_H_
#define LBMCALCULATION_H_


#include "defineCUDA.h"
#include <iostream>
#include <string>
#include "defineCal.h"
#include "definePrecision.h"
#include "defineLBM.h"
#include "Field.h"
#include "MPICommunication.h"
#include "MPICommunicationMG.h"
#include "MPICommEnsemble.h"


class  LBMCalculation {
private:
    const MPICommEnsemble comm_;
    MPICommunicationMG  mpiCommunicationMG_;
    MPICommunication    mpiCommunication_;

public:
    LBMCalculation() = delete;

    LBMCalculation(const MPICommEnsemble comm): 
    comm_(comm),
    mpiCommunicationMG_(comm),
    mpiCommunication_(comm)
    {}

    ~LBMCalculation(){}

public:
    void LBM_data_assimilation(int t, Field& field);
    void LBM_incompressible_flow(int t, Field& field);

public: // trick for __device__ lambda
    void NSUpdate(const int lv, Field& field, const int*  id_tasks, const int  num_tasks);

//private:
    void LBM_incompressible_flow_w_tb  (int t, Field& field);
    void LBM_incompressible_flow_wo_tb (int t, Field& field);

    void LBM_stream_collision_opt( const int lv, Field& field,
                                   const int*  id_tasks_wo_bc, const int  num_tasks_wo_bc,
                                   const int*  id_tasks_w_bc,  const int  num_tasks_w_bc);

    void LBM_to_euler_variables(const int lv, Field& field, const int*  id_tasks, const int num_tasks);
    void LBM_to_euler_variables_with_time_average(const int lv, Field& field, const int*  id_tasks, const int num_tasks);


    void LBM_L2F (const int lv, Field& field, const int*  id_tasks, const int num_tasks);
    void LBM_L2FA(const int lv, Field& field, const int*  id_tasks, const int num_tasks);


    void Val_L2F (const int lv, Field& field, const int*  id_tasks, const int num_tasks);
    void Val_L2FA(const int lv, Field& field, const int*  id_tasks, const int num_tasks);

    void LBM_F2L   (const int lv, Field& field, const int*  id_tasks, const int num_tasks);

    void Val_F2L   (const int lv, Field& field, const int*  id_tasks, const int num_tasks);

    void boundary_conditions(const int lv, Field& field, const int*  id_tasks, const int  num_tasks);

    void NSCommunications (const int lv, Field& field);
    void TCommunications (const int lv, Field& field);
    void ScalarCommunications (const int lv, Field& field);
    void TimeAverageCommunications (const int lv, Field& field);
    void LBMCommunications(const int lv, Field& field);

    void set_cuda_block(int block3[], int blockDimX, int blockDimY, int blockDimZ);
};


#endif
