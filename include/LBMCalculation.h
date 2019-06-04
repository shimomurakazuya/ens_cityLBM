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


class  LBMCalculation {
private:
    MPICommunication    mpiCommunication_;

public:
     LBMCalculation(){}
    ~LBMCalculation(){}

public:
    void LBM_data_assimilation(int t, Field& field);
    void LBM_incompressible_flow(int t, Field& field);

private:
    void LBM_incompressible_flow_w_tb  (int t, Field& field);
    void LBM_incompressible_flow_wo_tb (int t, Field& field);

    void LBM_stream_collision_opt( const int lv, Field& field,
                                   const int*  id_tasks_wo_bc, const int  num_tasks_wo_bc,
                                   const int*  id_tasks_w_bc,  const int  num_tasks_w_bc);

    void LBM_to_euler_variables(const int lv, Field& field, const int*  id_tasks, const int num_tasks);

    void ScalarAdvection(const int lv, Field& field, const int*  id_tasks, const int  num_tasks);

    void LBM_L2F (const int lv, Field& field, const int*  id_tasks, const int num_tasks);
    void LBM_L2FA(const int lv, Field& field, const int*  id_tasks, const int num_tasks);


    void Val_L2F (const int lv, Field& field, const int*  id_tasks, const int num_tasks);
    void Val_L2FA(const int lv, Field& field, const int*  id_tasks, const int num_tasks);

    void LBM_F2L   (const int lv, Field& field, const int*  id_tasks, const int num_tasks);

    void Val_F2L   (const int lv, Field& field, const int*  id_tasks, const int num_tasks);

    void boundary_conditions(const int lv, Field& field, const int*  id_tasks, const int  num_tasks);

    void ValuesCommunications(const int lv, Field& field);

    void set_cuda_block(int block3[], int blockDimX, int blockDimY, int blockDimZ);
};


#endif
