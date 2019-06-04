#include "TimeEvolutionLoop.h"


void  TimeEvolutionLoop::
time_evolution(
    Field&          field,
    WorkerThread&   workerThread
    )
{
    for (int t=start_step_; t<finish_step_; t++) {
        // Timer //
        field.allTimeInfo().StartTimer("AllCal");
        {
            if (rank_ == 0) { std::cout << "-" << std::flush; }

            step_ = t;


            // output //
            field.allTimeInfo().StartTimer("Output");
            {
                outputFunc_.OutputFluidData(t, field, workerThread);
            }
            field.allTimeInfo().StopTimer("Output");
            field.allTimeInfo().SubmitElapsedTimeInfo("Output", 0, field.taskID().num_task_lbm_per_step(), field.parameters().step_now());


            // data assimilation //
            field.allTimeInfo().StartTimer("DataAssimilation");
            {
                lbmCalculation_.LBM_data_assimilation(t, field);
            }
            field.allTimeInfo().StopTimer("DataAssimilation");
            field.allTimeInfo().SubmitElapsedTimeInfo("DataAssimilation", 0, field.taskID().num_task_lbm_per_step(), field.parameters().step_now());


            // calculation //
            field.allTimeInfo().StartTimer("LBM_Calculation");
            {
                lbmCalculation_.LBM_incompressible_flow(t, field);
            }
            field.allTimeInfo().StopTimer("LBM_Calculation");
            field.allTimeInfo().SubmitElapsedTimeInfo("LBM_Calculation", 0, field.taskID().num_task_lbm_per_step(), field.parameters().step_now());


            MPI_Barrier(MPI_COMM_WORLD);
            time_ += dt_;
        }
        field.allTimeInfo().StopTimer("AllCal");
        field.allTimeInfo().SubmitElapsedTimeInfo("AllCal", 0, field.taskID().num_task_lbm_per_step(), field.parameters().step_now());
    }
}
