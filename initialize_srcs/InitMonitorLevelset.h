#pragma once
#ifndef INITMONITORLEVELSET_H_
#define INITMONITORLEVELSET_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <mpi.h>
#include "defineAMR.h"
#include "definePrecision.h"
#include "Grid.h"
#include "FuncMapData.h"


class  InitMonitorLevelset {
private:
    enum  MonitorLevelset {
        monitor_map,
        monitor_cavity2d,
        monitor_nc2d,
        monitor_nc3d,
        monitor_flow_cube,
        monitor_other
    };

    MonitorLevelset   monitorLevelset_;

public:
    InitMonitorLevelset () :
        monitorLevelset_(monitor_map)
//        monitorLevelset_(monitor_cavity2d)
//        monitorLevelset_(monitor_nc2d)
//        monitorLevelset_(monitor_nc3d)
//        monitorLevelset_(monitor_flow_cube)
    {
    }

    ~InitMonitorLevelset () {}

public:
    std::vector<int> create_id_flags(const Grid* grids, const MapData& map);

//    std::vector<int> create_id_flags_cavity2d(const Grid* grids);
//    std::vector<int> create_id_flags_flow_around_cube(const Grid* grids);


private:
    void
    init_amr_level(
        std::vector<int>& amr_lv,
        const real  offset_x,
        const real  offset_y,
        const real  offset_z,
        const real  dx,
        const int   nx,
        const int   ny,
        const int   nz,
        const MapData& map
        );


    void
    check_amr_level(
        std::vector<int>&  amr_lv,
        const int   nx,
        const int   ny,
        const int   nz
        );


    int func_monitor_levelset_to_amr_level_mapfile(
        const real x,
        const real y,
        const real z,
        const real dx_fine,
        const MapData& map
    );


    int  index_amr(int dlv, int i, int j, int k, int nx, int ny, int nz)
    {
        int  tmp = pow(2, dlv);
        return  index(tmp*i, tmp*j, tmp*k, nx, ny, nz);
    }

    int  index(int i, int j, int k, int nx, int ny, int nz)
    {
        const int  ii = (i+nx)%nx;
        const int  jj = (j+ny)%ny;
        const int  kk = (k+nz)%nz;

        return  ii + nx*jj + nx*ny*kk;
    }

};


#endif
