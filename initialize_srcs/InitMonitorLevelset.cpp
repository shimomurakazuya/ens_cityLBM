#include "InitMonitorLevelset.h"
#include "FuncObj.h"
#include "FuncLoop.h"
#include "defineCal.h"
#include "mpi_wrapper.hpp"
#include "range.hpp"


std::vector<int> InitMonitorLevelset::
create_id_flags(
    const Grid* grids,
    const MapData& map
    )
{
    const auto& mpi_world = util::mpi(MPI_COMM_WORLD);
    if (mpi_world.rank() == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    // refinement condition //
    // return cal_flag //

    constexpr int  lv_max = DefAMR::LV_MAX;

    std::array<int, lv_max> nx, ny, nz;
    std::array<real, lv_max> offset_x, offset_y, offset_z, dx;

    FuncLoop::Loop1d  loop1d(0, lv_max);
    loop1d.for_each( [&nx, grids](int i){ nx.at(i) = grids[i].nx().x(); } );
    loop1d.for_each( [&ny, grids](int i){ ny.at(i) = grids[i].nx().y(); } );
    loop1d.for_each( [&nz, grids](int i){ nz.at(i) = grids[i].nx().z(); } );

    loop1d.for_each( [&offset_x, grids](int i){ offset_x.at(i) = grids[i].offset().x(); } );
    loop1d.for_each( [&offset_y, grids](int i){ offset_y.at(i) = grids[i].offset().y(); } );
    loop1d.for_each( [&offset_z, grids](int i){ offset_z.at(i) = grids[i].offset().z(); } );

    loop1d.for_each( [&dx, grids](int i){ dx.at(i) = grids[i].dx(); } );

    // lv //
    int  nx_fine = nx.at(lv_max-1);
    int  ny_fine = ny.at(lv_max-1);
    int  nz_fine = nz.at(lv_max-1);
    std::vector<int> amr_lv(nx_fine * ny_fine * nz_fine);

    init_amr_level(
        amr_lv,
        offset_x[lv_max-1] + 0.5*dx[lv_max-1],
        offset_y[lv_max-1] + 0.5*dx[lv_max-1],
        offset_z[lv_max-1] + 0.5*dx[lv_max-1],
        dx[lv_max-1],
        nx[lv_max-1],
        ny[lv_max-1],
        nz[lv_max-1],
        map
        );

    check_amr_level(
        amr_lv,
        nx[lv_max-1],
        ny[lv_max-1],
        nz[lv_max-1]
        );


    std::vector<int>  cal_flags;
    std::array<int, lv_max> count_cal;

    for (int lv=0; lv<lv_max; lv++) {
        count_cal.at(lv) = 0;
        for (int k=0; k<nz[lv]; k++) {
        for (int j=0; j<ny[lv]; j++) {
        for (int i=0; i<nx[lv]; i++) {
            const int dlv = (lv_max-1) - lv;
            const int id = index_amr(dlv, i,j,k, nx_fine, ny_fine, nz_fine);
            int  tmp_amr_lv = amr_lv[id];
            cal_flags.push_back( (tmp_amr_lv == lv) ); // true, false //
            if(tmp_amr_lv == lv) { ++count_cal.at(lv); }
        }
        }
        }
    }

    mpi_world.for_each_rank([&]() { 
        std::cout << "number of cal_true = " << std::flush;
        for(int lv=0; lv<lv_max; lv++) {
            const auto& count = count_cal.at(lv);
            const auto& cmax = nx.at(lv) * ny.at(lv) * nz.at(lv);
            std::cout << count << "/" << cmax << " " << std::flush;
        }
        std::cout << std::endl;
    });

    return cal_flags;
}

void InitMonitorLevelset::
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
    )
{
#pragma omp parallel for collapse(3)
    for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
        const int id = index(i,j,k, nx,ny,nz);

        const real xx = offset_x + dx*i;
        const real yy = offset_y + dx*j;
        const real zz = offset_z + dx*k;

        amr_lv.at(id) = func_monitor_levelset_to_amr_level_mapfile(xx, yy, zz, dx, map);
    }
    }
    }
}


void InitMonitorLevelset::
check_amr_level(
    std::vector<int>& amr_lv,
    const int   nx,
    const int   ny,
    const int   nz
    )
{
    const auto mpi_world = util::mpi(MPI_COMM_WORLD);
    if (mpi_world.rank() == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    constexpr int  lv_max = DefAMR::LV_MAX;

    // block //
    for (int _lv=0; _lv<lv_max; _lv++) {
#pragma omp parallel for collapse(3)
        for (int k=0; k<nz; k++) {
        for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {

            const int  id = index(i,j,k, nx,ny,nz);

            const int  tmp_amr_lv = amr_lv.at(id);
            if(tmp_amr_lv < 0) { continue; } // memreduce

            const int  tmp_dlv = (lv_max-1) - tmp_amr_lv;
            const int  tmp_nx = pow(2, tmp_dlv) * 2;

            const int  ist = ((int)(i/tmp_nx))*tmp_nx;
            const int  jst = ((int)(j/tmp_nx))*tmp_nx;
            const int  kst = ((int)(k/tmp_nx))*tmp_nx;

            for (int kk=0; kk<tmp_nx; kk++) {
            for (int jj=0; jj<tmp_nx; jj++) {
            for (int ii=0; ii<tmp_nx; ii++) {
                const int idn = index(ist+ii, jst+jj, kst+kk, nx,ny,nz);

                if(amr_lv.at(idn) < 0) { continue; } // memreduce
                const int tmp_amr_lv_max = std::max( tmp_amr_lv, amr_lv.at(idn) );
                amr_lv.at(idn) = tmp_amr_lv_max;

            }
            }
            }

        }
        }
        }

    }

    // neighbor //
#pragma omp parallel for collapse(3)
    for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
        const int id = index(i,j,k, nx,ny,nz);

        const int  tmp_amr_lv = amr_lv[id];
        if(tmp_amr_lv < 0) { continue; } // memreduce

        for (int kk=-1; kk<=1; kk++) {
        for (int jj=-1; jj<=1; jj++) {
        for (int ii=-1; ii<=1; ii++) {
            int in = i + ii;
            int jn = j + jj;
            int kn = k + kk;

            if ( !OuterBoundaryConditions::OuterPeriodicCommX ) {
//                continue;
                if ( in < 0 ) { in = 0; }   if ( in > nx-1 ) { in = nx-1; }
            }
            if ( !OuterBoundaryConditions::OuterPeriodicCommY ) {
//                continue;
                if ( jn < 0 ) { jn = 0; }   if ( jn > ny-1 ) { jn = ny-1; }
            }
            if ( !OuterBoundaryConditions::OuterPeriodicCommZ ) {
//                continue;
                if ( kn < 0 ) { kn = 0; }   if ( kn > nz-1 ) { kn = nz-1; }
            }

            const int idn = index(in, jn, kn, nx,ny,nz);
            const int tmp_amr_lvn = amr_lv[idn];
            if(tmp_amr_lvn < 0) { continue; } // memreduce

            if ( std::fabs(tmp_amr_lv - tmp_amr_lvn) >= 2 ) {
                std::cout << __PRETTY_FUNCTION__ << " : " << "error neighbor lv offset is larger than 2" << std::endl;
                exit(-1);
            }
        }
        }
        }
    }
    }
    }
}

int InitMonitorLevelset::
func_monitor_levelset_to_amr_level_mapfile(
    const real x,
    const real y,
    const real z,
    const real dx,
    const MapData& map
    )
{
#if defined(VVTARGET_TestOklahoma)
    constexpr int  lv_max = DefAMR::LV_MAX;

//    constexpr real height0 = 135.0 + 4.0 - 0.1;
//    constexpr real height1 = 768.0 + 4.0 - 0.1;

//    constexpr real height0 =  150.0 + 4.0 - 0.1;
//    constexpr real height1 =  700.0 + 4.0 - 0.1;

#if EXPECT_DX_MESH == 1
#ifdef REGRESSION_TEST
#error regression test should be run with EXPECT_DX_MESH=4
#endif
    // oklahoma 1m //
    constexpr real tmp_dx  =   1.0;
    constexpr real height0 =  64.0    -   4.0 - tmp_dx*3;
    constexpr real height1 =  height0 + 128.0 - tmp_dx*2*3;
#elif EXPECT_DX_MESH == 2
#ifdef REGRESSION_TEST
#error regression test should be run with EXPECT_DX_MESH=4
#endif
    // oklahoma 2m //
    constexpr real tmp_dx  =   2.0;
    constexpr real height0 =  96.0    -   8.0 - tmp_dx*3;
//    constexpr real height0 =  128.0    -   8.0 - tmp_dx*3;
    constexpr real height1 =  height0 + 256.0 - tmp_dx*2*3;
#elif EXPECT_DX_MESH == 4 || !defined(EXPECT_DX_MESH)
    // oklahoma 4m //
    constexpr real tmp_dx  =    4.0;
    constexpr real height0 =  256.0   -  16.0 - tmp_dx*3;
    constexpr real height1 =  height0 + 512.0 - tmp_dx*2*3;
#else
#error no matching for defined EXPECT_DX_MESH
#endif
    // citybox //
    constexpr real levelset_threash_fine = -99999;
    constexpr real box_radius = 99999; // meter

    // lv_max = 3 //
    if      ( fabs(z) <= height0 ) {
        #if defined(REGRESSION_TEST) || defined(NO_CITYBOX)
        return lv_max-1;
        #endif

        // citybox //
        if ( std::fabs(x) > box_radius || std::fabs(y) > box_radius) {
            return lv_max -2;
        }
        // amr with map //
        bool near_or_in_obj = false, has_fluid = false;
        constexpr int nx = DefAMR::NX_LEAF;
        const real dx_mesh = dx / nx;
        for(int k=-1; k<nx+1; k++) for(int j=-1; j<nx+1; j++) for(int i=-1; i<nx+1; i++) {
            const auto& xx = x + dx_mesh * i;
            const auto& yy = y + dx_mesh * j;
            const auto& zz = z + dx_mesh * k;
            const auto& levelset = map.get_levelset(xx, yy, zz, lv_max-1, dx_mesh);
            has_fluid |= (levelset <= 0);
            near_or_in_obj |= (levelset >= levelset_threash_fine);
        }
        if(near_or_in_obj) {
            if(has_fluid) { return lv_max -1; }
            return -404; // memreduce: if negative: no memory for this leaf
        } 
        return lv_max -1; // citybox: always finest where z <= height0
    }
    else if ( fabs(z) <= height1 ) {
        return lv_max -2;
    }
    else {
        return lv_max -3;
    }

//    if ( (x < 120 || x > 1800) && amr_lv == lv_max-1 ) {
//        amr_lv = lv_max - 2;
//    }
//
//    if ( (x < 60 || x > 1860) && amr_lv == lv_max-2 ) {
//        amr_lv = lv_max - 3;
//    }
#elif defined(VVTARGET_TestFlowAroundCube)
    constexpr int lv_max = DefAMR::LV_MAX;
    constexpr real H = 0.1;
    const real xx = std::fabs(x) / H;
    const real yy = std::fabs(y) / H;
    const real zz = std::fabs(z) / H;
    int amr_lv = -1;
    if       (lv_max == 3) {
        if     (xx < 20. && yy <  6. && zz <  8. ) { amr_lv = lv_max - 1; }
        else if(xx < 40. && yy < 12. && zz < 16. ) { amr_lv = lv_max - 2; }
        else { amr_lv = -1; }
    } else if(lv_max == 4) {
        if     (xx <  3. && yy <  2. && zz <  3. ) { amr_lv = lv_max - 1; }
        else if(xx < 20. && yy <  6. && zz <  8. ) { amr_lv = lv_max - 2; }
        else if(xx < 40. && yy < 12. && zz < 16. ) { amr_lv = lv_max - 3; }
        else { amr_lv = -1; }
    } else if(lv_max == 5) {
        if     (xx <  2. && yy <  1. && zz <  1.5) { amr_lv = lv_max - 1; }
        else if(xx <  4. && yy <  3. && zz <  4. ) { amr_lv = lv_max - 2; }
        else if(xx < 12. && yy <  6. && zz <  8. ) { amr_lv = lv_max - 3; }
        else if(xx < 30. && yy < 10. && zz < 10. ) { amr_lv = lv_max - 4; }
        else { amr_lv = 0; }
    }
    return std::max(0, amr_lv);
#elif  defined(VVTARGET_TestCavityFlow)
    constexpr int lv_max = DefAMR::LV_MAX;
    constexpr int nx_leaf = DefAMR::NX_LEAF;
    constexpr real dx_fine_leaf = 1./256. * nx_leaf;
    static_assert(lv_max == 2, "invalid configuration of AMR level");
    const real xx = x;
    const real yy = y;
    if(std::fabs(xx) >= 0.5 - 3. * dx_fine_leaf || std::fabs(yy) >= 0.5 - 3. * dx_fine_leaf) { return lv_max-1; }
    else { return 0; }
#endif
}


