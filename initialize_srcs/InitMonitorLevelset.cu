#include "InitMonitorLevelset.h"
#include "FuncObj.h"
#include "FuncLoop.h"
#include "defineCal.h"


std::vector<int> InitMonitorLevelset::
create_id_flags(
    const Grid* grids
    )
{
    int  rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    // refinement condition //
    // header : monitorLevelset_ //

    // return cal_flag //

    constexpr int  lv_max = DefAMR::LV_MAX;

    int  nx[lv_max];
    int  ny[lv_max];
    int  nz[lv_max];

    real offset_x[lv_max];
    real offset_y[lv_max];
    real offset_z[lv_max];

    real dx[lv_max];

    FuncLoop::Loop1d  loop1d(0, lv_max);
    loop1d.for_each( [&nx, grids](int i){ nx[i] = grids[i].nx().x(); } );
    loop1d.for_each( [&ny, grids](int i){ ny[i] = grids[i].nx().y(); } );
    loop1d.for_each( [&nz, grids](int i){ nz[i] = grids[i].nx().z(); } );

    loop1d.for_each( [&offset_x, grids](int i){ offset_x[i] = grids[i].offset().x(); } );
    loop1d.for_each( [&offset_y, grids](int i){ offset_y[i] = grids[i].offset().y(); } );
    loop1d.for_each( [&offset_z, grids](int i){ offset_z[i] = grids[i].offset().z(); } );

    loop1d.for_each( [&dx, grids](int i){ dx[i] = grids[i].dx(); } );


    // lv //
    int  nx_fine = nx[lv_max-1];
    int  ny_fine = ny[lv_max-1];
    int  nz_fine = nz[lv_max-1];
    int* amr_lv = new int [nx_fine*ny_fine*nz_fine];


    init_amr_level(
        amr_lv,
        offset_x[lv_max-1] + 0.5*dx[lv_max-1],
        offset_y[lv_max-1] + 0.5*dx[lv_max-1],
        offset_z[lv_max-1] + 0.5*dx[lv_max-1],
        dx[lv_max-1],
        nx[lv_max-1],
        ny[lv_max-1],
        nz[lv_max-1]
        );

    check_amr_level(
        amr_lv,
        nx[lv_max-1],
        ny[lv_max-1],
        nz[lv_max-1]
        );


//    std::cout << "initialize cal_flags\n";

    std::vector<int>  cal_flags;

//    int  nsum = 0;
//    for (int lv=0; lv<lv_max; lv++) { nsum += nx[lv]*ny[lv]*nz[lv]; }
//    cal_flags.reserve(nsum);

    int  id_tmp = 0;
    for (int lv=0; lv<lv_max; lv++) {
        for (int k=0; k<nz[lv]; k++) {
        for (int j=0; j<ny[lv]; j++) {
        for (int i=0; i<nx[lv]; i++) {
            const int dlv = (lv_max-1) - lv;
            const int id = index_amr(dlv, i,j,k, nx_fine, ny_fine, nz_fine);

            int  tmp_amr_lv = amr_lv[id];

            cal_flags.push_back( (tmp_amr_lv == lv) ); // true, false //
//            cal_flags[id_tmp] = (tmp_amr_lv == lv); // true, false //

            id_tmp++;
        }
        }
        }
    }

    // check //
#if 0
    id_tmp = 0;
    for (int lv=0; lv<lv_max; lv++) {
        std::cout << "lv = " << lv << std::endl;
        for (int k=0; k<nz[lv]; k++) {
        for (int j=0; j<ny[lv]; j++) {
        for (int i=0; i<nx[lv]; i++) {
            std::cout << cal_flags[id_tmp] << " ";

            id_tmp++;
        }
        std::cout << std::endl;
        }
        std::cout << std::endl;
        }
    }
    exit(-1);
#endif
    // check //

    delete [] amr_lv;

    return cal_flags;
}


void InitMonitorLevelset::
init_amr_level(
          int*  amr_lv,
    const real  offset_x,
    const real  offset_y,
    const real  offset_z,
    const real  dx,
    const int   nx,
    const int   ny,
    const int   nz
    )
{
    for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
        const int id = index(i,j,k, nx,ny,nz);

        const real xx = offset_x + dx*i;
        const real yy = offset_y + dx*j;
        const real zz = offset_z + dx*k;

        real tmp_levelset;
        int  tmp_amr_lv;
        if      (monitorLevelset_ == monitor_cavity2d) {
            tmp_levelset = func_monitor_levelset_cavity2d(xx, yy, zz);
            tmp_amr_lv   = func_monitor_levelset_to_amr_level_cavity2d( fabs(tmp_levelset), dx);
        }
        else if (monitorLevelset_ == monitor_nc2d) {
            tmp_levelset = func_monitor_levelset_nc2d(xx, yy, zz);
            tmp_amr_lv   = func_monitor_levelset_to_amr_level_nc2d( fabs(tmp_levelset), dx);
        }
        else if (monitorLevelset_ == monitor_nc3d) {
            tmp_levelset = func_monitor_levelset_nc3d(xx, yy, zz);
            tmp_amr_lv   = func_monitor_levelset_to_amr_level_nc3d( fabs(tmp_levelset), dx);
        }
        else if (monitorLevelset_ == monitor_map) {
            tmp_amr_lv = func_amr_level_flow_map(xx, yy, zz, dx);
        }
        else if (monitorLevelset_ == monitor_flow_cube) {
//            tmp_amr_lv = func_amr_level_flow_cube(xx, yy, zz, dx);
            tmp_amr_lv = func_amr_level_flow_map(xx, yy, zz, dx);
//            tmp_amr_lv = func_amr_level_channel_flow(xx, yy, zz, dx);
        }

        amr_lv[id] = tmp_amr_lv;
    }
    }
    }

    // check //
#if 0
    for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
        const int id = index(i,j,k, nx,ny,nz);

        const real xx = offset_x + dx*i;
        const real yy = offset_y + dx*j;
        const real zz = offset_z + dx*k;

        const real levelset = func_monitor_levelset(xx, yy, zz);
        std::cout << amr_lv[id] << " ";
    }
        std::cout << std::endl;
    }
        std::cout << std::endl;
    }
#endif
}


void InitMonitorLevelset::
check_amr_level(
          int*  amr_lv,
    const int   nx,
    const int   ny,
    const int   nz
    )
{
    int  rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    constexpr int  lv_max = DefAMR::LV_MAX;

    bool* check_list = new bool[nx*ny*nz];
#pragma omp parallel for collapse(3)
    for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
        const int  id = index(i,j,k, nx,ny,nz);
        check_list[id] = false;
    }
    }
    }


    // block //
    for (int _lv=0; _lv<lv_max; _lv++) {
        const int lv = lv_max-1 - _lv;

#pragma omp parallel for collapse(3) shared(check_list)
        for (int k=0; k<nz; k++) {
        for (int j=0; j<ny; j++) {
        for (int i=0; i<nx; i++) {

            const int  id = index(i,j,k, nx,ny,nz);
            if ( check_list[id] ) { continue; }

            const int  tmp_amr_lv = amr_lv[id];
            if (tmp_amr_lv != lv) { continue; }

            const int  tmp_dlv = (lv_max-1) - tmp_amr_lv;
            const int  tmp_nx = pow(2, tmp_dlv) * 2;

            const int  ist = ((int)(i/tmp_nx))*tmp_nx;
            const int  jst = ((int)(j/tmp_nx))*tmp_nx;
            const int  kst = ((int)(k/tmp_nx))*tmp_nx;

            for (int kk=0; kk<tmp_nx; kk++) {
            for (int jj=0; jj<tmp_nx; jj++) {
            for (int ii=0; ii<tmp_nx; ii++) {
                const int idn = index(ist+ii, jst+jj, kst+kk, nx,ny,nz);

                const int tmp_amr_lv_max = std::max( tmp_amr_lv, amr_lv[idn] );
                amr_lv[idn] = tmp_amr_lv_max;

                check_list[idn] = true;
            }
            }
            }

        }
        }
        }

    }

    delete [] check_list;

    // neighbor //
#pragma omp parallel for collapse(3)
    for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
        const int id = index(i,j,k, nx,ny,nz);

        const int  tmp_amr_lv = amr_lv[id];

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


void InitMonitorLevelset::
_check_amr_level(
          int*  amr_lv,
    const int   nx,
    const int   ny,
    const int   nz
    )
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;

    constexpr int  lv_max = DefAMR::LV_MAX;

    bool* check_list = new bool[nx*ny*nz];
#pragma omp parallel for collapse(3)
    for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
        const int  id = index(i,j,k, nx,ny,nz);
        check_list[id] = false;
    }
    }
    }


//    auto  cast_n = [](int i, int n){ return ((int)(i/n))*n; };

    // block //
//    std::cout << "check myself\n";
//#pragma omp parallel for collapse(3)
#pragma omp parallel for collapse(3) shared(check_list)
    for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
//        if (i==0 && j==0 && k%10==0) { std::cout << "k = " << k << " / " << nz << std::endl; }

        const int  id = index(i,j,k, nx,ny,nz);
//        if ( check_list[id] ) { continue; }

        const int  tmp_amr_lv = amr_lv[id];
        const int  tmp_dlv = (lv_max-1) - tmp_amr_lv;

//        const int  tmp_nx = pow(2.0, tmp_dlv) * 2;
        const int  tmp_nx = pow(2, tmp_dlv) * 2;

//        const int  ist = cast_n(i, tmp_nx);
//        const int  jst = cast_n(j, tmp_nx);
//        const int  kst = cast_n(k, tmp_nx);
        const int  ist = ((int)(i/tmp_nx))*tmp_nx;
        const int  jst = ((int)(j/tmp_nx))*tmp_nx;
        const int  kst = ((int)(k/tmp_nx))*tmp_nx;

        for (int kk=0; kk<tmp_nx; kk++) {
        for (int jj=0; jj<tmp_nx; jj++) {
        for (int ii=0; ii<tmp_nx; ii++) {
            const int idn = index(ist+ii, jst+jj, kst+kk, nx,ny,nz);

            const int tmp_amr_lv_max = std::max( tmp_amr_lv, amr_lv[idn] );
            amr_lv[idn] = tmp_amr_lv_max;

            check_list[idn] = true;

//            if (tmp_amr_lv_max == lv_max-1) { ii = tmp_nx; jj = tmp_nx; kk = tmp_nx; }
        }
        }
        }

    }
    }
    }

    delete [] check_list;

    // neighbor //
#pragma omp parallel for collapse(3)
    for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
    for (int i=0; i<nx; i++) {
        const int id = index(i,j,k, nx,ny,nz);

        const int  tmp_amr_lv = amr_lv[id];

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


real InitMonitorLevelset::
func_monitor_levelset_cavity2d(
    const real x,
    const real y,
    const real z
    )
{
    const real  xyz[] = { x, y, z };

    real  val;
    if ( monitorLevelset_ == monitor_cavity2d ) {
        constexpr real  box_min[] = {  -0.5,  -0.5,  -10.0 };
        constexpr real  box_max[] = {   0.5,   0.5,   10.0 };

        val = -FuncObj::length_from_box(xyz, box_min, box_max);
    }
    else {
        val = FuncObj::lv_fluid();
    }

    return  val;
}


real InitMonitorLevelset::
func_monitor_levelset_nc2d(
    const real x,
    const real y,
    const real z
    )
{
    const real  xyz[] = { x, y, z };

    real  val;
    if ( monitorLevelset_ == monitor_nc2d ) {
        constexpr real  box_min[] = {  -0.5,  -10.0,  -10.0 };
        constexpr real  box_max[] = {   0.5,   10.0,   10.0 };
//        constexpr real  box_min[] = {  -10.0,  -0.5,  -10.0 };
//        constexpr real  box_max[] = {   10.0,   0.5,   10.0 };

        val = -FuncObj::length_from_box(xyz, box_min, box_max);
    }
    else {
        val = FuncObj::lv_fluid();
    }

    return  val;
}


real InitMonitorLevelset::
func_monitor_levelset_nc3d(
    const real x,
    const real y,
    const real z
    )
{
    const real  xyz[] = { x, y, z };

    real  val;
    if ( monitorLevelset_ == monitor_nc3d ) {
        constexpr real  box_min[] = {  -0.3,  -0.3,   0.0 };
        constexpr real  box_max[] = {   0.3,   0.3,   0.8 };

        val = -FuncObj::length_from_box(xyz, box_min, box_max);
    }
    else {
        val = FuncObj::lv_fluid();
    }

    return  val;
}


int InitMonitorLevelset::
func_monitor_levelset_to_amr_level_cavity2d(
    const real  levelset,
    const real  dx_fine
    )
{
    constexpr int  lv_max = DefAMR::LV_MAX;

    real  ls_th[lv_max+1];
    ls_th[0] = 0.0;

    // max lv 5 //
//    const real dx_lv[] = { dx_fine, dx_fine*2, dx_fine, dx_fine, dx_fine };

    for (int i=1; i<=lv_max; i++) {
        // hand made parameter //
//        ls_th[i] = ls_th[i-1] + 0.5*dx_fine*pow(2.0, i);
        ls_th[i] = ls_th[i-1] + 1.0*dx_fine*pow(2.0, i);
//        ls_th[i] = ls_th[i-1] + 2.0*dx_fine*pow(2.0, i);
//        ls_th[i] = ls_th[i-1] + 2.5*dx_fine*pow(2.0, i);
    }


    int  tmp = lv_max-1;
    for (int i=0; i<lv_max; i++) {
        if ( levelset >= ls_th[i] && levelset < ls_th[i+1] ) {
            tmp = i;
            break;
        }
    }
    int  amr_lv = (lv_max-1) - tmp;


    return  amr_lv;
}


int InitMonitorLevelset::
func_monitor_levelset_to_amr_level_nc2d(
    const real  levelset,
    const real  dx_fine
    )
{
    constexpr int  lv_max = DefAMR::LV_MAX;

    real  ls_th[lv_max+1];
    ls_th[0] = 0.0;

    // max lv 5 //
//    const real dx_lv[] = { dx_fine, dx_fine*2, dx_fine, dx_fine, dx_fine };
    const real dx_lv[] = { 0.15, 0.15, 0.1, 0.1 };

    for (int i=1; i<=lv_max; i++) {
        // hand made parameter //
//        ls_th[i] = ls_th[i-1] + 0.5*dx_fine*pow(2.0, i);
//        ls_th[i] = ls_th[i-1] + 1.0*dx_fine*pow(2.0, i);
//        ls_th[i] = ls_th[i-1] + 2.0*dx_fine*pow(2.0, i);
        ls_th[i] = ls_th[i-1] + 5.0*dx_fine*pow(2.0, i);

        ls_th[i] = ls_th[i-1] + dx_lv[i-1];
    }


    int  tmp = lv_max-1;
    for (int i=0; i<lv_max; i++) {
        if ( levelset >= ls_th[i] && levelset < ls_th[i+1] ) {
            tmp = i;
            break;
        }
    }
    int  amr_lv = (lv_max-1) - tmp;


    return  amr_lv;
}


int InitMonitorLevelset::
func_monitor_levelset_to_amr_level_nc3d(
    const real  levelset,
    const real  dx_fine
    )
{
    constexpr int  lv_max = DefAMR::LV_MAX;

    real  ls_th[lv_max+1];
    ls_th[0] = 0.0;

//    const real dx_lv[] = { 0.05, 0.05, 0.05, 0.05 };
    const real dx_lv[] = { 0.1, 0.1, 0.1, 0.1 };

    for (int i=1; i<=lv_max; i++) {
        // hand made parameter //
//        ls_th[i] = ls_th[i-1] + 0.5*dx_fine*pow(2.0, i);
//        ls_th[i] = ls_th[i-1] + 1.0*dx_fine*pow(2.0, i);
//        ls_th[i] = ls_th[i-1] + 2.0*dx_fine*pow(2.0, i);
//        ls_th[i] = ls_th[i-1] + 5.0*dx_fine*pow(2.0, i);

        ls_th[i] = ls_th[i-1] + dx_lv[i-1];
    }


    int  tmp = lv_max-1;
    for (int i=0; i<lv_max; i++) {
        if ( levelset >= ls_th[i] && levelset < ls_th[i+1] ) {
            tmp = i;
            break;
        }
    }
    int  amr_lv = (lv_max-1) - tmp;


    return  amr_lv;
}


real InitMonitorLevelset::
func_amr_level_flow_cube(
    const real x,
    const real y,
    const real z,
    const real dx
    )
{
    constexpr int  lv_max = DefAMR::LV_MAX;
    constexpr real H = 0.1;

    int amr_lv;
#if 1
    // lv_max = 3 //
    if      ( fabs(x) < 20*H && fabs(y) <= 6*H && fabs(z) <= 8*H ) {
//    if      ( fabs(x) < 15*H && fabs(z) <= 7.5*H ) {
//    if      ( fabs(z) <= 7.5*H ) {
        amr_lv = lv_max - 1;
    }
    else if ( fabs(x) < 40*H && fabs(y) <= 12*H && fabs(z) <= 16*H ) {
//    else if ( fabs(x) < 40*H && fabs(z) <= 15*H ) {
//    else if ( fabs(z) <= 15*H ) {
        amr_lv = lv_max - 2;
    }
    else {
        amr_lv = lv_max - 3;
    }
#endif

#if 0
    // lv_max = 4 //
    if      ( fabs(x) < 3.0*H && fabs(y) <= 2.0*H && fabs(z) <= 3.0*H ) {
        amr_lv = lv_max - 1;
    }
    else if ( fabs(x) < 20.0*H && fabs(y) <= 6.0*H && fabs(z) <= 8.0*H ) {
        amr_lv = lv_max - 2;
    }
    else if ( fabs(x) < 40.0*H && fabs(y) <= 12.0*H && fabs(z) <= 16.0*H ) {
        amr_lv = lv_max - 3;
    }
    else {
        amr_lv = lv_max - 4;
    }
#endif

#if 0
    // lv_max = 5 //
    if      ( fabs(x) < 2.0*H && fabs(y) <= 1.0*H && fabs(z) <= 1.5*H ) {
        amr_lv = lv_max - 1;
    }
    else if ( fabs(x) < 4.0*H && fabs(y) <= 3.0*H && fabs(z) <= 4.0*H ) {
        amr_lv = lv_max - 2;
    }
    else if      ( fabs(x) < 12.0*H && fabs(y) <= 6.0*H && fabs(z) <= 8.0*H ) {
        amr_lv = lv_max - 3;
    }
    else if ( fabs(x) < 30.0*H && fabs(y) <= 10.0*H && fabs(z) <= 10.0*H ) {
        amr_lv = lv_max - 4;
    }
    else {
        amr_lv = lv_max - 5;
    }
#endif


    if (amr_lv < 0) { amr_lv = 0; }

    return  amr_lv;
}


real InitMonitorLevelset::
func_amr_level_flow_map(
    const real x,
    const real y,
    const real z,
    const real dx
    )
{
    constexpr int  lv_max = DefAMR::LV_MAX;

//    constexpr real height0 = 135.0 + 4.0 - 0.1;
//    constexpr real height1 = 768.0 + 4.0 - 0.1;

//    constexpr real height0 =  150.0 + 4.0 - 0.1;
//    constexpr real height1 =  700.0 + 4.0 - 0.1;

//    // oklahoma 1m //
//    constexpr real tmp_dx  =   1.0;
//    constexpr real height0 =  64.0    -   4.0 - tmp_dx*3;
//    constexpr real height1 =  height0 + 128.0 - tmp_dx*2*3;

    // oklahoma 2m //
    constexpr real tmp_dx  =   2.0;
    constexpr real height0 =  96.0    -   8.0 - tmp_dx*3;
//    constexpr real height0 =  128.0    -   8.0 - tmp_dx*3;
    constexpr real height1 =  height0 + 256.0 - tmp_dx*2*3;

//    // oklahoma 4m //
//    constexpr real tmp_dx  =    4.0;
//    constexpr real height0 =  256.0   -  16.0 - tmp_dx*3;
//    constexpr real height1 =  height0 + 512.0 - tmp_dx*2*3;


    int amr_lv;
    // lv_max = 3 //
    if      ( fabs(z) <= height0 ) {
        amr_lv = lv_max - 1;
    }
    else if ( fabs(z) <= height1 ) {
        amr_lv = lv_max - 2;
    }
    else {
        amr_lv = lv_max - 3;
    }

//    if ( (x < 120 || x > 1800) && amr_lv == lv_max-1 ) {
//        amr_lv = lv_max - 2;
//    }
//
//    if ( (x < 60 || x > 1860) && amr_lv == lv_max-2 ) {
//        amr_lv = lv_max - 3;
//    }


    if (amr_lv < 0) { amr_lv = 0; }

    return  amr_lv;
}


real InitMonitorLevelset::
func_amr_level_channel_flow(
    const real x,
    const real y,
    const real z,
    const real dx
    )
{
    constexpr int  lv_max = DefAMR::LV_MAX;
    int amr_lv;
    const real length = fabs( 1.0 - fabs(z) );


    amr_lv = 0;
#if 0
    // lv_max = 2 //
//    constexpr real height0 = 0.02;
    constexpr real height0 = 0.1;
//    constexpr real height0 = 0.2;
//    constexpr real height0 = 0.6;
    if      ( length <= height0 ) {
        amr_lv = lv_max - 1;
    }
    else {
        amr_lv = lv_max - 2;
    }
#else
    // lv_max = 3 //
    constexpr real height0 = 0.10;
    constexpr real height1 = 0.25;
    if      ( length <= height0 ) {
        amr_lv = lv_max - 1;
    }
    else if ( length <= height1 ) {
        amr_lv = lv_max - 2;
    }
    else {
        amr_lv = lv_max - 3;
    }
#endif


    if (amr_lv < 0) { amr_lv = 0; }

    return  amr_lv;
}
