#include "FuncMapData.h"
#include "defineAMR.h"
#include "defineFilenames.h"
#include "Index.h"
#include "FuncObj.h"
#include "FuncMath.h"


// now debugging //
namespace  FuncMapData {


void  levelset_map(
          real* lv_obs,
    const real* x,
    const real* y,
    const real* z,
    const int   nx_leaf,
    const int   lv,
    const real  dx,
    const int*  height_map,
    const int   mx,
    const int   my,
    const real  dx_map,
    const real  map_offset_x,
    const real  map_offset_y,
    const real  domain_min_x,
    const real  domain_min_y,
    const real  domain_max_x,
    const real  domain_max_y
    )
{
    // levelset //
    for (int k=0; k<nx_leaf; k++) {
    for (int j=0; j<nx_leaf; j++) {
    for (int i=0; i<nx_leaf; i++) {
        const int    id = Index::id(i,j,k);
        const double xx = x[id] + 0.5*dx; // cell based amr //
        const double yy = y[id] + 0.5*dx;
        const double zz = z[id] + 0.5*dx;

        lv_obs[id] = get_levelset(xx, yy, zz, lv, dx, height_map, mx, my, dx_map, map_offset_x, map_offset_y, domain_min_x, domain_min_y, domain_max_x, domain_max_y);

//        if ( check_index_in_map(ix, iy, mx, my) ) {
//            std::cout << x[id] << ", " << y[id] << ", " << z[id] << " : " << ix << ", " << iy << " = " << height_map[ix + mx*iy] << ", " << lv_obs[id] << std::endl;
//        }
    }
    }
    }
}


double  get_levelset(
    const double xc,
    const double yc,
    const double zc,
    const int    lv,
    const double dx,
    const int*   height_map,
    const int    mx,
    const int    my,
    const double dx_map,
    const double map_offset_x,
    const double map_offset_y,
    const double domain_min_x,
    const double domain_min_y,
    const double domain_max_x,
    const double domain_max_y
    )
{
    const double  ground_height  = 0.0;

    int  ix, iy; // map index @ lv_max //
    convert_geometory_to_index(
        ix, iy,
        xc, yc,
        dx_map,
        map_offset_x, map_offset_y
        );

    double  levelset;
    // ix < 0 or ix > nx-1 //
    if ( !check_index_in_map(ix, iy, mx, my) ) {
        return  ground_height - zc;
    }

    // out of domain //
    if ( xc < domain_min_x || xc > domain_max_x
      || yc < domain_min_y || yc > domain_max_y ) {
        return  ground_height - zc;
    }


    auto  levelset_neighbor_in_fluid = [&height_map, &mx, &my, &dx_map, &map_offset_x, &map_offset_y, &zc](double xx, double yy, double dx, double dd, double levelset)
    {
        int ixn, iyn;
        convert_geometory_to_index(
            ixn, iyn,
            xx, yy,
            dx_map,
            map_offset_x, map_offset_y
            );

        double  height_tmp = check_index_in_map(ixn, iyn, mx, my) ?  get_height(height_map, mx, my, ixn, iyn) : ground_height;

        double levelset_tmp;
        if   ( zc < height_tmp ) { levelset_tmp = -dd + 0.5*dx; }
        else                     { levelset_tmp = -sqrt( pow( zc - height_tmp, 2 ) + pow( -dd + 0.5*dx, 2 ) ); }

        if ( levelset_tmp >= 0.0 ) { std::cout << __PRETTY_FUNCTION__ << " : " << "error\n"; exit(-1); }

        return  FuncMath::max( levelset, levelset_tmp );
    };

    auto  levelset_neighbor_in_object = [&height_map, &mx, &my, &dx_map, &map_offset_x, &map_offset_y, &zc](double xx, double yy, double dx, double dd, double levelset)
    {
        int ixn, iyn;
        convert_geometory_to_index(
            ixn, iyn,
            xx, yy,
            dx_map,
            map_offset_x, map_offset_y
            );

        double  height_tmp = check_index_in_map(ixn, iyn, mx, my) ?  get_height(height_map, mx, my, ixn, iyn) : ground_height;

        double levelset_tmp;
        if   ( zc > height_tmp ) { levelset_tmp =  dd - 0.5*dx; }
//        else                     { levelset_tmp =  sqrt( pow( zc - height_tmp, 2) + dd*dd ); }
        else                     { levelset_tmp =  sqrt( pow( zc - height_tmp, 2 ) + pow( dd - 0.5*dx, 2 ) ); }

        if ( levelset_tmp < 0.0 ) { std::cout << __PRETTY_FUNCTION__ << " : " << "error\n"; exit(-1); }

        return  FuncMath::min( levelset, levelset_tmp );
    };


    // center //
    const double heightc = (double)get_height(height_map, mx, my, ix, iy);
    levelset = heightc - zc;
//    if (levelset < 0.0) { levelset = -0.5*dx; }
//    else                { levelset =  0.5*dx; }

//    std::cout << "lv : ix, iy, xc, yc, zc, heightc = " << lv << " : " << ix << ", " << iy << ", " << xc << ", " << yc << ", " << zc << ", " << heightc << std::endl;

    // neighbor //
    if ( FuncObj::is_fluid(levelset) ) {
        for (int jj=-1; jj<=1; jj++) {
        for (int ii=-1; ii<=1; ii++) {
            if (ii == 0 && jj == 0) { continue; }

            levelset = levelset_neighbor_in_fluid(xc+ii*dx, yc+jj*dx, dx, sqrt(ii*ii + jj*jj)*dx, levelset);
        }
        }
    }
    else {
        for (int jj=-1; jj<=1; jj++) {
        for (int ii=-1; ii<=1; ii++) {
            if (ii == 0 && jj == 0) { continue; }

            levelset = levelset_neighbor_in_object(xc+ii*dx, yc+jj*dx, dx, sqrt(ii*ii + jj*jj)*dx, levelset);
        }
        }

    }

//    std::cout << "lv : ix, iy, xc, yc, zc, heightc, levelset = " << lv << " : " << ix << ", " << iy << ", " << xc << ", " << yc << ", " << zc << ", " << heightc << ", " << levelset << std::endl;

    return levelset;
}


int  get_height(
    const int* height_map,
    const int mx,
    const int my,
    const int ix,
    const int iy
    )
{
    return  height_map[ix + mx*iy];
//    return  0.0;
}


bool  check_index_in_map(
    const int  ix,
    const int  iy,
    const int  nx,
    const int  ny
    )
{
    return ( ix < 0 || ix >= nx || iy < 0 || iy >= ny ) ? false : true;
}


void  convert_geometory_to_index(
          int& ix,
          int& iy,
    const double x,
    const double y,
    const double dx,
    const double map_offset_x,
    const double map_offset_y
    )
{
//    std::cout << x << ", " << map_offset_x << ", " << dx << std::endl;
    ix = std::floor( ( x - map_offset_x ) / dx );
    iy = std::floor( ( y - map_offset_y ) / dx );
}


};
