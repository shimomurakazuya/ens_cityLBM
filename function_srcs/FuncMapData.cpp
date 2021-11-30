#include "FuncMapData.h"
#include "defineAMR.h"
#include "defineFilenames.h"
#include "Index.h"
#include "FuncObj.h"
#include "FuncMath.h"
#include "runtime_error.hpp"


// now debugging //

void  MapData::levelset_map(
          real* lv_obj,
    std::function<double(int)> x,
    std::function<double(int)> y,
    std::function<double(int)> z,
    const int   nx_leaf,
    const int   lv,
    const real  dx
) const {
    // levelset //
    for (int k=0; k<nx_leaf; k++) {
    for (int j=0; j<nx_leaf; j++) {
    for (int i=0; i<nx_leaf; i++) {
        const int    id = Index::id(i,j,k);
        const double xx = x(id);
        const double yy = y(id);
        const double zz = z(id);
        try {
        lv_obj[id] = get_levelset(xx, yy, zz, lv, dx);
        } catch(std::out_of_range) {
            std::cerr << __PRETTY_FUNCTION__ << " map failed at " << xx << " " << yy << " " << zz <<  std::endl;
            lv_obj[id] = -zz;
        }
    }
    }
    }
}


double  MapData::get_levelset(
    const double xc,
    const double yc,
    const double zc,
    const int    lv,
    const double dx
) const {
    constexpr double  ground_height  = 0.0;

    // const auto[i,j] = convert_geometory_to_index(xc, yc); // C++17
    const auto ix = convert_geometory_to_index(xc, yc).first,
               iy = convert_geometory_to_index(xc, yc).second;

    double  levelset;
    // ix < 0 or ix > nx-1 //
    if ( !check_index_in_map(ix, iy) ) {
        return  ground_height - zc;
    }

    // out of domain //
    if ( xc < west_  || xc > east_
      || yc < south_ || yc > north_) {
        return  ground_height - zc;
    }


    auto  levelset_neighbor_in_fluid = [&](double xx, double yy, double dx, double dd, double levelset)
    {
        // const auto[ixn,jxn] = convert_geometory_to_index(xc, yc); // C++17
        const auto  ixn = convert_geometory_to_index(xx, yy).first,
                    iyn = convert_geometory_to_index(xx, yy).second;

        double  height_tmp = check_index_in_map(ixn, iyn) ?  get_height(ixn, iyn) : ground_height;

        double levelset_tmp;
        if   ( zc < height_tmp ) { levelset_tmp = -dd + 0.5*dx; }
        else                     { levelset_tmp = -sqrt( pow( zc - height_tmp, 2 ) + pow( -dd + 0.5*dx, 2 ) ); }

        runtime_assert(FuncObj::is_fluid(levelset_tmp), "InternalError: invalid levelset value; may be round error of float?");

        return  FuncMath::max( levelset, levelset_tmp );
    };

    auto  levelset_neighbor_in_object = [&](double xx, double yy, double dx, double dd, double levelset)
    {
        // const auto[ixn,jyn] = convert_geometory_to_index(xc, yc); // C++17
        const auto  ixn = convert_geometory_to_index(xx, yy).first,
                    iyn = convert_geometory_to_index(xx, yy).second;

        double  height_tmp = check_index_in_map(ixn, iyn) ?  get_height(ixn, iyn) : ground_height;

        double levelset_tmp;
        if   ( zc > height_tmp ) { levelset_tmp =  dd - 0.5*dx; }
//        else                     { levelset_tmp =  sqrt( pow( zc - height_tmp, 2) + dd*dd ); }
        else                     { levelset_tmp =  sqrt( pow( zc - height_tmp, 2 ) + pow( dd - 0.5*dx, 2 ) ); }

        runtime_assert(FuncObj::is_obj(levelset_tmp), "InternalError: invalid levelset value; maybe round error of float?");

        return  FuncMath::min( levelset, levelset_tmp );
    };


    // center //
    const double heightc = get_height(ix, iy);
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
