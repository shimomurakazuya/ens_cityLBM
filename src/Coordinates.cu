#include "Coordinates.h"
#include <algorithm>
#include "FuncAllocate.h"
#include "FuncMath.h"
#include "FuncLoop.h"


// public //
void  Coordinates::
init(const int  nn_max, const enum MemType  memType)
{
    memType_ = memType;

    allocate(nn_max);
}


void  Coordinates::
set_uniform(const int  lv, const Grid& grid, const Tree& tree)
{
    // initialize uniform mesh //
    int  id_offset = 0;

    auto func_init =  [lv, &id_offset, &grid, this](const Node& node)
    {
        if (lv == node.level()) {
            const real           dx = grid.dx() / (real)DefAMR::NX_LEAF;
            const Vector3d<int>  nx(DefAMR::NX_LEAF, DefAMR::NX_LEAF, DefAMR::NX_LEAF);
            const Vector3d<real> x0(grid.offset());

            const Vector3d<real> ii(node.x(), node.y(), node.z());
            const Vector3d<real> xx(x0 + ii*dx*DefAMR::NX_LEAF);

            set_xyz(id_offset, nx, xx, dx);

//            std::cout << "i, j, k = " << node.x() << ", " << node.y() << ", " << node.z() << std::endl;

            id_offset += pow(DefAMR::NX_LEAF, 3);
        }
    };

    std::for_each( tree.nodes().begin(), tree.nodes().end(), func_init );
}


void  Coordinates::
set_x_from_xp(real* x, const real* xp, const int nx, const int ny, const int nz)
const
{
    FuncLoop::Loop3d  loop3d(0,nx, 0,ny, 0,nz);
    loop3d.for_each(
            [&x, &xp, nx,ny,nz](const int i, const int j, const int k)
            {
                if (i < nx && j < ny && k < nz) {
                    const int  id  = i + nx*j + nx*ny*k;
                    const int  idp = i + (nx+1)*j + (nx+1)*(ny+1)*k;

                    x[id] = xp[idp];
                }
            }
            );

}


void  Coordinates::
copy (const int  nn_max, const Coordinates&  other)
{
    FuncMath::copy_array( x(), other.x(), nn_max );
    FuncMath::copy_array( y(), other.y(), nn_max );
    FuncMath::copy_array( z(), other.z(), nn_max );

    dx_ = other.dx();
    dy_ = other.dy();
    dz_ = other.dz();
}


// private //
void  Coordinates::
set_xyz(const int  offset, const Vector3d<int>&  nx, const Vector3d<real>&  x0, const real  dx)
{
    auto  index = [offset, &nx](const int i, const int j, const int k){ return  offset + i + nx.x()*j + nx.x()*nx.y()*k; };

    // uniform & cube //
    dx_ = dx;
    dy_ = dx;
    dz_ = dx;

    for (int k=0; k<nx.z(); k++) {
        for (int j=0; j<nx.y(); j++) {
            for (int i=0; i<nx.x(); i++) {
                x_[index(i,j,k)] = x0.x() + dx_*i;
                y_[index(i,j,k)] = x0.y() + dy_*j;
                z_[index(i,j,k)] = x0.z() + dz_*k;
            }
        }
    }
}


void  Coordinates::
allocate(const int nn_max)
{
    FuncAllocate::allocate_value<real>(&x_, nn_max, memType_);
    FuncAllocate::allocate_value<real>(&y_, nn_max, memType_);
    FuncAllocate::allocate_value<real>(&z_, nn_max, memType_);
}


void  Coordinates::
release ()
{
    FuncAllocate::release_value<real>(x_, memType_);
    FuncAllocate::release_value<real>(y_, memType_);
    FuncAllocate::release_value<real>(z_, memType_);
}
