#include "Grid.h"
#include "defineFilenames.h"


// public //
void  Grid::
copy(const Grid& grid)
{
    nx_     = grid.nx();
    offset_ = grid.offset();
    dx_     = grid.dx();
}


void  Grid::
writeGrid(const std::string  filename)
const
{
    const stIOGrid  stIOgrid{
                        nx_.x(), nx_.y(), nx_.z(),
                        offset_.x(), offset_.y(), offset_.z(),
                        dx_
                        };

    // write //
    std::ofstream  fout;
    fout.open(filename, std::ios::binary);
    fout.write( ( char * ) &stIOgrid, sizeof( stIOGrid ) );
    fout.close();
}


void  Grid::
readGrid (const std::string  filename)
{
    stIOGrid  stIOgrid;

    // read //
    std::ifstream  fin;
    fin.open(filename, std::ios::binary);
    fin.read( ( char * ) &stIOgrid, sizeof( stIOGrid ) );
    fin.close();


    // update //
    nx_.    set_xyz( stIOgrid.nx,       stIOgrid.ny,       stIOgrid.nz       );
    offset_.set_xyz( stIOgrid.offset_x, stIOgrid.offset_y, stIOgrid.offset_z );

    dx_ = stIOgrid.dx;
}


void  Grid::
cout(const std::string& str)
const
{
    const Vector3d<real>  dx(dx_,dx_,dx_);

    std::cout << str << std::endl;
    offset_. cout("offset = ");
    nx_.     cout("nx     = ");
    dx.      cout("dx     = ");
    length().cout("length = ");
}


std::string  Grid::
filename(const int lv, const int rank, const int step)
const
{
    return    Foldernames::io_folder + "/"
            + Filenames::grid_name0
            + "-lv"   + std::to_string(lv)
            + "-rank" + std::to_string(rank)
            + "-step" + std::to_string(step)
            + ".dat";
}
