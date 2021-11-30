#pragma once
#ifndef MESHVALUE_H_
#define MESHVALUE_H_


#include <iostream>
#include <cstdlib>

#include "defineMemory.h"
#include "Grid.h"
#include "Tree.h"
#include "Coordinates.h"
// values //
#include "ValueNS.h"
#include "ValueObjLS.h"
#include "ValueLBM.h"

class  MeshValue {
private:
    MemType memType_{MemType::NullPtr};

    int     nn_max_;
    int     n_scalars_;

    // coordinates xyz //
    Coordinates     coordinates_;

    // values //
    ValueObjLS      valueObjLS_;
    ValueNS         valueNS_;
    ValueLBM        valueLBM_;

public:
    MeshValue ()  {
        memType_ = MemType::Host;
#ifdef GPU_CALCULATION__
        memType_ = MemType::Managed;
#endif
    } // default //
    ~MeshValue () {}

public:
    int  nn_max() const { return  nn_max_; }
    int  n_scalars() const { return n_scalars_; }

    // coordinates xyz //
          Coordinates& coordinates()          { return  coordinates_; }
    const Coordinates& coordinates()    const { return  coordinates_; }

    // values //
          ValueObjLS& valueObjLS()          { return  valueObjLS_; }
    const ValueObjLS& valueObjLS()    const { return  valueObjLS_; }

          ValueNS& valueNS()          { return  valueNS_; }
    const ValueNS& valueNS()    const { return  valueNS_; }

          ValueLBM& valueLBM()          { return  valueLBM_; }
    const ValueLBM& valueLBM()    const { return  valueLBM_; }

public:
    void init(const int  nn_max, const int n_scalars);
    void set_coordinate(const int  lv, const Grid& grid, const Tree& tree);

    void copy_MeshValue(
        const MeshValue&    other,
        const Grid&         grid,
        const Tree&         tree,
        const bool          is_reset
        );

    void mask_meshValue_wo_halo(const Tree& tree, const int lv);
    void mask_meshValue_wo_obj (const Tree& tree, const int lv);

private:
    void  init_coordinates(const int  nn_max);
    void  init_mesh_values(const int  nn_max, const int n_scalars);


};


#endif
