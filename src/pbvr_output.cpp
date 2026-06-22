#include "pbvr_output.h"

// All PBVR / vismodule / FunctionParser headers are confined to this TU,
// which never includes cityLBM headers (avoids the global `Node` clash).
#include "kvs_wrapper.h"

void pbvr_generate_particles(
    int time_step,
    int num_ensemble,
    float dom_min_x, float dom_min_y, float dom_min_z,
    float dom_max_x, float dom_max_y, float dom_max_z,
    float** values,
    int nvariables,
    float* coordinates,
    int ncoords,
    unsigned int* connections,
    int ncells )
{
    domain_parameters_unstruct dom = {
        dom_min_x, dom_min_y, dom_min_z,
        dom_max_x, dom_max_y, dom_max_z
    };
    ensemble_generate_particles(
        time_step, num_ensemble, dom,
        values, nvariables,
        coordinates, ncoords,
        connections, ncells,
        vismodule::VolumeObjectBase::CellType::Hexahedra );
}
