#ifndef PBVR_OUTPUT_H_INCLUDED
#define PBVR_OUTPUT_H_INCLUDED

// Thin wrapper around PBVR InSituLib (ensemble_generate_particles).
// Confines vismodule / FunctionParser headers to pbvr_output.cpp so the PBVR
// `class Node` (FunctionParser) never collides with cityLBM's global
// `class Node` (include/Node.h). Only plain types appear in this header.
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
    int ncells );

#endif // PBVR_OUTPUT_H_INCLUDED
