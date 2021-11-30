#include "PostprocessTaylorGreen.h"
#include "Index.h"
#include "IndexLBM.h"
#include "defineFluidProperty.h"
#include "FuncLBM.h"
#include "InitTaylorGreen.h"


void  PostprocessTaylorGreen::OutputStatistics(int t, float time, Field& field, const MeshValue*  meshValues)
{
    if (rank_ == 0) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    std::cout << "t = " << time << std::endl;


    const double ReNum = air_property::ReNum;
    InitTaylorGreen initTaylorGreen(ReNum);

    const int  nQ  = LBM_velocity_model::nQ;

    const real  c_ref  = field.parameters().c_ref_lbm();

    const Tree& tree   = field.tree();
    const int   n_leaf = tree.number_of_nodes();

    double  energy_sum = 0.0;

    double  error_vel_sum = 0.0;
    double  norm_vel_sum  = 0.0;

    double  error_dp_sum = 0.0;
    double  norm_dp_sum  = 0.0;
//    int    num_grids = 0;

    const Parameters parameters = field.parameters();
    const double zwidth = meshValues[0].coordinates().dx() * DefAMR::NX_LEAF;

    for (int l=0; l<n_leaf; l++) {
        const Node*  node = tree.nodes(l);
        if ( !node->nodeCalFlags().Cal() ) { continue; }

        const int  lv     = tree.nodes(l)->level();
        const int  offset = tree.nodes(l)->mesh_offset();

        // coordinates xyz //
        const auto& coordinates = meshValues[lv].coordinates();
        const real   dx = meshValues[lv].coordinates().dx();

        // LBM //
        const real* f_lbm = field.meshValue(lv).valueLBM().f_lbm();

        for (int k=0; k<DefAMR::NX_LEAF; k++) {
        for (int j=0; j<DefAMR::NX_LEAF; j++) {
        for (int i=0; i<DefAMR::NX_LEAF; i++) {
            const int  id          = Index   ::id0      (i,j,k, offset);
            const int  id_lbm_base = IndexLBM::id0_base0(i,j,k, offset);

            real  fs[nQ];
            for (int kv=-1; kv<=1; kv++) {
                for (int jv=-1; jv<=1; jv++) {
                    for (int iv=-1; iv<=1; iv++) {
                        const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);
                        const int  id_lbm   = IndexLBM::id0(id_lbm_base, iv,jv,kv);

                        fs[idv_leaf] = f_lbm[id_lbm];
                    }
                }
            }

            real  rho, u,v,w;
            FuncLBM::velocity_lbm(fs, rho, u,v,w);
            const real dp = FuncLBM::dp_lbm(rho-1.0, c_ref);

            const auto xid = coordinates.xcell(id);
            const auto yid = coordinates.ycell(id);
            const real ua  = initTaylorGreen.ua(xid, yid, time)/c_ref;
            const real va  = initTaylorGreen.va(xid, yid, time)/c_ref;
            const real dpa = initTaylorGreen.dpa(xid, yid, time);

            energy_sum += 0.5*rho*(pow(u, 2) + pow(v, 2)) * c_ref*c_ref * pow(dx, 3) / (zwidth);

            error_vel_sum += ( pow(u-ua, 2) + pow(v-va, 2) ) * pow(dx, 3);
            norm_vel_sum  += ( pow(ua,   2) + pow(va,   2) ) * pow(dx, 3);

            error_dp_sum += ( pow(dp-dpa, 2) ) * pow(dx, 3);
            norm_dp_sum  += ( pow(dpa,    2) ) * pow(dx, 3);

//            num_grids++;
        }
        }
        }

    }

    const real  ep = 1.0e-14;
    const real  error_vel_norm = sqrt( error_vel_sum / (norm_vel_sum+ep) );
    const real  error_dp_norm  = sqrt( error_dp_sum  / (norm_dp_sum+ep) );

    std::cout << "time, energy_sum, error_vel_norm, error_dp_norm = " << time << ", " << energy_sum << ", " << error_vel_norm << ", " << error_dp_norm << std::endl;

    static bool initialized = false;

    std::ofstream fout;
    if (!initialized) { 
        fout.open(filename(t));
        fout << "time,energy_sum,error_vel_norm,error_dp_norm" << std::endl;
    }
    else {
        fout.open(filename(t), std::ofstream::app);
    }

    fout << time << "," << energy_sum << "," << error_vel_norm << "," << error_dp_norm << std::endl;

    fout.close();

    initialized = true;
}
