#include "PostprocessNaturalConvection2d.h"
#include "Index.h"
#include "IndexLBM.h"
#include "FuncObj.h"
#include "defineFluidProperty.h"
#include <algorithm>


void  PostprocessNaturalConvection2d::OutputStatistics(int t, Field& field, const MeshValue*  meshValues)
{
//    const Parameters parameters = field.parameters();

    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        dx_  [lv] = meshValues[lv].coordinates().dx();
    }

    double  uT  = 0.0;
    double  Tx  = 0.0;
    double  vol = 0.0;
    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        double  _uT;
        double  _Tx;
        double  _vol;
        Nusselt_number_lv( _uT, _Tx, _vol,  lv, field, meshValues, field.taskID().id_tasks_lbm(lv) );

        uT  += _uT;
        Tx  += _Tx;
        vol += _vol;
    }
    uT = uT/vol;
    Tx = Tx/vol;

    uT_ = uT;
    Tx_ = Tx;
    Nu_ = uT_ - Tx_;

    output_values(t);
}


// private //
void  PostprocessNaturalConvection2d::Nusselt_number_lv(
          double&    uT,
          double&    Tx,
          double&    vol,
    const int        lv,
    const Field&     field,
    const MeshValue* meshValues,
    const TaskID::vector_type& id_tasks
    )
{
    const real* u     = meshValues[lv].valueNS().u();
//    const real* v     = meshValues[lv].valueNS().v();
//    const real* w     = meshValues[lv].valueNS().w();

    const real* T     = meshValues[lv].valueNS().T();

    const real*  lv_obj    = field.meshValue(lv).valueObjLS().lv_obj();
    const real*  u_obj     = field.meshValue(lv).valueObjLS().u_obj();
    const real*  T_obj     = field.meshValue(lv).valueObjLS().T_obj();

    const Tree& tree = field.tree();
    const int*  mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();

    const Parameters& parameters = field.parameters();
    const double dx     = dx_[lv];
    const double c_ref  = parameters.c_ref_lbm();

    constexpr int  nQ = 27;

    // sum //
    double  vol_all = 0.0;
    double  _uT_all  = 0.0;
    double  _Tx_all  = 0.0;
//    double  _Nu_all  = 0.0;
    for (int l=0; l<(int)id_tasks.size(); l++) {
        if ( !tree.nodes(l)->nodeCalFlags().Cal() ) { continue; }

        const Array3D<int>  offsets = tree.nodes(l)->neighbor_mesh_offsets();

        for (int k=0; k<DefAMR::NX_LEAF; k++) {
        for (int j=0; j<DefAMR::NX_LEAF; j++) {
        for (int i=0; i<DefAMR::NX_LEAF; i++) {
            const int  idl = id_tasks[l];

            int  offset3d[27];
            for (int idv=0; idv<nQ; idv++) {
                offset3d[idv] = mesh_offsets3x3x3[idl*nQ + idv];
            }

            int  ids[27];
            for (int kv=-1; kv<=1; kv++) {
                for (int jv=-1; jv<=1; jv++) {
                    for (int iv=-1; iv<=1; iv++) {
                        const int  idv_leaf = IndexLBM::idv_lbm(iv, jv, kv);

                        ids[idv_leaf] = Index::id(i+iv,j+jv,k+kv, offset3d);
                    }
                }
            }


            // is object //
            const int  id0 = Index::id0(ids);
            if ( FuncObj::is_obj( lv_obj[id0] ) ) { continue; }


            // statistics //
            const double lv3[]       = { lv_obj[Index::id(ids, -1,  0,  0)], lv_obj[Index::id(ids,  0,  0,  0)], lv_obj[Index::id(ids,  1,  0,  0)] };

            // Tx //
//            const double  Tx = (double)( T[Index::id(ids,  1,  0,  0)] - T[Index::id(ids, -1,  0,  0)] ) / (2.0*dx);

            const double T3[]        = { T[Index::id(ids, -1,  0,  0)], T[Index::id(ids,  0,  0,  0)], T[Index::id(ids,  1,  0,  0)] };
            const double Tobj3[]     = { T_obj[Index::id(ids, -1,  0,  0)], T_obj[Index::id(ids,  0,  0,  0)], T_obj[Index::id(ids,  1,  0,  0)] };

            const double  _T[] = { (FuncObj::is_obj( lv3[0] )) ? Tobj3[0] : (T3[0] + T3[1])*0.5,
                                   (FuncObj::is_obj( lv3[2] )) ? Tobj3[2] : (T3[2] + T3[1])*0.5  };

            const double Tx = (double)(_T[1] - _T[0]) / (dx);
            const double  _Tx = air_property::xi*Tx;
//            const double  _Nu = ( us*Ts - air_property::xi*Tx );

            // uT //
//            const double  us =   u[Index::id(ids,  0,  0,  0)];
//            const double  us =   u[Index::id(ids,  0,  0,  0)] * c_ref;
//            const double  Ts =   T[Index::id(ids,  0,  0,  0)];
//            const double  _uT = us*Ts;

            const double u3[]        = { u[Index::id(ids, -1,  0,  0)], u[Index::id(ids,  0,  0,  0)], u[Index::id(ids,  1,  0,  0)] };
            const double uobj3[]     = { u_obj[Index::id(ids, -1,  0,  0)], u_obj[Index::id(ids,  0,  0,  0)], u_obj[Index::id(ids,  1,  0,  0)] };
            const double  _u[] = { FuncObj::is_obj( lv3[0] ) ? uobj3[0] : (u3[0] + u3[1])*0.5,
                                   FuncObj::is_obj( lv3[2] ) ? uobj3[2] : (u3[2] + u3[1])*0.5  };
//            const double _uT  = (_u[0]       + _u[1]      ) * (double)0.5 * c_ref;
            const double _uT  = (_u[0]*_T[0] + _u[1]*_T[1]) * (double)0.5 * c_ref;


            // total //
//            const double  sur = (dx*dx);
            const double  vol = (dx*dx*dx);
            _uT_all += _uT * vol;
            _Tx_all += _Tx * vol;
//            _Nu_all += _Nu * vol;
            vol_all += vol;
        }
        }
        }
    }
    double _uT_all_g = 0.0;
    double _Tx_all_g = 0.0;
    double vol_all_g = 0.0;

    MPI_Allreduce(&_uT_all, &_uT_all_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&_Tx_all, &_Tx_all_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&vol_all, &vol_all_g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    const double dT = 1.0;

    uT  = (_uT_all_g) / (dT*air_property::xi);
    Tx  = (_Tx_all_g) / (dT*air_property::xi);
    vol = vol_all_g;
}


void  PostprocessNaturalConvection2d::output_values(int step)
{
    if (rank_ != 0) { return; }

    // fout //
    std::ofstream fout;
    fout.open(filename(step));

    fout << "viscosity,Nu" << std::endl;
    fout << air_property::Viscosity << "," << Nu_ << std::endl;
    fout.close();


    std::cout << "***************************\n";
    std::cout << "viscosity, Nu ( uT, -Tx )" << std::endl;
    std::cout << air_property::Viscosity << ", " << Nu_ << " ( " << uT_ << ", " << -Tx_ << " ) " << std::endl;
    std::cout << "***************************\n";
}


