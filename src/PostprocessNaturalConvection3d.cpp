#include "PostprocessNaturalConvection3d.h"
#include "Index.h"
#include "defineFluidProperty.h"
#include <algorithm>


void  PostprocessNaturalConvection3d::OutputStatistics(int t, Field& field, const MeshValue*  meshValues)
{
    const Parameters parameters = field.parameters();
    const real  z_min = parameters.coefDomain().z_global_domain_min;
    const real  z_max = parameters.coefDomain().z_global_domain_max;

    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        dx_  [lv] = meshValues[lv].coordinates().dx();
        numz_[lv] = (z_max-z_min)/dx_[lv];

        if (rank_ == 0) {
            std::cout << "lv, numz, dx = " << lv << ", " << numz_[lv] << ", " << meshValues[lv].coordinates().dx() << std::endl;
        }
        values_  [lv].resize(numz_[lv]);
        values_g_[lv].resize(numz_[lv]);
    }
    init_value_cf();

    for (int lv=0; lv<DefAMR::LV_MAX; lv++) { exist_flag ( lv, field, meshValues, field.taskID().id_tasks_lbm(lv) ); }
    for (int lv=0; lv<DefAMR::LV_MAX; lv++) { mean_values( lv, field, meshValues, field.taskID().id_tasks_lbm(lv) ); }
    for (int lv=0; lv<DefAMR::LV_MAX; lv++) { rms_values ( lv, field, meshValues, field.taskID().id_tasks_lbm(lv) ); }

    output_values(t);
}


// private //
void  PostprocessNaturalConvection3d::exist_flag(
    const int        lv,
    const Field&     field,
    const MeshValue* meshValues,
    const TaskID::vector_type& id_tasks
    )
{
    const auto& coordinates = meshValues[lv].coordinates();
    const real*  lv_obj = field.meshValue(lv).valueObjLS().lv_obj();

    const Tree& tree = field.tree();
    const int*  mesh_offsets = tree.mesh_offsets();

    const Parameters& parameters = field.parameters();
    const real dx     = dx_[lv];
    const real c_ref  = parameters.c_ref_lbm();

    for (int l=0; l<(int)id_tasks.size(); l++) {
        for (int k=0; k<DefAMR::NX_LEAF; k++) {
        for (int j=0; j<DefAMR::NX_LEAF; j++) {
        for (int i=0; i<DefAMR::NX_LEAF; i++) {
            const int  idl    = id_tasks[l];
            const int  offset = mesh_offsets[idl];
            const int  id     = Index::id0(i,j,k, offset);

            const real xcell = coordinates.xcell(id);
            const real ycell = coordinates.ycell(id);
            const real zcell = coordinates.zcell(id);
            if ( zcell < z_min_nc3d_ || zcell > z_max_nc3d_ ) { continue; }
            if ( !is_center( xcell, ycell ) )                 { continue; }

            int  idz;
            real z_obj;
            get_idz_zobj_nc3d(idz, z_obj, zcell, dx, lv_obj[id]);
            if (idz >= numz_[lv]) { std::cout << __PRETTY_FUNCTION__ << " : error idz\n"; }

            add_exist_flag(lv, idz, z_obj);
        }
        }
        }
    }

    reduce_exist_flag(lv);
    reduce_z(lv);
    average_num(lv);
}


void  PostprocessNaturalConvection3d::mean_values(
    const int        lv,
    const Field&     field,
    const MeshValue* meshValues,
    const TaskID::vector_type& id_tasks
    )
{
    const real* u     = meshValues[lv].valueNS().u();
    const real* v     = meshValues[lv].valueNS().v();
    const real* w     = meshValues[lv].valueNS().w();
    const real* T     = meshValues[lv].valueNS().T();

    const auto& coordinates = meshValues[lv].coordinates();

    const real*  lv_obj = field.meshValue(lv).valueObjLS().lv_obj();

    const Tree& tree = field.tree();
    const int*  mesh_offsets = tree.mesh_offsets();

    const Parameters& parameters = field.parameters();
    const real dx     = dx_[lv];
    const real c_ref  = parameters.c_ref_lbm();

    for (int l=0; l<(int)id_tasks.size(); l++) {
        for (int k=0; k<DefAMR::NX_LEAF; k++) {
        for (int j=0; j<DefAMR::NX_LEAF; j++) {
        for (int i=0; i<DefAMR::NX_LEAF; i++) {
            const int  idl = id_tasks[l];

            const int  offset = mesh_offsets[idl];
            const int  id     = Index::id0(i,j,k, offset);

            const real xcell = coordinates.xcell(id);
            const real ycell = coordinates.ycell(id);
            const real zcell = coordinates.zcell(id);
            if (zcell < z_min_nc3d_ || zcell >= z_max_nc3d_) { continue; }
            if ( !is_center( xcell, ycell ) )                { continue; }

            int  idz;
            real z_obj;
            get_idz_zobj_nc3d(idz, z_obj, zcell, dx, lv_obj[id]);
            if (idz >= numz_[lv]) { std::cout << __PRETTY_FUNCTION__ << " : error idz\n"; }

            add_mean_values(lv, idz, z_obj, u[id]*c_ref, v[id]*c_ref, w[id]*c_ref, T[id]);
        }
        }
        }
    }

    average_mean_values(lv);
}


void  PostprocessNaturalConvection3d::rms_values(
    const int        lv,
    const Field&     field,
    const MeshValue* meshValues,
    const TaskID::vector_type& id_tasks
    )
{
    const real* u     = meshValues[lv].valueNS().u();
    const real* v     = meshValues[lv].valueNS().v();
    const real* w     = meshValues[lv].valueNS().w();

    const auto& coordinates = meshValues[lv].coordinates();

    const real*  lv_obj = field.meshValue(lv).valueObjLS().lv_obj();

    const Tree& tree = field.tree();
    const int*  mesh_offsets = tree.mesh_offsets();

    const Parameters& parameters = field.parameters();
    const real dx     = dx_[lv];
    const real c_ref  = parameters.c_ref_lbm();

    for (int l=0; l<(int)id_tasks.size(); l++) {
        for (int k=0; k<DefAMR::NX_LEAF; k++) {
        for (int j=0; j<DefAMR::NX_LEAF; j++) {
        for (int i=0; i<DefAMR::NX_LEAF; i++) {
            const int  idl = id_tasks[l];

            const int  offset = mesh_offsets[idl];
            const int  id     = Index::id0(i,j,k, offset);

            const real xcell = coordinates.xcell(id);
            const real ycell = coordinates.ycell(id);
            const real zcell = coordinates.zcell(id);
            if (zcell < z_min_nc3d_ || zcell >= z_max_nc3d_) { continue; }
            if ( !is_center( xcell, ycell ) )                { continue; }

            int  idz;
            real z_obj;
            get_idz_zobj_nc3d(idz, z_obj, zcell, dx, lv_obj[id]);

            add_rms_values(lv, idz, z_obj, u[id]*c_ref, v[id]*c_ref, w[id]*c_ref);
        }
        }
        }
    }

    average_rms_values(lv);
}


void  PostprocessNaturalConvection3d::init_value_cf()
{
    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        for (int i=0; i<numz_[lv]; i++) {
            values_[lv][i].is_exist  = false;
            values_[lv][i].is_existg = false;
            values_[lv][i].num       = 0;
            values_[lv][i].lv        = lv;

            values_[lv][i].z         = z_min_nc3d_ - 1.0;

            values_[lv][i].u         = 0.0;
            values_[lv][i].v         = 0.0;
            values_[lv][i].w         = 0.0;

            values_[lv][i].uu        = 0.0;
            values_[lv][i].vv        = 0.0;
            values_[lv][i].ww        = 0.0;

//            values_[lv][i].uv        = 0.0;
//            values_[lv][i].vw        = 0.0;
//            values_[lv][i].wu        = 0.0;
        }
    }
}


void  PostprocessNaturalConvection3d::add_exist_flag(const int lv, const int idz, const real z_obj)
{
    if ( values_[lv][idz].is_exist && values_[lv][idz].z != z_obj ) { std::cout << "error !! " << __PRETTY_FUNCTION__ << "values.z\n"; exit(0); }

    values_[lv][idz].is_exist  = true;
    values_[lv][idz].num      += 1;
    values_[lv][idz].z         = z_obj;
}


void  PostprocessNaturalConvection3d::add_mean_values(const int lv, const int idz, const real z_obj, const real u, const real v, const real w, const real T)
{
    if ( values_[lv][idz].is_exist ) {
        values_[lv][idz].u += u;
        values_[lv][idz].v += v;
        values_[lv][idz].w += w;

        values_[lv][idz].T += T;
    }
}


void  PostprocessNaturalConvection3d::add_rms_values(const int lv, const int idz, const real z_obj, const real u, const real v, const real w)
{
    if ( values_[lv][idz].is_exist ) {
        values_[lv][idz].uu += pow(u - values_[lv][idz].u, 2);
        values_[lv][idz].vv += pow(v - values_[lv][idz].v, 2);
        values_[lv][idz].ww += pow(w - values_[lv][idz].w, 2);

//        values_[lv][idz].uv += (u - values_[lv][idz].u) * (v - values_[lv][idz].v);
//        values_[lv][idz].vw += (v - values_[lv][idz].v) * (w - values_[lv][idz].w);
//        values_[lv][idz].wu += (w - values_[lv][idz].w) * (u - values_[lv][idz].u);
    }
}


void  PostprocessNaturalConvection3d::reduce_exist_flag(const int lv)
{
    int*   _is_exist  = new int[numz_[lv]];
    int*   _is_existg = new int[numz_[lv]];

    for (int i=0; i<numz_[lv]; i++) {
        _is_exist[i] = values_[lv][i].is_exist ? 1 : 0;
    }

    MPI_Allreduce(_is_exist, _is_existg, numz_[lv], MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    for (int i=0; i<numz_[lv]; i++) {
        values_[lv][i].is_existg = _is_existg[i];
    }

    delete [] _is_exist;
    delete [] _is_existg;

}


void  PostprocessNaturalConvection3d::reduce_z(const int lv)
{
    real*   _z  = new real[numz_[lv]];
    real*   _zg = new real[numz_[lv]];

    for (int i=0; i<numz_[lv]; i++) {
        _z[i] = values_[lv][i].z;
    }

    MPI_Allreduce(_z, _zg, numz_[lv], MFLOAT, MPI_MAX, MPI_COMM_WORLD);

    for (int i=0; i<numz_[lv]; i++) {
        values_[lv][i].z = _zg[i];
    }

    delete [] _z;
    delete [] _zg;

}


void  PostprocessNaturalConvection3d::average_num(const int lv)
{
    int*   _num  = new int[numz_[lv]];
    int*   _numg = new int[numz_[lv]];

    for (int i=0; i<numz_[lv]; i++) {
        if ( values_[lv][i].is_exist ) { _num[i] = values_[lv][i].num; }
        else                           { _num[i] = 0; }
    }

    MPI_Allreduce(_num, _numg, numz_[lv], MPI_INT, MPI_SUM, MPI_COMM_WORLD);


    for (int i=0; i<numz_[lv]; i++) {
        if ( values_[lv][i].is_existg ) {
            values_[lv][i].num = _numg[i];
        }
    }

    delete [] _num;
    delete [] _numg;

}


void  PostprocessNaturalConvection3d::average_mean_values(const int lv)
{
    real* _u = new real[numz_[lv]];
    real* _v = new real[numz_[lv]];
    real* _w = new real[numz_[lv]];
    real* _T = new real[numz_[lv]];

    real* _ug = new real[numz_[lv]];
    real* _vg = new real[numz_[lv]];
    real* _wg = new real[numz_[lv]];
    real* _Tg = new real[numz_[lv]];

    for (int i=0; i<numz_[lv]; i++) {
        if ( values_[lv][i].is_exist ) {
            _u[i] = values_[lv][i].u;
            _v[i] = values_[lv][i].v;
            _w[i] = values_[lv][i].w;
            _T[i] = values_[lv][i].T;
        }
        else {
            _u[i] = 0.0;
            _v[i] = 0.0;
            _w[i] = 0.0;
            _T[i] = 0.0;
        }
    }

    MPI_Allreduce(_u, _ug, numz_[lv], MFLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(_v, _vg, numz_[lv], MFLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(_w, _wg, numz_[lv], MFLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(_T, _Tg, numz_[lv], MFLOAT, MPI_SUM, MPI_COMM_WORLD);

    for (int i=0; i<numz_[lv]; i++) {
        if ( values_[lv][i].is_existg ) {
            values_[lv][i].u  = _ug[i] / values_[lv][i].num;
            values_[lv][i].v  = _vg[i] / values_[lv][i].num;
            values_[lv][i].w  = _wg[i] / values_[lv][i].num;
            values_[lv][i].T  = _Tg[i] / values_[lv][i].num;
        }
    }

    delete [] _u;
    delete [] _v;
    delete [] _w;
    delete [] _T;
    
    delete [] _ug;
    delete [] _vg;
    delete [] _wg;
    delete [] _Tg;
}


void  PostprocessNaturalConvection3d::average_rms_values(const int lv)
{
    real* _uu = new real[numz_[lv]];
    real* _vv = new real[numz_[lv]];
    real* _ww = new real[numz_[lv]];

//    real* _uv = new real[numz_[lv]];
//    real* _vw = new real[numz_[lv]];
//    real* _wu = new real[numz_[lv]];

    real* _uug = new real[numz_[lv]];
    real* _vvg = new real[numz_[lv]];
    real* _wwg = new real[numz_[lv]];

//    real* _uvg = new real[numz_[lv]];
//    real* _vwg = new real[numz_[lv]];
//    real* _wug = new real[numz_[lv]];
    for (int i=0; i<numz_[lv]; i++) {

        if ( values_[lv][i].is_exist ) {
            _uu[i] = values_[lv][i].uu;
            _vv[i] = values_[lv][i].vv;
            _ww[i] = values_[lv][i].ww;

//            _uv[i] = values_[lv][i].uv;
//            _vw[i] = values_[lv][i].vw;
//            _wu[i] = values_[lv][i].wu;
        }
        else {
            _uu[i] = 0.0;
            _vv[i] = 0.0;
            _ww[i] = 0.0;

//            _uv[i] = 0.0;
//            _vw[i] = 0.0;
//            _wu[i] = 0.0;
        }
    }

    MPI_Allreduce(_uu, _uug, numz_[lv], MFLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(_vv, _vvg, numz_[lv], MFLOAT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(_ww, _wwg, numz_[lv], MFLOAT, MPI_SUM, MPI_COMM_WORLD);

//    MPI_Allreduce(_uv, _uvg, numz_[lv], MFLOAT, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(_vw, _vwg, numz_[lv], MFLOAT, MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(_wu, _wug, numz_[lv], MFLOAT, MPI_SUM, MPI_COMM_WORLD);

    for (int i=0; i<numz_[lv]; i++) {
        if ( values_[lv][i].is_existg ) {
            const int num = values_[lv][i].num;

            values_[lv][i].uu  = _uug[i] / values_[lv][i].num;
            values_[lv][i].vv  = _vvg[i] / values_[lv][i].num;
            values_[lv][i].ww  = _wwg[i] / values_[lv][i].num;

//            values_[lv][i].uv  = _uvg[i] / values_[lv][i].num;
//            values_[lv][i].vw  = _vwg[i] / values_[lv][i].num;
//            values_[lv][i].wu  = _wug[i] / values_[lv][i].num;
        }
    }

    delete [] _uu;
    delete [] _vv;
    delete [] _ww;
    
//    delete [] _uv;
//    delete [] _vw;
//    delete [] _wu;
    
    delete [] _uug;
    delete [] _vvg;
    delete [] _wwg;

//    delete [] _uvg;
//    delete [] _vwg;
//    delete [] _wug;
}


void  PostprocessNaturalConvection3d::output_values(int step)
{
    if (rank_ != 0) { return; }

    // sort //
    std::vector<value_cf> values_statistic;

    for (int lv=0; lv<DefAMR::LV_MAX; lv++) {
        for (int i=0; i<numz_[lv]; i++) {
            if ( values_[lv][i].is_existg ) {
                values_statistic.emplace_back( values_[lv][i] );
            }
        }
    }

    std::sort( values_statistic.begin(), values_statistic.end(), [](value_cf a, value_cf b){ return (a.z < b.z); } );


    // fout //
    std::ofstream fout;
    fout.open(filename(step));

    fout << "z,lv,T,u,v,w,uu,vv,ww,num" << std::endl;
//    fout << "z,z+,lv,u,v,w,uu,vv,ww,uv,vw,wu,num" << std::endl;
    
    for (auto& elem : values_statistic) {
        fout << elem.z  << "," 
             << elem.lv << "," 
             << elem.T << "," 
             << elem.u  << "," << elem.v  << "," << elem.w << ","
             << elem.uu << "," << elem.vv << "," << elem.ww << ","
//             << elem.uv << "," << elem.vw << "," << elem.wu << ","
             << elem.num
             << std::endl;
    }

    fout.close();
}


void  PostprocessNaturalConvection3d::get_idz_zobj_nc3d(int& idz, real& z_obj, const real z, const real dx, const real lv_obj)
{
    const double z_floar     = z_min_nc3d_;
    const double tmp_nlv_obj = z - z_floar;

    idz   = floor(tmp_nlv_obj/dx);
    z_obj = (idz + 0.5)*dx;
}


bool PostprocessNaturalConvection3d::is_center(const real x, const real y)
{
    return (   x >= x_min_nc3d_ && x <= x_max_nc3d_ 
            && y >= y_min_nc3d_ && y <= y_max_nc3d_ ) ? true : false;
}
