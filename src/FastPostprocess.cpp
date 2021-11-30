#include "FastPostprocess.h"

#include "defineFilenames.h"
#include "defineLBM.h"

#include "Index.h"
#include "IndexLBM.h"

#include "range.hpp"
#include "foreach.h"

#include <string>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

void FastPostprocess::
readdata()
{
#ifdef NO_POSTPROCESS_MONITOR
if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skiped within -D NO_POSTPROCESS_MONITOR" << std::endl; }
return;
#endif

    heads_input.clear();
    x.clear();
    y.clear();
    z.clear();

    while(true) {
        // open input-*.csv
        const auto fname = Foldernames::input_folder + "/monitor/" + in_filename_ + 
            std::to_string(int(heads_input.size())) + ".csv";
        if(comm_.is_rank0()) { 
            std::cout << __PRETTY_FUNCTION__ << ": " << fname << std::endl; 
        }
        std::ifstream fin(fname, std::ios::in);
        if (!fin) { 
            if(comm_.is_rank0()) { 
                std::cout << __PRETTY_FUNCTION__ << ": finish scaning input-*.csv" << std::endl; 
                std::cout << " idx_total = " << int(x.size()) << std::endl;
            }
            break; 
        }

        heads_input.push_back(x.size());
        if(comm_.is_rank0()) { std::cout << "read index = " << int(heads_input.back()) << std::endl; }

        // read input-*.csv
        auto split = [](const std::string& input, char delimiter) {
            std::istringstream stream(input);
            std::string field;
            std::vector<std::string> result;
            while (std::getline(stream, field, delimiter)) {
                result.push_back(field);
            }
            return result;
        };
        std::string line;
        while(std::getline(fin, line)) {
            auto&& str_data = split(line, ','); // x, y, z
            x.push_back(std::stof(str_data.at(0)));
            y.push_back(std::stof(str_data.at(1)));
            z.push_back(std::stof(str_data.at(2)));
        }

        fin.close();
    } // while true; for input-*.csv

}

void FastPostprocess::
mallocdata(int n_scalars)
{
#ifdef NO_POSTPROCESS_MONITOR
if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skiped within -D NO_POSTPROCESS_MONITOR" << std::endl; }
return;
#endif

    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const auto nn = x.size();
    ids.resize(nn);
    calcflg.resize(nn);
    d_rank.resize(nn);
    u.resize(nn);
    v.resize(nn); 
    w.resize(nn);
    xitp.resize(nn);
    yitp.resize(nn);
    zitp.resize(nn);
    levelset_obj.resize(nn);
    scalars.resize(nn * n_scalars);
    T.resize(nn);

#ifdef USE_VALUE_STAT
    u_stat.resize(nn);
    v_stat.resize(nn); 
    w_stat.resize(nn);
    sc_stat.resize(nn * n_scalars);
    T_stat.resize(nn);
    uu_stat.resize(nn);
    vv_stat.resize(nn); 
    ww_stat.resize(nn);
    TT_stat.resize(nn);
#endif // USE_VALUE_STAT

    for(auto&& ii: ids) { 
        ii.i = ii.j = ii.k = ii.l = ii.lv = -1;
        ii.xx = ii.yy = ii.zz = NAN;
    }
    for(auto&& cc: calcflg) { cc = 0; }
    for(auto&& dd: d_rank) { dd = -1; }
}


// search cell ids from coordinates //
void FastPostprocess::
setup_coordinate_to_cell(
const Tree& tree,
const MeshValue* meshValues_ptr[]
) {
#ifdef NO_POSTPROCESS_MONITOR
if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skiped within -D NO_POSTPROCESS_MONITOR" << std::endl; }
return;
#endif

    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    
    const auto rank = comm_.world().rank();
    const auto n_leaf = tree.number_of_nodes();
    const auto nx_leaf = DefAMR::NX_LEAF;

    for(auto&& l: util::irange(n_leaf)) {
        const auto& node = tree.nodes().at(l);
        if(node.nodeCalFlags().Cal() == false ) { continue; }
        const int lv = node.level();

        const auto& meshValue = *(meshValues_ptr[lv]);
        const auto& coordinates = meshValue.coordinates();
        const auto& offsets = node.neighbor_mesh_offsets();
        const auto offset0 = offsets.offset0();

        // coordinates
        const auto dx_cell = coordinates.dx();

        const auto x0_leaf = coordinates.xnode(offset0);
        const auto y0_leaf = coordinates.ynode(offset0);
        const auto z0_leaf = coordinates.znode(offset0);
        const auto x1_leaf = coordinates.xnode(offset0, nx_leaf);
        const auto y1_leaf = coordinates.ynode(offset0, nx_leaf);
        const auto z1_leaf = coordinates.znode(offset0, nx_leaf);

        // search monitor in the leaf
        st_ids* ids = this->ids.data();
        char* calcflg = this->calcflg.data();
        int* d_rank = this->d_rank.data();
        float* x = this->x.data();
        float* y = this->y.data();
        float* z = this->z.data();
        const auto n_monitor = this->x.size();
        foreach::exec_1d<foreach::opti>(
            n_monitor,
            [=] __HD__ () {
                FOR_EACH1D_XX(im, n_monitor) {
                    // skip `out of range` or `already found`
                    if(im >= n_monitor) { SKIP_FOR(); }
                    if(calcflg[im] != 0) { SKIP_FOR(); }

                    // search leaf
                    const auto xc = x[im], 
                               yc = y[im], 
                               zc = z[im];
                    if(! (      x0_leaf -epsilon*dx_cell <= xc && xc < epsilon*dx_cell + x1_leaf
                             && y0_leaf -epsilon*dx_cell <= yc && yc < epsilon*dx_cell + y1_leaf
                             && z0_leaf -epsilon*dx_cell <= zc && zc < epsilon*dx_cell + z1_leaf
                    )) { SKIP_FOR(); }

                    // determine nearest cell; since cell nearest, maybe -0.5 < xx < 0.5
                    const real ii = (xc - x0_leaf) / dx_cell,
                               jj = (yc - y0_leaf) / dx_cell,
                               kk = (zc - z0_leaf) / dx_cell;
                    int i = int(ii), j = int(jj), k = int(kk);
                    /*sanity*/ if(i<0||i>=nx_leaf||j<0||j>=nx_leaf||k<0||k>=nx_leaf) { SKIP_FOR(); }
                    const auto xx = ii - i - 0.5, yy = jj - j - 0.5, zz = kk - k - 0.5; // 0.5 == cell-center offset

                    // store result
                    ids[im].i = i;
                    ids[im].j = j;
                    ids[im].k = k;
                    ids[im].l = l;
                    ids[im].lv = lv;
                    ids[im].xx = xx;
                    ids[im].yy = yy;
                    ids[im].zz = zz;
                    calcflg[im] = 1;
                    d_rank[im] = rank;
                } // [=] __HD__ ()
        }); // foreach
    }

    // check if all monitor points are found. if not, throw runtime_error
    int count = 0;
    for(auto&& im: util::irange(monitorcount())) {
        if(calcflg[im] > 0) { count++; }
    }

    // to print which station includes missing monitor points
    for(auto&& fileno: util::irange(heads_input.size())) {
        const auto n_monitor = this->x.size();
        const auto im_start = heads_input.at(fileno);
        const auto im_end = (fileno+1 == heads_input.size()) ? n_monitor : heads_input.at(fileno+1);
        int count_in_station = 0;
         
        for(auto&& im: util::irange(im_start, im_end)) {
            if(calcflg[im] > 0) count_in_station++;
        }
        auto count_sum = comm_.col_vector().reduce_sum(count_in_station);
        auto expected = im_end - im_start;
        if(comm_.is_rank0() && (count_sum < expected) ) {
            std::cout << __PRETTY_FUNCTION__ << ": station not found @ fileno = " << fileno << ", count_sum: " << count_sum <<     ", expected: " << expected << std::endl;
        }
    }

    auto count_sum = comm_.col_vector().reduce_sum(count);
    if(comm_.is_rank0()) { 
        std::cout << __PRETTY_FUNCTION__ << ": duplicative search finished. " << std::endl;
        std::cout << __PRETTY_FUNCTION__ << ": n monitor tmp total = " << count_sum << "/" << monitorcount() << std::endl;
    }
    #ifdef FASTPOST_IGNORE_MISSING
    if(count_sum < monitorcount() && comm_.is_rank0()) {
        std::cerr << __PRETTY_FUNCTION__ << ": warning: some monitor point(s) were not found. continue anyway." << std::endl;
    }
    #else
    runtime_assert(count_sum >= monitorcount(), "ValueError: some monitor point(s) were not found.");
    #endif

    // drop duplicates
    const auto mpirank = comm_.col_vector().rank();
    const auto mpisize = comm_.col_vector().size();
    std::vector<int> owner(monitorcount()), owner_g(monitorcount());
    for(auto&& im: util::irange(monitorcount())) {
        owner[im] = (calcflg[im] > 0) ? mpirank : mpisize;
    }
    MPI_Allreduce(owner.data(), owner_g.data(), monitorcount(), MPI_INT, MPI_MIN, comm_.col_vector().comm());
    for(auto&& im: util::irange(monitorcount())) {
        if(calcflg[im] > 0 && (owner_g[im] != owner[im])) {
            calcflg[im] = 0;
        }
    }

    // check if all monitor points are found && duplicates are droped
    count = 0;
    for(auto&& im: util::irange(monitorcount())) {
        if(calcflg[im] > 0) { count++; }
    }
    count_sum = comm_.col_vector().reduce_sum(count);
    if(comm_.is_rank0()) { 
        std::cout << __PRETTY_FUNCTION__ << ": drop_duplicates finished. " << std::endl;
        std::cout << __PRETTY_FUNCTION__ << ": n monitor = " << count_sum << "/" << monitorcount() << std::endl;
    }
    #ifdef FASTPOST_IGNORE_MISSING
    if(count_sum != monitorcount() && comm_.is_rank0()) {
        std::cerr << __PRETTY_FUNCTION__ << ": warning: drop_duplicates got invalid result. continue anyway." << std::endl;
    }
    #else
    runtime_assert(count_sum == monitorcount(), "InternalError: drop_duplicates got invalid result");
    #endif
}

void FastPostprocess::
setupdata(
const Tree& tree,
const MeshValue* meshValues_ptr[]
) {
#ifdef NO_POSTPROCESS_MONITOR
if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skiped within -D NO_POSTPROCESS_MONITOR" << std::endl; }
return;
#endif

    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }

    const auto w_comm = comm_.world().comm();
    const auto w_rank = comm_.world().rank();
    const auto w_size = comm_.world().size();

    // read files on each
    for(auto&& rank: util::irange(w_size)) {
        if(rank == w_rank) {
            readdata();
            if(comm_.is_rank0()) { std::cout << " monitorcount = " << monitorcount() << std::endl; }
        }
        MPI_Barrier(w_comm);
    }

    const int n_scalars = meshValues_ptr[0]->n_scalars();
    mallocdata(n_scalars);
    setup_coordinate_to_cell(tree, meshValues_ptr);

    // output dir
    ::mkdir(out_prefix0_().c_str(), 0775);
    ::mkdir(out_prefix_ens_().c_str(), 0775);
}


void FastPostprocess::
OutputMonitorData(
int step,
const Tree& tree,
const Parameters& parameters, 
const MeshValue* meshValues_ptr[]
#ifdef USE_VALUE_STAT
, const ValueStat* valueStat
#endif // USE_VALUE_STAT
)
{
#ifdef NO_POSTPROCESS_MONITOR
if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << ": skiped within -D NO_POSTPROCESS_MONITOR" << std::endl; }
return;
#endif
    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << std::endl; }
    if(comm_.is_rank0()) { std::cout << " monitor count: " << monitorcount() << std::endl; }

    // scan monitor data
    const auto n_monitor = this->x.size();
    const int n_scalars = meshValues_ptr[0]->n_scalars();
    const st_ids* ids = this->ids.data();
    const char* calcflg = this->calcflg.data();
    const int* d_rank = this->d_rank.data();
    float* u = this->u.data();
    float* v = this->v.data();
    float* w = this->w.data();
    float* xitp = this->xitp.data();
    float* yitp = this->yitp.data();
    float* zitp = this->zitp.data();
    float* levelset_obj = this->levelset_obj.data();
    float* scalar = this->scalars.data();
    float* T = this->T.data();

#ifdef USE_VALUE_STAT
    float* u_stat  = this->u_stat.data();
    float* v_stat  = this->v_stat.data();
    float* w_stat  = this->w_stat.data();
    float* sc_stat = this->sc_stat.data();
    float* T_stat  = this->T_stat.data();
    float* uu_stat = this->uu_stat.data();
    float* vv_stat = this->vv_stat.data();
    float* ww_stat = this->ww_stat.data();
    float* TT_stat = this->TT_stat.data();
#endif // USE_VALUE_STAT

    const float* x = this->x.data();
    const float* y = this->y.data();
    const float* z = this->z.data();
    const auto* mesh_offsets3x3x3 = tree.mesh_offsets3x3x3();
    for(auto&& lv: util::irange(DefAMR::LV_MAX)) {
        const auto& meshValue = *(meshValues_ptr[lv]);
        const auto& coordinates = meshValue.coordinates();
        const auto dx_cell = coordinates.dx();
        const int   nn_max    = meshValue.nn_max();
        const auto* u_mesh = meshValue.valueNS().u();
        const auto* v_mesh = meshValue.valueNS().v();
        const auto* w_mesh = meshValue.valueNS().w();
        const auto* lv_obj_mesh = meshValue.valueObjLS().lv_obj();
        const auto* scalar_mesh = meshValue.valueNS().scalar();
        const auto* T_mesh = meshValue.valueNS().T();

#ifdef USE_VALUE_STAT
        const auto  stat_count = valueStat[lv].t_count();
        const auto  inv_stat_count = 1.f/stat_count;

        const auto* u_sum  = valueStat[lv].u_sum();
        const auto* v_sum  = valueStat[lv].v_sum();
        const auto* w_sum  = valueStat[lv].w_sum();
        const auto* sc_sum = valueStat[lv].sc_sum();
        const auto* T_sum  = valueStat[lv].T_sum();
        const auto* uu_sum = valueStat[lv].uu_sum();
        const auto* vv_sum = valueStat[lv].vv_sum();
        const auto* ww_sum = valueStat[lv].ww_sum();
        const auto* TT_sum = valueStat[lv].TT_sum();
#endif // USE_VALUE_STAT

        const auto c_ref = parameters.c_ref_lbm();

        foreach::exec_1d<foreach::opti>( 
            n_monitor,
            [=] __HD__ () {
                FOR_EACH1D_XX(im, n_monitor) {
                    if(calcflg[im] == 0) { return; } // skip noncal leaf

                    const auto id = ids[im];
                    if(id.lv != lv) { return; } // kernel called for each level

                    // mesh id
                    int offset3d[27];
                    for(int idv=0; idv<27; idv++) {
                        offset3d[idv] = mesh_offsets3x3x3[27*id.l + idv];
                    }

                    // interpolate val
                    float _u=0, _v=0, _w=0, _xitp=0, _yitp=0, _zitp=0, _l=0, _T=0;
#ifdef USE_VALUE_STAT
                    float _u_sum=0, _v_sum=0, _w_sum=0, _T_sum=0;
                    float _uu_sum=0, _vv_sum=0, _ww_sum=0, _TT_sum=0;
#endif // USE_VALUE_STAT
                    for(int kku=0; kku<2; kku++) {
                    for(int jju=0; jju<2; jju++) {
                    for(int iiu=0; iiu<2; iiu++) { // iiu == 0: this cell, iiu == 1: neighbor cell
                        const auto ii = (id.xx < 0 ? -iiu : iiu),
                                   jj = (id.yy < 0 ? -jju : jju),
                                   kk = (id.zz < 0 ? -kku : kku); // direction of neighbor cell
                        const auto axx = abs(id.xx), ayy = abs(id.yy), azz = abs(id.zz);  // abs of sub_cell offset
                        const auto wx = axx*iiu + ((iiu+1)%2)*(1.f-axx),
                                   wy = ayy*jju + ((jju+1)%2)*(1.f-ayy),
                                   wz = azz*kku + ((kku+1)%2)*(1.f-azz); // axx: weight for neighbor cell, (1-axx): weight for this cell
                        const auto weight = wx*wy*wz;

                        const auto id_mesh = Index::id(id.i+ii, id.j+jj, id.k+kk, offset3d);
                        _u += weight * u_mesh[id_mesh];
                        _v += weight * v_mesh[id_mesh];
                        _w += weight * w_mesh[id_mesh];
                        _xitp += weight * coordinates.xcell(id_mesh);
                        _yitp += weight * coordinates.ycell(id_mesh);
                        _zitp += weight * coordinates.zcell(id_mesh);
                        _l += weight * lv_obj_mesh[id_mesh];
                        _T += weight * T_mesh[id_mesh];

#ifdef USE_VALUE_STAT
                        _u_sum  += weight * u_sum[id_mesh];
                        _v_sum  += weight * v_sum[id_mesh];
                        _w_sum  += weight * w_sum[id_mesh];
                        _T_sum  += weight * T_sum[id_mesh];
                        _uu_sum += weight * uu_sum[id_mesh];
                        _vv_sum += weight * vv_sum[id_mesh];
                        _ww_sum += weight * ww_sum[id_mesh];
                        _TT_sum += weight * TT_sum[id_mesh];
#endif // USE_VALUE_STAT
                    }}}

                    // scalar vars separated
                    for(int n=0; n<n_scalars; n++) {
                        float _s=0;
#ifdef USE_VALUE_STAT
                        float _sc_sum=0;
#endif // USE_VALUE_STAT
                        for(int kk=0; kk<2; kk++) {
                        for(int jj=0; jj<2; jj++) {
                        for(int ii=0; ii<2; ii++) {
                            const auto wx = id.xx*ii + ((ii+1)%2)*(1.f-id.xx);
                            const auto wy = id.yy*jj + ((jj+1)%2)*(1.f-id.yy);
                            const auto wz = id.zz*kk + ((kk+1)%2)*(1.f-id.zz);
                            const auto weight = wx*wy*wz;

                            const auto id_mesh = nn_max * n + Index::id(id.i+ii, id.j+jj, id.k+kk, offset3d); // to compute 2D index
                            _s += weight * scalar_mesh[id_mesh];

#ifdef USE_VALUE_STAT
                            _sc_sum += weight * sc_sum[id_mesh];
#endif // USE_VALUE_STAT
                        }}}
                        scalar[n_monitor * n + im] = _s;
#ifdef USE_VALUE_STAT
                        sc_stat[n_monitor * n + im] = (inv_stat_count) * _sc_sum;
#endif // USE_VALUE_STAT
                    }

                    // write final
                    u[im] = _u * c_ref;
                    v[im] = _v * c_ref;
                    w[im] = _w * c_ref;
                    xitp[im] = _xitp;
                    yitp[im] = _yitp;
                    zitp[im] = _zitp;
                    levelset_obj[im] = _l;
                    T[im] = _T;

#ifdef USE_VALUE_STAT
                    u_stat[im]  = (inv_stat_count) * _u_sum * c_ref;
                    v_stat[im]  = (inv_stat_count) * _v_sum * c_ref;
                    w_stat[im]  = (inv_stat_count) * _w_sum * c_ref;
                    T_stat[im]  = (inv_stat_count) * _T_sum;
                    uu_stat[im] = (inv_stat_count) * _uu_sum * c_ref * c_ref;
                    vv_stat[im] = (inv_stat_count) * _vv_sum * c_ref * c_ref;
                    ww_stat[im] = (inv_stat_count) * _ww_sum * c_ref * c_ref;
                    TT_stat[im] = (inv_stat_count) * _TT_sum;
#endif // USE_VALUE_STAT
                } // FOR_EACH1D_XX
            } // [=] __HD__ ()
        ); // foeach::exec_1d
    } // for lv

    /// fprintf
    for(auto&& fileno: util::irange(heads_input.size())) {
        const auto filename = out_filename_(-1, fileno);
        if(step == 0) { // make new
            std::ofstream fout(filename, std::ios::trunc);
            std::stringstream ss;
            for(int i=0; i<n_scalars; i++) {
                ss << "scalar" + std::to_string(i) << ',';
            } 
#ifdef USE_VALUE_STAT
            std::stringstream ss2;
            for(int i=0; i<n_scalars; i++) {
                ss2 << "scalar_stat" + std::to_string(i) << ',';
            } 
#endif // USE_VALUE_STAT

            runtime_assert(fout, "FileIOError");
            fout << "step" << ','
                 << "calcflg" << ','
                 << "d_rank"  << ','
                 << "u" << ',' << "v" << ',' << "w" << ','
                 << "xitp" << ',' << "yitp" << ',' << "zitp" << ','
                 << "levelset_obj" << ','
                 << ss.str()
                 << "T" << ','
#ifdef USE_VALUE_STAT
                 << "u_stat" << ',' << "v_stat" << ',' << "w_stat" << ','
                 << ss2.str()
                 << "T_stat" << ','
                 << "uu_stat" << ',' << "vv_stat" << ',' << "ww_stat" << ','
                 << "TT_stat" << ','
#endif // USE_VALUE_STAT
                 << "x" << ',' << "y" << ',' << "z" << std::endl;
            fout.close();
        }
        std::ofstream fout(filename, std::ios::app | std::ios::ate);
        runtime_assert(fout, "FileIOError");
        const auto im_start = heads_input.at(fileno);
        const auto im_end = (fileno+1 == heads_input.size()) ? n_monitor : heads_input.at(fileno+1);
        for(auto&& im: util::irange(im_start, im_end)) {
            if(!calcflg[im]) { continue; } // will be on other rank
            std::stringstream ss;
            for(int i=0; i<n_scalars; i++) {
                ss << scalar[n_monitor * i + im] << ',';
            } 
#ifdef USE_VALUE_STAT
            std::stringstream ss2;
            for(int i=0; i<n_scalars; i++) {
                ss2 << sc_stat[n_monitor * i + im] << ',';
            } 
#endif // USE_VALUE_STAT
            fout << step << ','
                << int(calcflg[im]) << ','
                << d_rank[im] << ','
                << u[im] << ','
                << v[im] << ','
                << w[im] << ','
                << xitp[im] << ','
                << yitp[im] << ','
                << zitp[im] << ','
                << levelset_obj[im] << ','
                << ss.str()
                << T[im] << ','
#ifdef USE_VALUE_STAT
                << u_stat[im] << ','
                << v_stat[im] << ','
                << w_stat[im] << ','
                << ss2.str()
                << T_stat[im] << ','
                << uu_stat[im] << ','
                << vv_stat[im] << ','
                << ww_stat[im] << ','
                << TT_stat[im] << ','
#endif // USE_VALUE_STAT
                << x[im] << ','
                << y[im] << ','
                << z[im] << std::endl;
        }

        fout.close();
    }
    if(comm_.is_rank0()) { std::cout << __PRETTY_FUNCTION__ << "finished." << std::endl; }
}
